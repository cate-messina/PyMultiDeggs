#!/usr/bin/env python
# coding: utf-8


import os
import re
import json
import warnings

from datetime import datetime
from io import StringIO

import pandas as pd
import requests
from importlib.resources import files

from bs4 import BeautifulSoup


def _default_network_path():
    """Get path to the bundled OmniPath network data file.
    
    Returns the file system path for the default human protein-protein interaction
    network included with the package. This network is used when no user cache exists
    and is always available offline.
    
    Returns:
        Path object pointing to the bundled defaultNetowrk_human.csv file.
    """
    return files("multideggs.data").joinpath("defaultNetowrk_human.csv")

def _user_cache_paths(organism: int = 9606, archive_version: str = None):
    """Construct paths for user-cached network data and metadata.
    
    Generates file paths for storing user-downloaded network data (CSV) and its
    associated metadata (JSON). Creates the cache directory (~/.cache/multideggs) if
    needed. Cached networks take priority over bundled networks in load_network().
    
    Args:
        organism (int, optional): NCBI taxonomy ID for the organism network.
            Defaults to 9606 (human).
        
        archive_version (str, optional): Specific OmniPath archive version (format: "YYYYMMDD").
            If provided, used in the filename to distinguish frozen builds.
            Defaults to None.
    
    Returns:
        tuple: (csv_path, meta_path) where:
            - csv_path (str): Path to cached network CSV file; includes archive_version suffix if specified
            - meta_path (str): Path to metadata JSON file with download information
    """
    home_dir = os.path.expanduser("~")
    cache_dir = os.path.join(home_dir, ".cache", "multideggs")
    
    suffix = f"_{archive_version}" if archive_version else ""
    
    # Create cache directory if it doesn't exist
    os.makedirs(cache_dir, exist_ok=True)
    
    # Define file paths
    csv_path  = os.path.join(cache_dir, f"defaultNetowrk_human_{organism}_user{suffix}.csv")
    meta_path = os.path.join(cache_dir, f"defaultNetowrk_human_{organism}_user{suffix}_meta.json")
    
    return csv_path, meta_path

def _check_for_updates(organism: int = 9606) -> tuple[bool, str]:
    
    """Check if OmniPath has released a newer database build.
    
    Parses the OmniPath archive index page to find the most recent
    interactions file, identified by the pattern YYYYMMDD-YYYYMMDD in
    the filename. The second date (the "to" date) represents the actual
    biological content cutoff and is compared against the locally cached
    version to determine if an update is available.
    
    Args:
        organism (int, optional): NCBI organism ID. Defaults to 9606 (human).
    
    Returns:
        tuple: (is_newer, remote_version) where:
            - is_newer (bool): True if remote version is newer than local cache
            - remote_version (str): Remote build date as YYYYMMDD, "No local cache",
              or "Unable to verify server status" if the check fails
    """
    
    csv_path, meta_path = _user_cache_paths(organism=organism, archive_version=None)
    
    if not os.path.exists(meta_path):
        return True, "No local cache"
    
    with open(meta_path) as f:
        local_meta = json.load(f)
    local_version = local_meta.get("remote_version", "")
    
    try:
        response = requests.get("https://archive.omnipathdb.org/", timeout=10)
        soup = BeautifulSoup(response.text, "html.parser")
        
        pattern = re.compile(r"omnipath_webservice_interactions__(\d{8})-(\d{8})\.tsv\.xz")
        
        dates_to = []
        for link in soup.find_all("a", href=True):
            m = pattern.search(link["href"])
            if m:
                dates_to.append(m.group(2))  # data "to"
        
        if not dates_to:
            return False, "Unable to verify server status"
        
        remote_version = max(dates_to)  # es. "20250813"
        is_newer = remote_version > local_version
        return is_newer, remote_version
        
    except Exception:
        return False, "Unable to verify server status"

def _fetch_omnipath_raw(organism: int = 9606, archive_version: str = None) -> pd.DataFrame:
    """Fetch raw gene-gene or protein-protein interaction data from OmniPath database.
    
    Retrieves human (organism 9606) interactions with reference and curation effort information, using gene symbols as identifiers.
    
    Args:
        organism (int, optional): NCBI taxonomy ID for interaction orthology.
            Defaults to 9606 (human). Note: Archive versions contain human data only.
        
        archive_version (str, optional): Specific OmniPath archive build (format: "YYYYMMDD").
            If None, fetches latest live data from web service. Defaults to None.
    
    Returns:
        pd.DataFrame: Raw interaction data with columns:
            - source_genesymbol: Source gene symbol
            - target_genesymbol: Target gene symbol
            - references: Semicolon-separated PubMed IDs
            - curation_effort: Curation score per interaction source
    
    Raises:
        RuntimeError: If connection to OmniPath fails or archive version is unavailable.
    """
    if archive_version:
        if organism != 9606:
            # Warn user that archive is human-only; live service supports orthologs
            warnings.warn(
                f"[MultiDEGGs] OmniPath archive files contain human data only (organism 9606). "
                f"Orthology-translated data for organism {organism} is only available via the live web service. "
                f"To download mouse/rat data, use update_network(organism={organism}) without archive_version."
            )
        print(f"[MultiDEGGs] Check the package README for the correct archive_version format (expected: 'YYYYMMDD-YYYYMMDD').")
        # Construct URL to frozen archive version, if the user want a specific version
        url = f"https://archive.omnipathdb.org/omnipath_webservice_interactions__{archive_version}.tsv.xz"
    else:
        # Construct URL to live web service with species and data source filters
        url = (
        f"https://omnipathdb.org/interactions"
        f"?organisms={organism}"
        f"&fields=references,curation_effort,type"  
        f"&genesymbols=1"
        f"&datasets=omnipath,kinaseextra,ligrecextra,pathwayextra,"
        f"mirnatarget,dorothea,collectri,tf_target,lncrna_mrna,tf_mirna"
    )
    
    try:
        # Fetch data with timeout to prevent hanging
        response = requests.get(url, timeout=90)
        response.raise_for_status()
    
    except Exception as e:
        # Provide helpful error messages based on request type
        if archive_version:
            raise RuntimeError(
                f"Unable to download archive version {archive_version}: {e}. "
                f"Check that it exists at https://archive.omnipathdb.org/"
            )
        raise RuntimeError(f"Connection error to OmniPath web service: {e}")
    
    # Parse TSV response and return as DataFrame
    return pd.read_csv(StringIO(response.text), sep="\t")


def _build_network_df(raw_df: pd.DataFrame) -> pd.DataFrame:
    """Build and filter the interaction network from raw OmniPath data.
    
    Constructs a clean network with columns: from, to, n_references and curation_effort. Removes entries with
    invalid gene symbols (containing underscores), and aggregates reference counts
    per gene pair.
    
    Args:
        raw_df (pd.DataFrame): Raw OmniPath interaction data with columns:
            - source_genesymbol: Source gene symbol
            - target_genesymbol: Target gene symbol
            - references: Semicolon-separated PubMed IDs (may be null)
            - curation_effort: Curation effort score (may be null)
    
    Returns:
        pd.DataFrame: Processed network with deduplicated edges and columns:
            - from: Source gene symbol
            - to: Target gene symbol
            - n_references: Count of PubMed references per edge
            - curation_effort: Aggregated curation score per edge
    """
    df = pd.DataFrame({
        "from": raw_df["source_genesymbol"],
        "to": raw_df["target_genesymbol"],
        "references": raw_df["references"].fillna(""),
        "curation_effort": raw_df["curation_effort"].fillna(0).astype(int),
        "is_stimulation": raw_df["is_stimulation"].fillna(0).astype(int),
        "is_inhibition": raw_df["is_inhibition"].fillna(0).astype(int),
        "interaction_type": raw_df["type"].fillna(""),
        "is_directed": raw_df["is_directed"].fillna(0).astype(int),
        "consensus_stimulation": raw_df["consensus_stimulation"].fillna(0).astype(int),
        "consensus_inhibition": raw_df["consensus_inhibition"].fillna(0).astype(int),
        "consensus_direction": raw_df["consensus_direction"].fillna(0).astype(int),
    })
    df["PMDI_interaction"] = raw_df["references"].fillna("").apply(
    lambda x: ";".join(sorted(set(idx.strip() for idx in str(x).split(";") if idx.strip())))
    )

    # Remove rows with missing gene symbols
    df = df.dropna(subset=["from", "to"])
    
    # Filter out genes with underscores in names 
    mask = (
        ~df["from"].str.contains("_", na=False) &
        ~df["to"].str.contains("_", na=False)
    )
    df = df[mask].copy()

    # Count PubMed references per edge (extract 4+ digit PubMed IDs)
    df["n_references"] = df["references"].apply(
        lambda x: len(re.findall(r"\d{4,}", str(x)))
    )

    # Aggregate duplicate gene pairs (sum reference counts and curation effort)
    df = df.groupby(["from", "to"], as_index=False).agg({
    "n_references":        "sum",
    "curation_effort":     "sum",
    "interaction_type":    lambda x: ";".join(sorted(set(x))),
    "is_stimulation":      "max",
    "is_inhibition":       "max",
    "is_directed":         "max",
    "consensus_stimulation": "max",    
    "consensus_inhibition":  "max",
    "consensus_direction":   "max",
    "PMDI_interaction":    lambda x: ";".join(sorted(set(";".join(x).split(";")) - {""})),
})
        
    return df[["from", "to", "n_references", "curation_effort",
           "interaction_type",
           "is_stimulation", "consensus_stimulation",
           "is_inhibition",  "consensus_inhibition",
           "is_directed",    "consensus_direction",
           "PMDI_interaction"]].reset_index(drop=True)    


def update_network(organism: int = 9606, archive_version: str = None, verbose: bool = True) -> pd.DataFrame:
    """Download and cache the latest OmniPath network data.
    
    Fetches fresh protein-protein interaction data from OmniPath database and stores
    it in the user's cache directory (~/.cache/multideggs/) with metadata. Creates a
    reproducible frozen build if archive_version is specified.
    
    Args:
        organism (int, optional): NCBI taxonomy ID. Defaults to 9606 (human).
        
        archive_version (str, optional): Specific OmniPath archive build (format: "YYYYMMDD").
        If None, fetches current live data. When set, creates a frozen reproducible build
        saved to a separate cache file to avoid overwriting the live cache.
        Defaults to None.
        
        verbose (bool, optional): If True, prints status messages about download progress.
            Defaults to True.
    
    Returns:
        pd.DataFrame: Processed network with columns from, to, n_references, curation_effort.
    
    Raises:
        RuntimeError: If network download fails due to connectivity or invalid archive version.
    """
    if verbose:
        if archive_version:
            print(f"[MultiDEGGs] Downloading frozen network build {archive_version} from OmniPath Archive...")
        else:
            print(f"[MultiDEGGs] Downloading latest live network from OmniPath...")

    # Fetch and process raw data
    raw = _fetch_omnipath_raw(organism=organism, archive_version=archive_version)
    df  = _build_network_df(raw)

    # Save processed network to cache
    csv_path, meta_path = _user_cache_paths(organism=organism, archive_version=archive_version)
    df.to_csv(csv_path, index=False)
    
    # Determine the build ID to write in metadata
    if archive_version:
        current_remote_version = archive_version
    else:
        # Try to retrieve current server version, fallback to download date
        _, current_remote_version = _check_for_updates(organism=organism)
        if current_remote_version in ["No local cache", "Unable to verify server status"]:
            current_remote_version = datetime.now().strftime("%Y%m%d")
    
    # Write metadata sidecar with download information
    with open(meta_path, "w") as f:
        json.dump({
            "downloaded_at": datetime.now().isoformat(),
            "organism":      organism,
            "n_edges":       len(df),
            "remote_version": current_remote_version,
            "is_frozen":     archive_version is not None,
        }, f, indent=2)

    if verbose:
        print(f"[MultiDEGGs] Network saved: {len(df):,} edges  →  {csv_path}")

    return df


def load_network(organism: int = 9606, archive_version: str = None, verbose: bool = True) -> pd.DataFrame:
    """Load the OmniPath network, preferring cached data when available.
    
    Loads the best available network following this priority order:
    1. User cache (created by update_network()) if present
    2. Bundled network shipped with package (always offline-available)
    
    Note: User-supplied networks passed via the network parameter in get_diffNetworks()
    bypass this function entirely.
    
    Args:
        organism (int, optional): NCBI taxonomy ID. Defaults to 9606 (human).
        
        archive_version (str, optional): Specific frozen build (format: "YYYYMMDD").
            If not cached locally, raises FileNotFoundError. Defaults to None.
        
        verbose (bool, optional): If True, prints information about which network
            version is being used and alerts about available updates. Defaults to True.
    
    Returns:
        pd.DataFrame: Network with columns from, to, n_references, curation_effort.
    
    Raises:
        FileNotFoundError: If archive_version specified but not found in cache,
            or if non-human organism requested without cache available.
    """
    csv_path, meta_path = _user_cache_paths(organism=organism)

    # Try to use user cache if available
    if os.path.exists(csv_path):
        meta = {}
        
        # Load metadata if available
        if os.path.exists(meta_path):
            with open(meta_path) as f:
                meta = json.load(f)
        
        # Extract metadata for status reporting
        local_version = meta.get("remote_version", "unknown")
        n_edges = meta.get("n_edges", "?")
        is_frozen = meta.get("is_frozen", False)
        
        if verbose:
            # Format version date (YYYYMMDD -> YYYY-MM-DD) for readability
            if len(local_version) == 8:
                local_v_fmt = f"{local_version[:4]}-{local_version[4:6]}-{local_version[6:]}"
            else:
                local_v_fmt = local_version

            if is_frozen:
                # Reproducible build from specific archive version
                print(f"[MultiDEGGs] Using FROZEN user cache (Reproducible Build: {local_v_fmt}, edges: {n_edges:,})")
            else:
                # Current live data version
                print(f"[MultiDEGGs] Using user cache (OmniPath Build: {local_v_fmt}, edges: {n_edges:,})")

        # Check for updates if not a frozen build
        if not is_frozen:
            has_update, remote_version = _check_for_updates(organism=organism)
            if verbose:
                if has_update and remote_version not in ["No local cache", "Unable to verify server status"]:
                    # Format remote version for readability
                    remote_v_fmt = f"{remote_version[:4]}-{remote_version[4:6]}-{remote_version[6:]}"
                    print(f"\n[NETWORK UPDATE AVAILABLE]")
                    print(f"  • Local cache build: {local_v_fmt}")
                    print(f"  • Latest OmniPath build: {remote_v_fmt}")
                    print(f"  • Run multideggs.update_network() to download the latest version.\n")
                elif remote_version == "Unable to verify server status":
                    print("[MultiDEGGs] Could not check for updates (offline or server unavailable).")

        # Load and return cached network
        return pd.read_csv(csv_path)

    # Fallback to bundled network for human organism with no archive version specified
    if organism == 9606 and archive_version is None:
        if verbose:
            print(f"[MultiDEGGs] No cache found for organism {organism}.")
            print("[MultiDEGGs] Using bundled human network (offline).")
        ref = _default_network_path()
        with ref.open("r") as f:
            return pd.read_csv(f)
    else:
        # If a frozen version or specific organism is missing, the user must first run update_network
        ver_str = f" build {archive_version}" if archive_version else ""
        raise FileNotFoundError(
            f"No cached data found for organism {organism}{ver_str}. "
            f"Please run update_network(organism={organism}, archive_version={repr(archive_version)}) first while online."
        )
    