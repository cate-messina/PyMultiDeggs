#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np

from IPython.display import display, HTML

# Import functions from the core
from .core_functions import get_diffNetworks, get_sig_degg, get_multiOmics_diffNetworks
from .final_table import MultiDEGGsTable, MultiDEGGsTableDict, show_html_results

def run_multideggs(
    assayData: pd.DataFrame | np.ndarray | list | dict,
    metadata: pd.DataFrame,
    category_variable: str,
    category_subset: list = None,
    padj_method: str = "bonferroni",
    sig_threshold: float = 0.05,
    percentile_vector: list | np.ndarray = None,
    regression_method: str = "lm",
    verbose: bool = True,
    organism: int = 9606,
    archive_version: str = None
):
    """Run the complete MultiDEGGs pipeline for one or more omics layers.

    This wrapper orchestrates metadata validation, network loading, differential
    network inference, and significance filtering in a single call. It supports
    both single-omics and multi-omics inputs.

    Args:
        assayData (pd.DataFrame | np.ndarray | list | dict): Omics expression data. This can be a single matrix or multiple datasets.
        
        metadata (pd.DataFrame): Sample annotations indexed by sample IDs.
        
        category_variable (str): Column name in `metadata` containing the sample grouping variable.
        
        category_subset (list, optional): Subset of categories to retain.
            Defaults to None.
            
        padj_method (str, optional): Multiple-testing adjustment method.
            Defaults to "bonferroni".
            
        sig_threshold (float, optional): Significance cutoff for retained edges.
            Defaults to 0.05.
            
        percentile_vector (list | np.ndarray, optional): Percentile thresholds to evaluate during percolation analysis. 
            Defaults to None.
            
        regression_method (str, optional): Regression method for edge testing.
            Defaults to "lm".
            
        verbose (bool, optional): Whether to print workflow progress and warnings.
            Defaults to True.
            
        organism (int, optional): Taxonomic ID of the interaction network.
            Defaults to 9606.
            
        archive_version (str, optional): OmniPath archive version to use. 
            If `None`, the cached or bundled network is used.

    Returns:
        tuple: A `(results, final_outputs)` tuple where: 
        - `results` contains the full differential-network analysis 
        - `final_outputs` contains the significant edges in a display-friendly table.

    Raises:
        ValueError: If required metadata or sample group information is missing or
            inconsistent.
        KeyError: If the grouping variable is absent from the metadata.
        FileNotFoundError: If the selected network cannot be loaded.
    """
    if verbose:
        print("\n" + "="*60)
        print("[Pipeline] Starting complete MultiDEGGs workflow...")
        print("="*60)

    from .omnipath_network import load_network
    network = load_network(organism=organism, archive_version=archive_version, verbose=verbose)

    results = get_diffNetworks(
        assayData=assayData,
        metadata=metadata,
        category_variable=category_variable,
        category_subset=category_subset,
        padj_method=padj_method,
        percentile_vector=percentile_vector,
        regression_method=regression_method,
        verbose=verbose,
        network=network 
    )
    
    # Determine analysis pathway: detect if single-omic or multi-omic dataset by counting the number of individual assay layers in the results
    internal_assays = results["assayData"]
    
    if len(internal_assays) == 1:
        # Single-omic pathway: extract significant edges from one assay layer
        assay_name = list(internal_assays.keys())[0]
        if verbose:
            print(f"\n[Pipeline] Detected single-omic dataset ('{assay_name}').")
            print(f"[Pipeline] Extracting significant edges (p-adj threshold: {sig_threshold})...")
            
        final_outputs = MultiDEGGsTable(get_sig_degg(
            degg=results, 
            assayDataName=assay_name, 
            sig_threshold=sig_threshold,
            network=network
        ))
    else:
        # Multi-omic pathway: merge and aggregate significant edges across assay layers
        if verbose:
            print(f"\n[Pipeline] Detected multi-omic dataset ({len(internal_assays)} layers).")
            print(f"[Pipeline] Merging layers and extracting significant edges (p-adj threshold: {sig_threshold})...")
            
        final_outputs = MultiDEGGsTableDict(get_multiOmics_diffNetworks(
            degg=results, 
            sig_threshold=sig_threshold,
            network=network
        ))
        
    if verbose:
        print("="*60)
        print("[Pipeline] Workflow completed successfully!")
        print("="*60 + "\n")
    
    # Display results: render DataFrames as formatted tables (HTML when available)
    print("\n[Pipeline] Significant edges found:")
    display(final_outputs)
  
    return results, final_outputs