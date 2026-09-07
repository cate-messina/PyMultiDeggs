#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from itertools import combinations
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests
from importlib.resources import files
from joblib import Parallel, delayed 
import warnings 

# Load OmniPath network
from .omnipath_network import load_network


# In[ ]:

def _make_gene_url(gene):
    return f"https://www.ncbi.nlm.nih.gov/gene/?term={gene}"


def tidy_metadata(metadata: pd.DataFrame,
                  category_variable: str,  
                  category_subset: list = None,  
                  verbose: bool = True):
    """Validate and filter sample metadata for differential network analysis.
    
    This function prepares metadata for the differential analysis by optionally subsetting to
    specified categories, removing samples with missing category values, filtering out
    categories with insufficient sample size (< 5 samples), and removing unused categorical
    levels to ensure data integrity.

    Args:
        metadata (pd.DataFrame): Sample metadata indexed by sample IDs, with annotation columns.
        
        category_variable (str): Column name in metadata containing category/group labels
            for sample grouping.
        
        category_subset (list, optional): Specific category values to retain in metadata.
            If None, all categories are retained after quality filtering. Defaults to None.
        
        verbose (bool, optional): If True, prints warnings about removed samples and groups.
            Defaults to True.

    Returns:
        pd.DataFrame: Cleaned metadata containing only valid samples and categories with
            sufficient sample sizes (≥5 samples) for differential analysis.

    Raises:
        ValueError: If fewer than 2 categories remain with at least 5 samples each after filtering.
    """
    # If a subset of categories is provided, keep only those rows
    if category_subset is not None:
        metadata = metadata[metadata[category_variable].isin(category_subset)].copy()

        # If the column is categorical, remove unused categories after subsetting
        if isinstance(metadata[category_variable].dtype, pd.CategoricalDtype):
            metadata[category_variable] = metadata[category_variable].cat.remove_unused_categories()

    # Remove samples with NA in the category column and warn if any were removed
    na_count = metadata[category_variable].isna().sum()
    metadata = metadata[metadata[category_variable].notna()].copy()
    if na_count > 0 and verbose:
        warnings.warn(f"Category column contains {na_count} NaN values; samples were removed.")

    category_vector = metadata[category_variable]

    # Count number of samples in each category
    counts = category_vector.value_counts()

    # Identify categories with insufficient sample size for statistical analysis
    small_groups = counts[counts < 5].index.tolist()

    if small_groups:
        warnings.warn(f"Removing categories with <5 samples: {small_groups}")
        # Filter out samples belonging to categories with too few samples
        metadata = metadata[~metadata[category_variable].isin(small_groups)].copy()

        # If the column is categorical, remove unused categories after removal
        if isinstance(metadata[category_variable].dtype, pd.CategoricalDtype):
            metadata[category_variable] = metadata[category_variable].cat.remove_unused_categories()

        category_vector = metadata[category_variable]

    # Verify that at least two groups remain for differential analysis
    if category_vector.nunique() < 2:
        raise ValueError(
            f"To perform differential analysis at least 2 categories are needed "
            f"with 5+ samples. Found: {category_vector.nunique()}"
        )

    return metadata


# In[ ]:

def get_diffNetworks_singleOmic(
    assayData: pd.DataFrame,
    assayDataName: str,
    category_variable: str,
    metadata: pd.DataFrame,
    padj_method: str = "bonferroni",
    regression_method: str = "lm",
    percentile_vector = None,
    verbose = True,
    remove_shared_edges = True,
    network = None,
    organism: int = 9606):
    """Generate differential networks for a single omics layer.

    This function builds category-specific biological networks, applies a gene
    expression percentile filter, estimates edge significance with regression
    models, and selects the percolation threshold that maximizes the number of
    significant interactions.

    Args:
        assayData (pd.DataFrame): Features (rows) x Samples (columns) matrix.
        
        assayDataName (str): Name of the assay for logging purposes.
        
        category_variable (str): Column name in metadata defining sample groups.
        
        metadata (pd.DataFrame): Sample metadata indexed by sample IDs.
        
        padj_method (str, optional): P-value adjustment method. 
            Defaults to "bonferroni".
        
        regression_method (str, optional): Regression method for differential analysis. 
            Defaults to "lm".
        
        percentile_vector (array-like, optional): Percentiles to test in percolation analysis. 
            Defaults to np.arange(0.35, 0.981, 0.05).
        
        verbose (bool, optional): Enable warning messages. 
            Defaults to True.
        
        remove_shared_edges (bool, optional): Remove edges shared across category contrasts. 
            Defaults to True.
        
        network (pd.DataFrame, optional): Pre-loaded biological network with 'from' and 'to' columns representing gene interactions. If None, loads the OmniPath network (cached or bundled). 
            Defaults to None.
        
        organism (int, optional): NCBI organism ID for OmniPath network.
            Defaults to 9606 (human).
    
    Returns:
        dict: Dictionary containing:
            - Category names as keys mapping to edge DataFrames with 'p.value' and 'p.adj' columns
            - "_best_percentile": Optimal percolation threshold (0-1) that maximizes significant edges
            - "_sig_pvalues_count": Number of significant edges (p < 0.05) at best threshold
            - "_all_percentiles_results": Complete results for all tested percentile thresholds
    
    Raises:
        ValueError: If assayData is not a DataFrame, sample IDs don't match,
            no differential analysis is possible, or network validation fails.
    """
    
    if verbose:
        print(f"Processing {assayDataName} differential network...")

    # Validate input data format
    if not isinstance(assayData, pd.DataFrame):
        raise ValueError(f"{assayDataName} must be a pandas DataFrame.")
    
    # Identify samples present in both metadata and assayData
    common_samples = assayData.columns.intersection(metadata.index)
    
    if len(common_samples) == 0:
        raise ValueError(f"No matching sample IDs between metadata and {assayDataName}.")

    # Keep only samples present in metadata
    assayData = assayData[common_samples]
    
    # Filter out samples with excessive missing values (< 2 NA values)
    assayData = assayData.loc[:, assayData.notna().sum() >= 2]

    # Identify and warn about samples in metadata but missing from assayData
    missing_metadata = metadata.index.difference(assayData.columns)
    if len(missing_metadata) != 0 and verbose:
        warnings.warn(f"The following samples IDs are missing in {assayDataName}:\n{', '.join(missing_metadata)}")

    # Re-align metadata to match filtered assayData samples
    metadata = metadata.loc[assayData.columns]
    if metadata[category_variable].nunique() == 1:
        raise ValueError(f"{assayDataName}: All samples belong to one category; differential analysis requires ≥2 groups.")

    # Select biological network
    if network is not None:
        network_to_use = network          
    else:
        network_to_use = load_network(verbose=verbose, organism=organism)  # cached or bundled
        
    # Filter edges: 'from' and 'to' must be present in assayData rownames
    edges = network_to_use[
        network_to_use['from'].isin(assayData.index) & 
        network_to_use['to'].isin(assayData.index)].copy()

    if edges.empty:
        raise ValueError(f"Rownames of {assayDataName} don't match any link of the biological network.")

    # Remove duplicate gene pair edges
    edges = edges.drop_duplicates(subset=['from', 'to'])

    # Filter assayData to keep only nodes present in the network
    nodes = pd.unique(edges[['from', 'to']].values.ravel())
    assayData = assayData.loc[assayData.index.isin(nodes)]

    # Identify all unique sample categories/groups
    categories = metadata[category_variable].unique()
    
    # Generate all pairwise category comparisons
    contrasts = list(combinations(categories, 2))
    
    # Pre-compute median gene expression per category for percentile filtering
    category_median_dict = {}
    for cat in categories:
        sample_ids = metadata[metadata[category_variable] == cat].index
        category_median_dict[cat] = assayData[sample_ids].median(axis=1, skipna=True)

    if percentile_vector is None:
        percentile_vector = np.arange(0.35, 0.981, 0.05)  # 0.35 to 0.98 percentiles

    all_results = percolation_analysis(
        assayData_clean=assayData,
        category_median_dict=category_median_dict,
        edges=edges,
        contrasts=contrasts,
        metadata=metadata,
        category_variable=category_variable,
        padj_method=padj_method,
        regression_method=regression_method,
        remove_shared_edges=remove_shared_edges,
        percentile_vector=percentile_vector
    )

    # Identify optimal percolation threshold that maximizes significant edges
    best = max(all_results, key=lambda x: x["n_significant_arch"])

    if verbose:
        warnings.warn(
            f"Percolation analysis: genes below the "
            f"{best['percentile'] * 100:.0f}th percentile removed from networks."
        )
    
    output_dict = best["edges_pvalue"].copy()
    
    # Append analysis metadata (prefixed with '_' to distinguish from sample categories)
    output_dict["_best_percentile"] = best["percentile"]
    output_dict["_sig_pvalues_count"] = best["n_significant_arch"]
    output_dict["_all_percentiles_results"] = all_results

    return output_dict


# In[ ]:


def calc_pvalue_network(edge_df,
                        assayData_clean,
                        metadata,
                        category_variable,
                        padj_method,
                        regression_method = "lm"
                        ):
    """Calculate p-values for a set of gene-gene edges.

    For two groups, the function fits a linear model with a gene-by-group
    interaction term. For three or more groups, it uses ANOVA on the interaction
    term. The resulting p-values are optionally adjusted with a multiple-testing
    correction.

    Args:
        edge_df (pd.DataFrame): DataFrame containing `from` and `to` columns for
            gene pairs.
            
        assayData_clean (pd.DataFrame): Expression matrix with genes as rows and
            samples as columns.
            
        metadata (pd.DataFrame): Sample metadata indexed by sample IDs.
        
        category_variable (str): Column name in metadata that defines sample groups.
        
        padj_method (str): P-value adjustment method (e.g.,'bonferroni','fdr_bh').
            Use "none" to skip adjustment.
        
        regression_method (str, optional): Regression method to use. Supported values are "lm" and "rlm". 
            Defaults to "lm"
    
    Returns:
        pd.DataFrame: Input edge_df with added columns:
            - "p.value": Interaction p-value for each edge
            - "p.adj": Adjusted p-value (if padj_method != "none")
            - "edge_id": Unique edge identifier (set as index)
        
        Returns string message if edge_df is empty or not a DataFrame.
        
        Edges where regression failed (e.g. zero-variance genes) receive NaN in "p.value" and are excluded from multiple testing correction for "p.adj".
    
    Raises:
        NotImplementedError: If `regression_method` is neither "lm" nor "rlm".
    """
    
    if not isinstance(edge_df, pd.DataFrame) or edge_df.empty:
        return "No valid edges to test in this category."
    
    if regression_method not in ("lm", "rlm"):
        raise NotImplementedError(
            f"regression_method='{regression_method}' not supported. Use 'lm' or 'rlm'."
        )
    
    n_categories = metadata[category_variable].nunique()
    
    pvalue_list = []
    
    # Ensure metadata rows align with assayData columns before regression loop
    metadata = metadata.loc[assayData_clean.columns]
    
    for geneA, geneB in edge_df[["from", "to"]].values:
        # Construct data frame for regression with gene pair and category information
        df_linregress = pd.DataFrame({
            "gene_A": assayData_clean.loc[geneA].values,
            "gene_B": assayData_clean.loc[geneB].values,
            "category": metadata[category_variable].values})

        if regression_method == "rlm":
            model = smf.rlm("gene_B ~ gene_A * category", data=df_linregress).fit()
            
            mask = (
                model.pvalues.index.str.contains("gene_A", regex=False) &
                model.pvalues.index.str.contains("category", regex=False)
            )
            interaction_terms = model.pvalues[mask]
            interaction_pval = interaction_terms.iloc[0] if not interaction_terms.empty else np.nan

        elif regression_method == "lm":
            # Standard linear regression with interaction term
            model = smf.ols("gene_B ~ gene_A * category", data=df_linregress).fit()

            if n_categories == 2:
                # 2-group case: extract p-value from interaction coefficient
                interaction_terms = model.pvalues.filter(like="gene_A:category")
                interaction_pval = interaction_terms.iloc[0] if not interaction_terms.empty else np.nan

            elif n_categories >= 3:
                # Multi-group case: use ANOVA F-test for interaction effect
                anova_table = anova_lm(model, typ=1)
                interaction_row = anova_table.filter(like="gene_A:category", axis=0)
                interaction_pval = interaction_row["PR(>F)"].iloc[0] if not interaction_row.empty else np.nan

        pvalue_list.append(interaction_pval)   

    edge_df = edge_df.copy()
    edge_df["p.value"] = pvalue_list       

    if padj_method != "none":
        pvalue_array = np.array(pvalue_list, dtype=float)
        valid_mask = ~np.isnan(pvalue_array)
        edge_df["p.adj"] = np.nan
        if valid_mask.any():
            _, p_adj_valid, _, _ = multipletests(pvalue_array[valid_mask], method=padj_method)
            edge_df.loc[valid_mask, "p.adj"] = p_adj_valid

    # Create unique edge identifier and set as index
    edge_df["edge_id"] = edge_df["from"] + "-" + edge_df["to"]
    edge_df = edge_df.set_index("edge_id")

    return edge_df
        
                    


# In[ ]:


def calc_pvalue_percentile(assayData_clean,
                           p,
                           category_median_dict,
                           edges,
                           contrasts,
                           remove_shared_edges = True):
    """Filter edges using a gene-expression percentile threshold.

    For each category, genes with median expression below the threshold are
    removed, and only edges connecting retained genes are kept. Optional shared-
    edge removal prevents duplicates across category contrasts.

    Args:
        assayData_clean (pd.DataFrame): Expression matrix with genes as rows and
            samples as columns.
            
        p (float): Percentile threshold in the [0, 1] interval.
        
        category_median_dict (dict): Mapping between category names and gene-wise
            median expression values.
            
        edges (pd.DataFrame): Network edge table containing `from` and `to`.
        
        contrasts (list): Category pairs to compare.
        
        remove_shared_edges (bool, optional): Removes edges shared between contrasted categories. 
            Defaults to True.

    Returns:
        dict: Category-to-edge-DataFrame mapping after percentile filtering.
    """
    
    # Compute percentile-based expression threshold across all genes and samples
    threshold_p = np.nanquantile(assayData_clean.to_numpy(), p)
    
    # Identify genes with median expression above percentile threshold per category
    keep_genes = {}
    for category, medians in category_median_dict.items():
        keep_genes[category] = []
        for gene, median in medians.items():
            # Include genes with median above percentile threshold
            if median > threshold_p:
                keep_genes[category].append(gene)
    
    # Filter edges: retain only edges where both nodes are above threshold
    best_edge_results = {}
    for category, genes_above in keep_genes.items():
        # Keep edges where both genes meet the percentile threshold
        filtered = edges[
            edges['from'].isin(genes_above) &
            edges['to'].isin(genes_above)].copy()
        best_edge_results[category] = filtered
    
    # Optionally remove edges shared across category contrasts
    if remove_shared_edges:
        # Identify edges that are common across contrasted category pairs
        common_edges = set()
        for cat1, cat2 in contrasts:
            set1 = set(zip(best_edge_results[cat1]['from'], best_edge_results[cat1]['to']))
            set2 = set(zip(best_edge_results[cat2]['from'], best_edge_results[cat2]['to']))
            common_edges |= set1 & set2
        
        # Remove shared edges from all categories
        common_keys = {f"{a} {b}" for a, b in common_edges}
        for category, edge_df in best_edge_results.items():
            edge_keys = edge_df['from'] + " " + edge_df['to']
            best_edge_results[category] = edge_df[~edge_keys.isin(common_keys)].copy()
    
    return best_edge_results
                      


# In[ ]:


def single_percentile_calc(assayData_clean,
                           p,
                           category_median_dict,
                           edges,
                           contrasts,
                           metadata,
                           category_variable,
                           padj_method,
                           regression_method="lm",
                           remove_shared_edges=True):
    """Calculate p-values and count significant edges for a single percentile threshold.
    
    The function filters the network by the requested expression threshold,computes interaction p-values for each retained edge, and counts how many edges remain significant after multiple-testing correction.
    
    Args:
        assayData_clean (pd.DataFrame): Expression matrix with genes as rows and samples as columns.
        
        p (float): Percentile threshold used to retain highly expressed genes.
        
        category_median_dict (dict): Per-category median expression values.
        
        edges (pd.DataFrame): Network edge table containing `from` and `to`.
        
        contrasts (list): Category pairs used for the pairwise comparison.
        
        metadata (pd.DataFrame): Sample metadata indexed by sample IDs.
        
        category_variable (str): Grouping column in `metadata`.
        
        padj_method (str): Multiple-testing correction method.
        
        regression_method (str, optional): Regression method for the interaction test. 
            Defaults to "lm".
            
        remove_shared_edges (bool, optional): Whether to remove shared edges across category contrasts. 
            Defaults to True.

    Returns:
        dict: Dictionary containing:
            - "percentile": The percentile threshold used (rounded to 2 decimals)
            - "n_significant_arch": Count of significant edges (p < 0.05) across all categories
            - "edges_pvalue": Mapping of category names to DataFrames with p-values and adjusted p-values
    """
    
    # Filter edges based on gene expression percentile threshold by category
    filtered_edge = calc_pvalue_percentile(assayData_clean,
                                           p,
                                           category_median_dict,
                                           edges,
                                           contrasts,
                                           remove_shared_edges)
    
    # Compute p-values for statistical interaction of filtered edges
    pvalue_edge = {}
    for category, edge_df in filtered_edge.items():
        if edge_df is None or edge_df.empty:
            pvalue_edge[category] = pd.DataFrame()
        else:
            pvalue_edge[category] = calc_pvalue_network(
                edge_df, assayData_clean, metadata, category_variable, padj_method, regression_method)
    
    # Determine which p-value column to use for significance filtering
    if padj_method != "none":
        sig_col = "p.adj"
    else:
        sig_col = "p.value"

    # Count significant edges across all categories
    significant_arch = 0
    for edge_df in pvalue_edge.values():
        if isinstance(edge_df, pd.DataFrame) and sig_col in edge_df.columns:
            p_values = pd.to_numeric(edge_df[sig_col], errors='coerce')
            significant_arch += (p_values < 0.05).sum()
            
    # Return results for this percentile
    return {
        "percentile": round(p, 2),
        "n_significant_arch": significant_arch,
        "edges_pvalue": pvalue_edge
    }
    


# In[ ]:


def percolation_analysis(assayData_clean,
                         category_median_dict,
                         edges,
                         contrasts,
                         metadata,
                         category_variable,
                         padj_method,
                         regression_method="lm",
                         remove_shared_edges=True,
                         percentile_vector=None):
    """Perform percolation analysis across multiple percentile thresholds.
    
    This function iteratively evaluates several gene-expression percentile filters in parallel
    and returns one result object per threshold. These results can be compared to
    select the threshold that maximizes the number of significant edges
    
    Args:
        assayData_clean (pd.DataFrame): Expression matrix with genes as rows and samples as columns.
        
        category_median_dict (dict): Per-category median expression values.
        
        edges (pd.DataFrame): Network edge table containing `from` and `to`.
        
        contrasts (list): Category pairs used for the contrast analysis.
        
        metadata (pd.DataFrame): Sample metadata indexed by sample IDs.
        
        category_variable (str): Grouping column in `metadata`.
        
        padj_method (str): Multiple-testing adjustment method.
        
        regression_method (str, optional): Regression method used for edge testing.
            Defaults to "lm".
            
        remove_shared_edges (bool, optional): Whether to drop shared edges across contrasts. 
            Defaults to True.
            
        percentile_vector (array-like, optional): Thresholds to test. If `None`, the default range is used.

    Returns:
        list: A list of result dictionaries, one for each percentile threshold.
    """
    
    # Prepare percentile threshold values for parallel processing
    range_percentile = percentile_vector 
    
    # Compute p-values in parallel across percentile thresholds (n_jobs=-2: all but one CPU)
    all_results = Parallel(n_jobs=-2, verbose=10)(
        delayed(single_percentile_calc)(
            assayData_clean, p, category_median_dict, edges, 
            contrasts, metadata, category_variable, padj_method, 
            regression_method, remove_shared_edges
        ) 
        for p in range_percentile
    )

    return all_results
    


# In[ ]:


def get_diffNetworks(
    assayData: pd.DataFrame | np.ndarray | list | dict,
    metadata: pd.DataFrame,
    category_variable: str,
    category_subset: list = None,
    padj_method: str = "bonferroni",
    regression_method: str = "lm",
    percentile_vector: list | np.ndarray = None,
    verbose: bool = True,
    network = None,
    organism: int = 9606
):
    """Prepare and validate multi-omics data and metadata for differential network analysis.
    
    This function normalizes assay data into a standardized dictionary format, validates
    metadata structure, cleans sample annotations, and aligns sample IDs across all data layers.
    It serves as the primary entry point for preparing data before differential network generation.
    
    Args:
        assayData (pd.DataFrame | np.ndarray | list | dict): Omics data in DataFrame, array, list, or dictionary format.
        
        metadata (pd.DataFrame): Sample metadata indexed by sample IDs.
        
        category_variable (str): Name of the column defining the biological groups.
        
        category_subset (list, optional): Restriction to selected category values.
            If `None`, all valid categories are retained.
            
        padj_method (str, optional): Multiple-testing adjustment method. 
            Defaults to "bonferroni".
            
        regression_method (str, optional): Regression method used for differential testing. 
            Defaults to "lm".
            
        percentile_vector (list | np.ndarray, optional): Percentile thresholds to test. 
            If `None`, the default scan range is used.
            
        verbose (bool, optional): Whether to print or emit warnings during preprocessing. 
            Defaults to True.
            
        network (pd.DataFrame, optional): Preloaded interaction network with `from` and `to` columns.
            
        organism (int, optional): NCBI organism ID used to load the OmniPath network. 
            Defaults to 9606.

    Returns:
        dict: Dictionary containing:
            - "assayData": Normalized dict mapping dataset names to DataFrames (features x samples)
            - "metadata": Cleaned metadata aligned with assayData sample IDs
            - "padj_method": The p-value adjustment method specified
            - "category_variable": Column name of the grouping variable
            - "category_subset": Subset of categories used (if any)
            - "diffNetworks": Differential networks for each assay dataset with results for all categories
    
    Raises:
        ValueError: If the assay format is invalid or the metadata does not
            contain valid samples for the given categories.
        TypeError: If the metadata object is not DataFrame-compatible.
    """
    
    # Convert input data to standardized dictionary format (name -> DataFrame)
    if isinstance(assayData, (pd.DataFrame, np.ndarray)):
        assayData_dic = {"assayData1": pd.DataFrame(assayData)}
    elif isinstance(assayData, list):
        assayData_dic = {f"assayData{i+1}": pd.DataFrame(v) for i, v in enumerate(assayData)}
    elif isinstance(assayData, dict):
        assayData_dic = {k: pd.DataFrame(v) for k, v in assayData.items()}
    else:
        raise ValueError("assayData must be a DataFrame, ndarray, list, or dict.")

    # Convert metadata to DataFrame if necessary
    if isinstance(metadata, (np.ndarray, list)):
        metadata = pd.DataFrame(metadata)
    elif not isinstance(metadata, pd.DataFrame):
        raise TypeError("metadata must be a pd.DataFrame, np.ndarray, or list.")

    # Validate category_variable specification and presence
    if category_variable is None:
        raise ValueError("category_variable must be specified.")
    if category_variable not in metadata.columns:
        raise ValueError(f"The category variable column '{category_variable}' is missing in metadata.")

    # Clean metadata: delegate to tidy_metadata for sample and category filtering
    metadata_clean = tidy_metadata(
        metadata=metadata,
        category_variable=category_variable,
        category_subset=category_subset,
        verbose=verbose
    )

    # Align all assay datasets with cleaned metadata (remove unmapped samples)
    for k, v in list(assayData_dic.items()):
        common_ids = v.columns.intersection(metadata_clean.index)
        if len(common_ids) == 0:
            raise ValueError(f"No matching sample IDs between metadata and '{k}'.")
        missing = set(v.columns) - set(metadata_clean.index)
        if missing and verbose:
            warnings.warn(f"{k}: {len(missing)} samples removed (missing metadata).")
        # Retain only samples with valid metadata
        assayData_dic[k] = v[common_ids]

    if network is not None:
        if not isinstance(network, pd.DataFrame):
            raise TypeError("network must be a pd.DataFrame.")
        if not {"from", "to"}.issubset(network.columns):
            raise ValueError("network must have columns 'from' and 'to'.")

    # Process each assay dataset to generate differential networks
    diffNetworks = {}
    for name, data in assayData_dic.items():
        diffNetworks[name] = get_diffNetworks_singleOmic(
            assayData=data,
            assayDataName=name,
            category_variable=category_variable,
            metadata=metadata_clean,
            padj_method=padj_method,
            regression_method=regression_method,
            percentile_vector=percentile_vector,
            verbose=verbose,
            remove_shared_edges=True,
            network=network,
            organism=organism,
        )

    
    return {
        "assayData": assayData_dic,
        "metadata": metadata_clean,
        "padj_method": padj_method,
        "category_variable": category_variable,
        "category_subset": category_subset,
        "diffNetworks": diffNetworks
    }


# In[ ]:


def get_sig_degg(degg,
                 assayDataName=None,
                 sig_threshold=0.05,
                 network=None):
    """Extract significant edges from the differential network results.

    The function retains all edges whose adjusted or unadjusted p-value is below
    the selected threshold and aggregates the output across sample categories.

    Args:
        degg (dict): Differential network object returned by `get_diffNetworks`.
        
        assayDataName (str, optional): Name of the assay to extract. 
            If `None`, the first available assay is used.
        
        sig_threshold (float, optional): Significance threshold for retention.
            Defaults to 0.05.
        
        network (pd.DataFrame, optional): Reference interaction network used to add metadata columns such as interaction type and PMDI annotations.

    Returns:
        pd.DataFrame: A concatenated table of significant edges.

    Raises:
        ValueError: If `degg` is not a dictionary or does not contain the required `diffNetworks` key.
    """
    
    if not isinstance(degg, dict):
        raise ValueError("Degg is not a dictionary")
    
    if "diffNetworks" not in degg or not degg["diffNetworks"]:
        raise ValueError("'diffNetworks' key missing or empty in degg object.")
    
    # Use first assay if none specified
    if assayDataName is None:
        assayDataName = list(degg["diffNetworks"].keys())[0]
    
    # Select significance column based on adjustment method
    if degg["padj_method"] == "none":
        sig_var = "p.value"
    else:
        sig_var = "p.adj"
        
    networks = degg["diffNetworks"][assayDataName]
    
    sig_list = []
    
    # Extract and aggregate significant edges from each category
    for category, edge_df in networks.items():
        if not isinstance(edge_df, pd.DataFrame) or edge_df.empty:
            continue
        if sig_var not in edge_df.columns:
            warnings.warn(f"{sig_var} column missing for category '{category}'.")
            continue
        
        # Filter edges below significance threshold
        sig_edges = edge_df[edge_df[sig_var] < sig_threshold].copy()
        
        if sig_edges.empty:
            continue
        
        # Add category column for tracking
        sig_edges["category"] = category
        sig_list.append(sig_edges)
    
    # Return empty DataFrame if no significant edges found at any category
    if not sig_list:
        warnings.warn(f"No significant interactions found with {sig_var} < {sig_threshold}.")
        return pd.DataFrame()
    
    # Combine all categories into a single DataFrame
    sig_deggs = pd.concat(sig_list, ignore_index=True)
    
    if network is not None and "n_references" in network.columns:
        if "n_references" not in sig_deggs.columns:
            net_sub = network[["from", "to", "n_references", "curation_effort",
                               "interaction_type", "is_stimulation", "consensus_stimulation",
                               "is_inhibition", "consensus_inhibition",
                               "is_directed", "consensus_direction", "PMDI_interaction"]]
            sig_deggs = sig_deggs.merge(net_sub, on=["from", "to"], how="left")

            # Check both directions for reciprocal gene pairs
            missing = sig_deggs["n_references"].isna()
            if missing.any():
                net_rev = net_sub.rename(columns={"from": "to", "to": "from"})
                filled = sig_deggs[missing].drop(columns=["n_references", "curation_effort",
                                   "interaction_type", "is_stimulation", "consensus_stimulation",
                                   "is_inhibition", "consensus_inhibition",
                                   "is_directed", "consensus_direction", "PMDI_interaction"]).merge(
                    net_rev, on=["from", "to"], how="left"
                )
                for col in ["n_references", "curation_effort", "interaction_type",
                            "is_stimulation", "consensus_stimulation", "is_inhibition",
                            "consensus_inhibition", "is_directed", "consensus_direction",
                            "PMDI_interaction"]:
                    sig_deggs.loc[missing, col] = filled[col].values

    # Add NCBI Gene URLs for interactive exploration
    sig_deggs["url_from"] = sig_deggs["from"].apply(_make_gene_url)
    sig_deggs["url_to"]   = sig_deggs["to"].apply(_make_gene_url)

    return sig_deggs
            
        


# In[ ]:


def get_multiOmics_diffNetworks(degg,
                                sig_threshold = 0.05,
                                network=None):
    """Aggregate differential networks across multiple omics layers by sample category.
    
    This function combines network results from multiple omics datasets (multi-layer networks)
    for each sample category, filters edges by significance threshold, and organizes them
    with layer annotations.
    
    Args:
        degg (dict): Differential networks object containing:
            - "assayData": Dictionary of omics datasets
            - "diffNetworks": Dictionary mapping omics names to their network results
            - "padj_method": P-value adjustment method used
            - "metadata": Sample metadata
            - "category_variable": Column name defining sample groups
        
        sig_threshold (float, optional): Significance threshold for p-values. Edges with p-values below this threshold are retained. 
            Defaults to 0.05.

        network (pd.DataFrame, optional): Reference network used to enrich the
                    output with interaction metadata when available.
    
    Returns:
        dict: Dictionary mapping category names (keys) to aggregated edge DataFrames with columns:
            - All columns from original edge DataFrames (e.g., "from", "to", "p.value", "p.adj");
                row index is reset (integer, not edge_id) to avoid duplicate index conflicts across layers.
            - "layer": Omics dataset name for each edge.
            - Empty DataFrame for categories with no significant edges across any layer.
    
    Raises:
        ValueError: If degg is not a dictionary, lacks required keys, or contains only a single omics dataset.
    """
    
    if not isinstance(degg, dict):
        raise ValueError("Input must be a deggs object")
    
    if len(degg["assayData"]) == 1:
        raise ValueError(
            "Only 1 omics layer detected; use get_sig_degg() for single-omic analysis. "
            "This function is for multi-omic (≥2 layers) scenarios only."
        )
    
    # Extract unique sample categories/groups
    categories = degg["metadata"][degg["category_variable"]].unique().tolist()
    
    # Retrieve names of all omics datasets in analysis
    omicDatasets = list(degg["assayData"].keys())
    
    # Determine which p-value column to use for significance thresholding
    if degg["padj_method"] == "none":
        sig_var = "p.value"
    else:
        sig_var = "p.adj"
    
    multilayer_networks = {}
    
    for category in categories:
        category_networks = []
        
        for omicDataset in omicDatasets:
            # Fetch network edges for this omics layer and sample category
            edge_network = degg["diffNetworks"].get(omicDataset, {}).get(category)
            
            if not isinstance(edge_network, pd.DataFrame) or edge_network.empty:
                warnings.warn(f"{category} in {omicDataset} is empty or not valid")
                continue
            
            # Add layer annotation and filter by significance threshold
            edge_network = edge_network.copy()
            edge_network["layer"] = omicDataset
            # Retain only edges passing significance threshold
            edge_network = edge_network[edge_network[sig_var] < sig_threshold]
            category_networks.append(edge_network)
        
        # Aggregate significant edges from all omics layers for this category
        if category_networks:
            merged = pd.concat(category_networks, ignore_index=True)
            # Convert layer column to categorical for efficient storage
            merged["layer"] = merged["layer"].astype("category")
            multilayer_networks[category] = merged
        else:
            warnings.warn(f"No significant edges for category: {category}")
            multilayer_networks[category] = pd.DataFrame()

    for category, df in multilayer_networks.items():
        if not df.empty:
            out = df
            if network is not None and "n_references" in network.columns:
                if "n_references" not in out.columns:
                    out = out.merge(
                        network[["from", "to", "interaction_type", "n_references", "curation_effort",
                                "is_stimulation", "consensus_stimulation",
                                "is_inhibition",  "consensus_inhibition",
                                "is_directed",    "consensus_direction",
                                "PMDI_interaction"]],
                        on=["from", "to"], how="left"
)
            out = out.copy()
            out["url_from"] = out["from"].apply(_make_gene_url)
            out["url_to"]   = out["to"].apply(_make_gene_url)
            multilayer_networks[category] = out
    return multilayer_networks
    
             

