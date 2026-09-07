#!/usr/bin/env python
# coding: utf-8

import base64
import io
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pyvis.network import Network


# Colour-blind friendly palette (Okabe-Ito)
LAYER_COLOURS    = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
CATEGORY_COLOURS = ["#2166AC", "#D6604D", "#4DAC26", "#7B2D8B"]
NODE_COLOUR = "#AED6F1"
NODE_BORDER = "#1A5276"


def _layer_colour_map(layers):
    """Create a color mapping for omics layers to distinguishable hex codes.
    
    Maps layer names to colors from a color-blind-friendly palette (Okabe-Ito),
    cycling through colors if more layers exist than palette options.
    
    Args:
        layers (iterable): Layer names to map to colours.
    
    Returns:
        dict: Mapping of layer name (str) to color hex code (str).
    """
    return {layer: LAYER_COLOURS[i % len(LAYER_COLOURS)]
            for i, layer in enumerate(sorted(set(layers)))}


def _sig_col(padj_method):
    """Determine p-value column name based on multiple testing correction method.
    
    Returns the appropriate column name to use for significance filtering and plotting.
    
    Args:
        padj_method (str): P-value adjustment method name ('bonferroni', 'fdr_bh', 'none', etc.).
    
    Returns:
        str: Column name 'p.value' if method is 'none' (unadjusted), otherwise 'p.adj' (adjusted).
    """
    return "p.value" if padj_method == "none" else "p.adj"


def _fit_regression(x, y):
    """Fit linear regression line to valid (non-NaN) data points.
    
    Computes linear regression coefficients only if sufficient valid data exists.
    
    Args:
        x (array-like): X-axis values (independent variable). May contain NaN.
        y (array-like): Y-axis values (dependent variable). May contain NaN.
    
    Returns:
        array | None: Polynomial coefficients [slope, intercept] if ≥3 valid samples exist,
                     None if insufficient valid data.
    """
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 3:
        return None
    return np.polyfit(x[mask], y[mask], 1)

# Note: Currently unused — kept for potential future use (e.g., network edge hover tooltips)
def _regression_png(results, gene_A, gene_B, assayDataName):
    """Generate a regression scatter plot as a base64-encoded PNG image.
    
    Creates a multi-panel matplotlib figure with scatter plots and fitted regression lines
    for each sample category, then encodes as PNG in base64 format (suitable for embedding
    in HTML tooltips or interactive visualizations).
    
    Args:
        results (dict): Output from get_diffNetworks() containing:
            - "assayData": Dictionary mapping dataset names to expression DataFrames
            - "metadata": Sample metadata indexed by sample IDs
            - "category_variable": Column name defining sample groups
        
        gene_A (str): Gene symbol for X-axis variable.
        
        gene_B (str): Gene symbol for Y-axis variable.
        
        assayDataName (str): Name of the assay dataset in results["assayData"].
    
    Returns:
        str: Base64-encoded PNG image string (ready for embedding in HTML/JSON).
    """
    assay_df   = results["assayData"][assayDataName]
    metadata   = results["metadata"]
    cat_var    = results["category_variable"]
    categories = metadata[cat_var].unique()
    col_map    = {cat: CATEGORY_COLOURS[i % len(CATEGORY_COLOURS)]
                  for i, cat in enumerate(categories)}

    fig, axes = plt.subplots(1, len(categories),
                             figsize=(4 * len(categories), 4),
                             sharey=True)
    if len(categories) == 1:
        axes = [axes]

    # Generate scatter plots and regression lines for each category
    for ax, cat in zip(axes, categories):
        ids = metadata[metadata[cat_var] == cat].index.intersection(assay_df.columns)
        x   = assay_df.loc[gene_A, ids].values.astype(float)
        y   = assay_df.loc[gene_B, ids].values.astype(float)
        col = col_map[cat]

        ax.scatter(x, y, color=col, alpha=0.7, s=40)

        fit = _fit_regression(x, y)
        if isinstance(fit, np.ndarray):
            m, b = fit
            xl = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            ax.plot(xl, m * xl + b, color=col, linewidth=2)

        ax.set_title(cat)
        ax.set_xlabel(gene_A)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_ylabel(gene_B)
    fig.suptitle(f"{gene_A} × {gene_B}", fontweight="bold")
    plt.tight_layout()

    # Encode plot as base64 PNG for embedding in HTML
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _build_legend(category, col_map):
    """Build an HTML legend for the network visualization.
    
    Generates styled HTML legend showing the sample category and layer color codes
    for a network visualization.
    
    Args:
        category (str): Sample category name to display in legend.
        
        col_map (dict): Mapping of layer names to color hex codes.
    
    Returns:
        str: HTML-formatted legend string with category label and colored layer indicators.
    """
    dot  = lambda c, l: f'<span style="background:{c};width:12px;height:12px;border-radius:50%;display:inline-block;margin-right:4px"></span>{l}'
    items = " &nbsp; ".join(dot(c, l) for l, c in sorted(col_map.items()))
    # Format legend HTML with category and layer colors
    return (
        f'<div style="font-family:Arial;padding:8px;background:#F0F3F4;border-radius:6px">'
        f'<b>Category:</b> {category} &nbsp; <b>Layers:</b> {items}<br>'
        #f'<small>Hover over edges to see the regression plot. Edge width = number of references.</small>'
        f'</div>'
    )


def save_network_html(results, edge_df, category,
                 assayDataName=None,
                 output_file="multideggs_network.html"):
    """Generate an interactive network visualization for differential edges in a category.
    
    Creates a directed network graph using Pyvis showing gene interactions (edges) filtered
    by sample category. Edge widths are proportional to literature reference counts.
    
    Args:
        results (dict): Output dictionary from get_diffNetworks() containing:
            - "assayData": Dictionary mapping dataset names to expression DataFrames
            - "metadata": Sample metadata indexed by sample IDs
            - "category_variable": Column name defining sample groups
            - "padj_method": P-value adjustment method used in analysis
        
        edge_df (pd.DataFrame): Network edges with required columns:
            - "from": Source gene symbol
            - "to": Target gene symbol
            - "p.value" or "p.adj": Statistical p-value
            - "n_references": Number of literature references (optional)
            - "curation_effort": Curation effort score (optional)
            - "layer": Omics dataset name (optional; added if missing)
        
        category (str): Sample category name for visualization labeling.
        
        assayDataName (str, optional): Name of the assay dataset (used for single-omic datasets).
            Defaults to None.
        
        output_file (str, optional): Path to save HTML file. 
            Defaults to "multideggs_network.html".
    
    Returns:
        str: Path to the saved HTML file, or None if edge_df is empty.
    
    Raises:
        Warning: Emitted if edge_df is empty.
    """
    if edge_df.empty:
        warnings.warn(f"No significant edges found for category '{category}'.")
        return None

    # Add layer column if missing (single-omic case)
    if "layer" not in edge_df.columns:
        edge_df = edge_df.copy()
        edge_df["layer"] = assayDataName or list(results["assayData"].keys())[0]
    
    sig_col = _sig_col(results["padj_method"])
    col_map = _layer_colour_map(edge_df["layer"].astype(str).unique())

    net = Network(height="620px", width="100%", bgcolor="#ffffff",
                  font_color="#2C3E50", directed=True)
    # Barnes-Hut physics simulation for force-directed layout
    net.barnes_hut(gravity=-2000, central_gravity=0.1,
                   spring_length=170, spring_strength=0.02)
    
    gene_layers = {}
    for _, row in edge_df.iterrows():
        layer = str(row["layer"])
        for gene in (row["from"], row["to"]):
            gene_layers.setdefault(gene, set()).add(layer)

    added_nodes = set()
    seen_edges = {}
    
    for _, row in edge_df.iterrows():
        gA, gB = row["from"], row["to"]
        layer  = str(row["layer"])
        colour = col_map[layer]
        pval   = row.get(sig_col, float("nan"))
        pval_s = f"{pval:.3e}" if pd.notna(pval) else "n/a"
        nref   = int(row.get("n_references", 0) or 0)
        ceff   = int(row.get("curation_effort", 0) or 0)
        width  = max(1.5, min(nref * 0.4, 8))
        
        for gene in (gA, gB):
            if gene not in added_nodes:
                gene_layer_list = gene_layers[gene]
                if len(gene_layer_list) > 1:
                    node_bg     = "#B0BEC5"
                    node_border = "#37474F"
                    node_title  = f"{gene} (layers: {', '.join(sorted(gene_layer_list))})"
                else:
                    single_layer = next(iter(gene_layer_list))
                    node_bg      = col_map[single_layer]
                    node_border  = NODE_BORDER
                    node_title   = f"{gene} ({single_layer})"
                net.add_node(gene, label=gene,
                            color={"background": node_bg, "border": node_border,
                                    "highlight": {"background": "#F9E79F", "border": "#B7950B"}},
                            size=22, title=node_title)
                added_nodes.add(gene)

        edge_key = tuple(sorted([gA, gB]))
        count = seen_edges.get(edge_key, 0)
        seen_edges[edge_key] = count + 1
        roundness = 0.2 + count * 0.25  # 0.2, 0.45, 0.7 for successive edges
        
        net.add_edge(gA, gB, color=colour, width=width,
                    smooth={"type": "curvedCW", "roundness": roundness},
                    title=f"{gA} → {gB} | p.adj = {pval_s} | Ref: {nref} | Curation: {ceff}")

    # Write HTML and inject legend for layer colors
    net.write_html(output_file, open_browser=False)
    legend = _build_legend(category, col_map)
    # Inject legend into HTML body
    with open(output_file, "r", encoding="utf-8") as f:
        html = f.read()
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html.replace("<body>", f"<body>\n{legend}", 1))

    print(f"[MultiDEGGs] Network saved → {output_file}")
    return output_file


def plot_significant_edges(final_outputs, sig_threshold=0.05,
                   title="Significant edges by category"):
    """Generate a horizontal bar chart of significant edges ranked by -log₁₀(p-value).
    
    Creates an interactive Plotly bar chart displaying statistically significant edges,
    grouped by sample category and colored by omics layer. Bars are sorted by significance.
    
    Args:
        final_outputs (dict | pd.DataFrame): Network results with edges containing columns:
            - "from": Source gene symbol
            - "to": Target gene symbol  
            - "p.value" or "p.adj": Statistical p-value or adjusted p-value
            - "category": Sample category/group name (required for dict input)
            - "layer": Omics dataset name (optional, for coloring)
            - "n_references": Number of literature references (optional)
            - "curation_effort": Curation effort score (optional)
        
        sig_threshold (float, optional): P-value significance threshold for display.
            Only edges with p-value below threshold are shown. Defaults to 0.05.
        
        title (str, optional): Figure title. Defaults to "Significant edges by category".
    
    Returns:
        plotly.graph_objects.Figure: Interactive bar chart with:
            - X-axis: -log₁₀(p-value)
            - Y-axis: Edge labels grouped by category
            - Colors: Mapped to omics layers (if present)
            - Hover: Full statistics including p-value, references, curation effort
            - Threshold line: Vertical dashed line at -log₁₀(sig_threshold)
            Returns empty figure if no significant edges found.
    
    Raises:
        Warning: If no edges exist or if no edges pass significance threshold.
    """
    # Normalise to a single DataFrame with a 'category' column
    if isinstance(final_outputs, dict):
        frames = []
        for cat, df in final_outputs.items():
            if isinstance(df, pd.DataFrame) and not df.empty:
                tmp = df.copy()
                tmp["category"] = cat
                frames.append(tmp)
        if not frames:
            warnings.warn("No significant edges to plot.")
            return go.Figure()
        combined = pd.concat(frames, ignore_index=True)
    else:
        combined = final_outputs.copy() if isinstance(final_outputs, pd.DataFrame) \
                   else pd.DataFrame()

    if combined.empty:
        warnings.warn("No significant edges to plot.")
        return go.Figure()

    sig_col  = "p.adj" if "p.adj" in combined.columns else "p.value"
    combined = combined[combined[sig_col] < sig_threshold].copy()

    if combined.empty:
        warnings.warn(f"No edges below threshold {sig_threshold}.")
        return go.Figure()

    combined["neg_log_p"]  = -np.log10(combined[sig_col].clip(lower=1e-300))
    combined["edge_label"] = combined["from"] + "  →  " + combined["to"]

    has_layer  = "layer" in combined.columns
    col_map    = _layer_colour_map(combined["layer"].astype(str).unique()) \
                 if has_layer else {"": LAYER_COLOURS[0]}
    categories = combined["category"].unique().tolist() \
                 if "category" in combined.columns else [""]

    fig  = go.Figure()
    seen = set()

    for cat in categories:
        sub = combined[combined["category"] == cat].sort_values("neg_log_p", ascending=True)
        for _, row in sub.iterrows():
            layer  = str(row["layer"]) if has_layer else ""
            colour = col_map.get(layer, LAYER_COLOURS[0])
            nref   = int(row.get("n_references", 0) or 0)
            ceff   = int(row.get("curation_effort", 0) or 0)
            # Format hover information with full edge statistics
            hover  = (
                f"<b>{row['edge_label']}</b><br>"
                f"Category: {cat}" + (f" | Layer: {layer}" if layer else "") + "<br>"
                f"{sig_col} = {row[sig_col]:.3e} | -log10(p) = {row['neg_log_p']:.2f}<br>"
                f"Ref: {nref} | Curation: {ceff}<extra></extra>"
                )
            # Add horizontal bar for each edge
            fig.add_trace(go.Bar(
                x=[row["neg_log_p"]], y=[f"[{cat}]  {row['edge_label']}"],
                orientation="h", marker_color=colour,
                name=layer or "edges", legendgroup=layer,
                showlegend=layer not in seen,
                hovertemplate=hover,
                text=f"  {row[sig_col]:.2e}", textposition="outside",
                ))
            seen.add(layer)
    
    # Add vertical threshold line
    fig.add_vline(x=-np.log10(sig_threshold), line_dash="dash",
                  line_color="#E74C3C",
                  annotation_text=f"p = {sig_threshold}",
                  annotation_position="top right",
                  annotation_font_color="#E74C3C")

    fig.update_layout(
        title=title,
        xaxis_title="-log₁₀(p.adj)",
        bargap=0.25,
        height=max(300, 36 * len(combined) + 100),
    )
    fig.update_xaxes(showgrid=True, zeroline=False)
    fig.update_yaxes(showgrid=False, autorange="reversed")

    return fig


def plot_gene_pair_regression(results, gene_A, gene_B, assayDataName=None):
    """Generate multi-panel regression scatter plots across sample categories.
    
    Creates a Plotly figure with one subplot per category, showing the relationship
    between two genes with scatter points and fitted regression lines.
    
    Args:
        results (dict): Output dictionary from get_diffNetworks() containing:
            - "assayData": Dictionary mapping dataset names to expression DataFrames
                          (features × samples)
            - "metadata": Sample metadata DataFrame indexed by sample IDs
            - "category_variable": Column name in metadata defining sample groups
        
        gene_A (str): Gene symbol for X-axis (predictor variable). 
            Must exist in assayData rows.
        
        gene_B (str): Gene symbol for Y-axis (response variable). 
            Must exist in assayData rows.
        
        assayDataName (str, optional): Name of the assay dataset to use. 
            If None, uses the first available dataset. Defaults to None.
    
    Returns:
        plotly.graph_objects.Figure: Multi-panel scatter plot with:
            - Subplots: One per category (rows=1, cols=number of categories)
            - Markers: Individual samples colored by category
            - Regression lines: Linear fit line per category (if ≥3 valid samples)
            - Shared Y-axis: For easy comparison across categories
            - Title: Formatted as "gene_A × gene_B  [assayDataName]"
    
    Raises:
        ValueError: If gene_A or gene_B not found in the assayData rows.
    """
    if assayDataName is None:
        assayDataName = list(results["assayData"].keys())[0]

    assay_df   = results["assayData"][assayDataName]
    metadata   = results["metadata"]
    cat_var    = results["category_variable"]
    categories = metadata[cat_var].unique().tolist()

    # Validate that requested genes exist in the data
    for gene in (gene_A, gene_B):
        if gene not in assay_df.index:
            raise ValueError(f"Gene '{gene}' not found in '{assayDataName}'.")

    col_map = {cat: CATEGORY_COLOURS[i % len(CATEGORY_COLOURS)]
               for i, cat in enumerate(categories)}

    # Create subplots - one per category with shared Y-axis
    fig = make_subplots(rows=1, cols=len(categories),
                        subplot_titles=categories,
                        shared_yaxes=True)

    # Add scatter plot and regression line for each category
    for i, cat in enumerate(categories, start=1):
        ids    = metadata[metadata[cat_var] == cat].index.intersection(assay_df.columns)
        x      = assay_df.loc[gene_A, ids].values.astype(float)
        y      = assay_df.loc[gene_B, ids].values.astype(float)
        colour = col_map[cat]

        # Add scatter points for this category
        fig.add_trace(
            go.Scatter(x=x, y=y, mode="markers",
                       marker=dict(color=colour, size=7, opacity=0.75,
                                   line=dict(color="white", width=0.5)),
                       name=cat, showlegend=True,
                       hovertemplate=(
                           f"{gene_A}: %{{x:.2f}}<br>"
                           f"{gene_B}: %{{y:.2f}}"
                           f"<extra>{cat}</extra>"
                       )),
            row=1, col=i
        )

        # Add regression line if sufficient data points
        fit = _fit_regression(x, y)
        if isinstance(fit, np.ndarray):
            m, b = fit
            xl = np.linspace(np.nanmin(x), np.nanmax(x), 100)
            fig.add_trace(
                go.Scatter(x=xl, y=m * xl + b, mode="lines",
                           line=dict(color=colour, width=2.5),
                           showlegend=False, hoverinfo="skip"),
                row=1, col=i
            )

        fig.update_xaxes(title_text=gene_A, row=1, col=i,
                         showgrid=True, zeroline=False)

    fig.update_layout(
        title=f"{gene_A} × {gene_B}  [{assayDataName}]",
        height=420,
    )
    fig.update_yaxes(title_text=gene_B, row=1, col=1,
                     showgrid=True, zeroline=False)

    return fig

def show_plot(analysis_output, plot_type="summary", category=None,
              genes=None, gene_A=None, gene_B=None, **kwargs):
    """Display different types of analytical plots based on user request.
    
    Unified interface for generating various visualization types from MultiDEGGs analysis:
    summary (edge bar chart), regression (scatter with regression lines), or network
    (interactive Pyvis graph).
    
    Args:
        analysis_output (tuple): Output from run_multideggs() containing (results, final_outputs).
        
        plot_type (str, optional): Type of visualization to create:
            - "summary": Horizontal bar chart of significant edges by -log₁₀(p-value)
            - "regression": Multi-panel scatter plots with regression lines
            - "network": Interactive network graph visualization
            Defaults to "summary".
        
        category (str, optional): Sample category for network visualization.
            Required for plot_type="network". Defaults to None.
        
        genes (str | list, optional): Specific genes to focus on in network.
            Generates sub-network including only these genes. Defaults to None.
        
        gene_A (str, optional): Gene symbol for X-axis in regression plot.
            Required for plot_type="regression". Defaults to None.
        
        gene_B (str, optional): Gene symbol for Y-axis in regression plot.
            Required for plot_type="regression". Defaults to None.
        
        **kwargs: Additional keyword arguments passed to plot functions
            (e.g., output_file, sig_threshold, title, assayDataName).
    
    Returns:
        plotly.graph_objects.Figure | str | None: Depending on plot_type:
            - "summary": Plotly Figure object
            - "regression": Plotly Figure object
            - "network": Path to saved HTML file, or None if no edges found
    
    Raises:
        ValueError: If required parameters missing for selected plot_type or invalid plot_type.
        Warning: If no edges found for network visualization or genes.
    """
    results, final_outputs = analysis_output

    # Route to appropriate plotting function based on plot_type
    if plot_type == "summary":
        return plot_significant_edges(final_outputs, **kwargs)

    elif plot_type == "regression":
        if not gene_A or not gene_B:
            raise ValueError("For 'regression' plot, specify both 'gene_A' and 'gene_B'.")
        return plot_gene_pair_regression(results, gene_A, gene_B,
                                         assayDataName=kwargs.get("assayDataName"))

    elif plot_type == "network":
        if category is None:
            raise ValueError("Specify the 'category' parameter (e.g. category='Tumor').")

        # Extract edges for the specified category
        if isinstance(final_outputs, dict):
            edge_df = final_outputs.get(category, pd.DataFrame()).copy()
        else:
            edge_df = final_outputs[final_outputs["category"] == category].copy()

        if edge_df.empty:
            warnings.warn(f"No significant edges found for category '{category}'.")
            return None

        # Focus mode: filter to specific genes
        if genes is not None:
            if isinstance(genes, str):
                genes = [genes]
            edge_df = edge_df[edge_df["from"].isin(genes) | edge_df["to"].isin(genes)].copy()
            if edge_df.empty:
                print(f"[MultiDEGGs] No interactions found for genes {genes} in category '{category}'.")
                return None
            output_file = kwargs.pop("output_file", f"multideggs_network_{category}_focus.html")
            print(f"[MultiDEGGs] Focus mode → generating sub-network for: {genes}")
            path = save_network_html(results, edge_df, category=category,
                                     output_file=output_file, **kwargs)
            try:
                from IPython.display import IFrame, display
                display(IFrame(path, width="100%", height=650))
            except Exception:
                pass
            return path

        # Global mode: full network for the category
        else:
            output_file = kwargs.pop("output_file", f"multideggs_network_{category}_global.html")
            print(f"[MultiDEGGs] Global mode → generating full network ({len(edge_df)} edges).")
            path = save_network_html(results, edge_df, category=category,
                                     output_file=output_file, **kwargs)
            try:
                from IPython.display import IFrame, display
                display(IFrame(path, width="100%", height=650))
            except Exception:
                pass
            return path

    else:
        raise ValueError("Invalid plot_type. Choose from: 'summary', 'network', 'regression'.")