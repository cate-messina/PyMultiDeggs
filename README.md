# PyMultiDEGGs

## Overview

**PyMultiDEGGs** is a Python implementation of the original **multiDEGGs R package** by Elisabetta Sciacca et al.

The package performs **differential network analysis** to identify gene–gene and, more generally, entity–entity interactions whose association structure differs between sample categories. It can be applied to both single-omic and multi-omic datasets.

Rather than focusing only on individual molecular features that change between conditions, PyMultiDEGGs investigates how relationships between molecular entities change across groups of samples. These differential interactions can provide insights into biological mechanisms and can also be used as engineered features for downstream predictive modeling.

---

## Installation

### Install from GitHub

The package can be installed directly from GitHub:

```bash
pip install git+https://github.com/<your-username>/PyMultiDEGGs.git
```

Alternatively, clone the repository and install it locally:

```bash
git clone https://github.com/<your-username>/PyMultiDEGGs.git
cd PyMultiDEGGs
pip install .
```

---

## Requirements

* Python ≥ 3.10

The package relies on the following main Python libraries:

* pandas
* numpy
* statsmodels
* joblib
* requests
* beautifulsoup4
* matplotlib
* plotly
* pyvis
* IPython

Dependencies are automatically installed when installing PyMultiDEGGs through `pip`.

---

## Quick Start

### Import

```python
import multideggs as md
```

> **Note:** `metadata` must be provided as a pandas DataFrame with sample IDs as the index. For a single-omic dataset, `assayData` contains features (e.g., genes) in rows and samples in columns. The sample IDs must match between `assayData` and `metadata`.

### Single-omic analysis

For a single-omic dataset, provide the molecular data as a DataFrame:

```python
analysis = md.run_multideggs(
    assayData=rnaseq_df,
    metadata=metadata,
    category_variable="response",
    sig_threshold=0.05,
    verbose=True,
)
```

The analysis returns two objects:

```python
results, final_outputs = analysis
```

* `results` contains the complete differential-network analysis.
* `final_outputs` contains the significant differential interactions selected using the specified `sig_threshold`.

### Multi-omic analysis

Multiple omic layers can be provided as a dictionary of datasets:

```python
assayData = {
    "RNAseq": rnaseq_df,
    "Proteomics": proteomics_df,
    "Olink": olink_df,
}
```

The same analysis function can then be used:

```python
results, final_outputs = md.run_multideggs(
    assayData=assayData,
    metadata=metadata,
    category_variable="response",
    sig_threshold=0.05,
    verbose=True,
)
```

PyMultiDEGGs automatically detects whether the input contains a single omic layer or multiple layers. For multi-omic analyses, the resulting interactions retain information about their corresponding omic layer.

---

## Main Analysis Function

The main entry point of PyMultiDEGGs is:

```python
analysis = md.run_multideggs(
    assayData=assayData,
    metadata=metadata,
    category_variable="response",
    sig_threshold=0.05,
)
```

### Parameters

| Parameter           | Description                                                                        |        Default |
| ------------------- | ---------------------------------------------------------------------------------- | -------------: |
| `assayData`         | Omics data provided as a single dataset or as multiple omic layers.                |              — |
| `metadata`          | Sample annotations indexed by sample IDs.                                          |              — |
| `category_variable` | Column in `metadata` defining the sample groups.                                   |              — |
| `category_subset`   | Optional subset of categories to retain for the analysis.                          |         `None` |
| `padj_method`       | Multiple-testing adjustment method.                                                | `"bonferroni"` |
| `sig_threshold`     | Adjusted p-value threshold used to select significant interactions.                |         `0.05` |
| `percentile_vector` | Percentile thresholds evaluated during the percolation analysis.                   |         `None` |
| `regression_method` | Regression method used for edge testing: `"lm"` or `"rlm"`.                        |         `"lm"` |
| `verbose`           | Whether to print workflow progress and warnings.                                   |         `True` |
| `organism`          | Taxonomic ID of the interaction network.                                           |         `9606` |
| `archive_version`   | OmniPath archive version to use. If `None`, the cached or bundled network is used. |         `None` |

### Example

```python
analysis = md.run_multideggs(
    assayData=assayData,
    metadata=metadata,
    category_variable="response",
    category_subset=None,
    padj_method="bonferroni",
    sig_threshold=0.05,
    percentile_vector=None,
    regression_method="lm",
    verbose=True,
    organism=9606,
    archive_version=None,
)
```

---

## Interaction Network

PyMultiDEGGs uses **OmniPath** as the source of molecular interaction information.

The `organism` parameter specifies the taxonomic identifier used to retrieve the interaction network. The default is human:

```python
organism=9606
```

For example, mouse can be specified with:

```python
analysis = md.run_multideggs(
    assayData=assayData,
    metadata=metadata,
    category_variable="response",
    organism=10090,
)
```

### OmniPath Archive Versions

A specific OmniPath archive version can be selected using the `archive_version` parameter:

```python
analysis = md.run_multideggs(
    assayData=assayData,
    metadata=metadata,
    category_variable="response",
    organism=9606,
    archive_version="20250601",
)
```

If `archive_version=None`, PyMultiDEGGs uses the cached or bundled interaction network.

The interaction data retrieved from OmniPath provide information about properties such as interaction type, directionality, regulatory effect, literature support, and curation.

> **Note:** OmniPath archive schemas have changed over time. Older archive versions may not contain fields introduced in later versions, such as `consensus_direction`, `consensus_stimulation`, `consensus_inhibition`, and `curation_effort`. In these cases, the corresponding columns are **absent from the network data** rather than being populated with missing values.

The interaction network can also be explicitly updated using:

```python
md.update_network(
    organism=9606,
    archive_version="20250601",
)
```

---

## Results

The analysis returns a tuple containing two objects:

```python
results, final_outputs = analysis
```

### `results`

`results` contains the complete differential-network analysis, including the tested interactions and their statistical results.

### `final_outputs`

`final_outputs` contains the significant differential interactions selected using the specified adjusted p-value threshold.

For multi-omic analyses, the results also retain information about the corresponding omic layer.

### Final Results Table

The final output combines the statistical results of the differential analysis with interaction-level annotations obtained from OmniPath.

| Column                  | Description                                                                      |
| ----------------------- | -------------------------------------------------------------------------------- |
| `category`              | Sample category associated with the differential interaction.                    |
| `from`                  | Source entity of the interaction.                                                |
| `to`                    | Target entity of the interaction.                                                |
| `interaction_type`      | Type or category of the interaction represented in the network.                  |
| `n_references`          | Number of distinct literature references associated with the interaction.        |
| `curation_effort`       | Number of unique database–citation pairs supporting the interaction in OmniPath. |
| `is_stimulation`        | Indicates whether the interaction is annotated as stimulatory.                   |
| `consensus_stimulation` | Consensus annotation regarding stimulation.                                      |
| `is_inhibition`         | Indicates whether the interaction is annotated as inhibitory.                    |
| `consensus_inhibition`  | Consensus annotation regarding inhibition.                                       |
| `is_directed`           | Indicates whether the interaction is directed from the source to the target.     |
| `consensus_direction`   | Consensus information regarding the directionality of the interaction.           |
| `PMDI_interaction`      | Representation or identifier of the interaction used by PyMultiDEGGs.            |
| `p.value`               | P-value from the statistical test for differential association.                  |
| `p.adj`                 | P-value adjusted for multiple testing.                                           |
| `layer`                 | Omic layer associated with the interaction.                                      |
| `url_from`              | NCBI Gene URL associated with the source entity.                                 |
| `url_to`                | NCBI Gene URL associated with the target entity.                                 |

The OmniPath-derived fields describe properties of the underlying biological interactions, whereas `p.value`, `p.adj`, `category`, and `layer` describe the PyMultiDEGGs analysis.

---

## Visualization

PyMultiDEGGs provides functions for visualizing the analysis results.

### Summary Plot

```python
md.show_plot(
    analysis,
    plot_type="summary",
)
```

### Differential Network

To visualize the network for a specific category:

```python
md.show_plot(
    analysis,
    plot_type="network",
    category="Responder",
)
```

Specific entities can be highlighted:

```python
md.show_plot(
    analysis,
    plot_type="network",
    category="Responder",
    genes=["TNF", "TNFRSF1A"],
)
```

### Regression Plot

The relationship between two interacting entities can be inspected using:

```python
md.show_plot(
    analysis,
    plot_type="regression",
    gene_A="TNF",
    gene_B="TNFRSF1A",
    assayDataName="Proteomics",
)
```

### Interactive HTML Results

Interactive results can be generated with:

```python
md.show_html_results(analysis)
```

---

## Validation

The Python implementation was validated against the original R implementation of multiDEGGs using the same analysis workflow and input data.

The resulting differential network analyses were compared to assess consistency between the two implementations.

---

## Project Structure

```text
PyMultiDEGGs/
├── pyproject.toml
├── README.md
├── LICENSE
├── examples/
├── tests/
└── src/
    └── multideggs/
        ├── __init__.py
        ├── ...
        └── data/
```

---

## Citation

If you use PyMultiDEGGs, please cite the original publications:

**Sciacca E, Alaimo S, Silluzio G, Ferro A, Latora V, Pitzalis C, Pulvirenti A, Lewis MJ (2023).**
*DEGGs: an R package with shiny app for the identification of differentially expressed gene–gene interactions in high-throughput sequencing data.*
Bioinformatics, 39(4), btad192.
https://doi.org/10.1093/bioinformatics/btad192

**Sciacca E, Wang SS, Pitzalis C, Lewis MJ (2026).**
*multiDEGGs: Single or Multiomic Differential Network Analysis for Biomarker Discovery and Feature Engineering for Predictive Modeling.*
Computational and Structural Biotechnology Journal, 35(1), 0001.
https://doi.org/10.34133/csbj.0001

The Python implementation is based on the original multiDEGGs R package.

A citation for the Python implementation will be added once the repository is publicly released.

---

## License

This project is distributed under the **GNU General Public License v3.0 (GPL-3.0)**.

See `LICENSE` for the full license text.
