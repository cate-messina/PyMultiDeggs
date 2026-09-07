# MultiDEGGs Documentation

Python implementation of differential network analysis for single and multi-omic data, based on the original [multiDEGGs R package](https://github.com/elisabettasciacca/multiDEGGs) by Elisabetta Sciacca et al.

See the [README](https://github.com/cate-messina/PyMultiDEGGs) for installation and quickstart examples.

## Pipeline

::: multideggs.wrapper_function.run_multideggs

## Core analysis functions

::: multideggs.core_functions.tidy_metadata
::: multideggs.core_functions.get_diffNetworks_singleOmic
::: multideggs.core_functions.get_diffNetworks
::: multideggs.core_functions.get_multiOmics_diffNetworks
::: multideggs.core_functions.get_sig_degg
::: multideggs.core_functions.calc_pvalue_network
::: multideggs.core_functions.calc_pvalue_percentile
::: multideggs.core_functions.single_percentile_calc
::: multideggs.core_functions.percolation_analysis

## Network management

::: multideggs.omnipath_network.load_network
::: multideggs.omnipath_network.update_network

## Visualization

::: multideggs.plot_functions.show_plot
::: multideggs.plot_functions.plot_significant_edges
::: multideggs.plot_functions.plot_gene_pair_regression
::: multideggs.plot_functions.save_network_html

## Results tables

::: multideggs.final_table.show_html_results
::: multideggs.final_table.MultiDEGGsTable
::: multideggs.final_table.MultiDEGGsTableDict