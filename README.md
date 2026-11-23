
# basicBulk

<!-- badges: start -->

<!-- badges: end -->

basicBulk provides tidyverse functions, tables, and plots for end-to-end analyses of (pseudo)bulk RNA-seq data within the Bioconductor ecosystem, with a focus on comparing contrasts to each other. Starting from a SummarizedExperiment of counts, it covers normalization and quality control, differential expression, and gene set enrichment. Alongside per-contrast volcano, MA, and heatmap plots, it summarizes the overlap of differentially expressed genes across contrasts (UpSet plots), compares pairs of contrasts gene by gene (4-way plots), correlates their logFC profiles, and compares top pathways. basicBulk also includes a Claude Code skill that enables agentic end-to-end analyses.

## Installation

basicBulk can be installed from [Bioconductor](https://bioconductor.org/packages/basicBulk):

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("basicBulk")
```

To install the development version directly from GitHub:

```r
if(!requireNamespace("remotes", quietly = TRUE)){
    install.packages("remotes")
  }

remotes::install_github("Genentech/basicBulk")
```

## Examples

See the package [vignette](https://bioconductor.org/packages/release/bioc/vignettes/basicBulk/inst/doc/basicBulk.html).

## Claude Code plugin

basicBulk also provides a Claude Code plugin with a skill for bulk RNA-seq analyses that uses the package's functions and can generate a complete R Markdown report for a dataset. The skill runs R, so basicBulk, augere.de, and augere.gsea need to be installed.

To install the plugin, add basicBulk's Bioconductor GitHub mirror as a plugin marketplace and then install the plugin from within a Claude Code session:

```
/plugin marketplace add bioc/basicbulk
/plugin install basicbulk@basicbulk
```

The same plugin is also available from the development repository (`/plugin marketplace add genentech/basicbulk`). Add only one of the two marketplaces, since both are named `basicbulk`. It's recommended to enable auto-update for the `basicbulk` marketplace (`/plugin`, under Marketplaces) so new versions of the plugin install automatically.
