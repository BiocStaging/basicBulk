#' basicBulk
#'
#' basicBulk provides basic functions that leverage the
#'  tidyverse for analyses of (pseudo)bulk omic data.
#' @keywords internal
#' @importFrom augere.de runVoom
#' @importFrom augere.gsea runContrast
"_PACKAGE"

.onAttach <- function(libname = find.package("basicBulk"), pkgname = "basicBulk") {
  packageStartupMessage(paste0("basicBulk v", getNamespaceVersion(pkgname)," has loaded."))
}

utils::globalVariables(c(
  ".", ":=", "-log10.FDR", ".feature", "AveExpr", "Contrast", "Correlation",
  "count", "DEGs", "Database", "Direction", "Down", "FDR", "Group",
  "HGNC.symbol", "ID", "Loading", "LogFC", "PC", "PC1", "PC2", "Pathway",
  "Sample", "Term", "Up", "color", "contrast", "description", "direction",
  "gene_id", "gene_name", "group1", "group2", "hjust", "name", "p.adj",
  "p.adj.signif", "r", "rowname", "set", "symbol", "term", "total", "value",
  "vjust", "x", "y", "y.position", "ensembl_gene", "gs_name"
))

## usethis namespace: start
## usethis namespace: end
