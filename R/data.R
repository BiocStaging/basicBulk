#' @title airwayDGE example data
#' @description Example `DGEdata` (the tidy differential-expression object
#'  produced by [makeDGEdata()]) for the basicBulk workflow, built by testing
#'  three cell lines against a shared reference (`N61311`, `N061011`, and
#'  `N080611` vs `N052611`) on the expression-filtered genes of the full
#'  normalized airway data, then subset to the most variable retained genes plus
#'  a few genes used in examples. It is a
#'  named list of one tibble per contrast, with columns including `ID`,
#'  `symbol`, `LogFC`, `AveExpr`, `t`, `PValue`, `FDR`, and `B`.
#' @format A named `list` of three tibbles, one per contrast.
#' @docType data
#' @return A named list of tibbles.
#' @usage data(airwayDGE)
#' @source \link[airway]{airway}
#'
"airwayDGE"

#' @title airwayGSEA example data
#' @description Example tidied GSEA results (the object produced by
#'  [tidyGSEAlist()]) for the basicBulk workflow, built by running
#'  camera on the expression-filtered genes of the full normalized airway
#'  data against the Reactome gene sets
#'  from [getReactome()] for each of the three contrasts. It is a named list of
#'  one tibble per contrast, with columns `Database`, `Term`, `NumGenes`,
#'  `Direction`, `PValue`, and `FDR`.
#' @format A named `list` of three tibbles, one per contrast.
#' @docType data
#' @return A named list of tibbles.
#' @usage data(airwayGSEA)
#' @source \link[airway]{airway}
#'
"airwayGSEA"
