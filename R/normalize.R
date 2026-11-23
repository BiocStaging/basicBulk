#' @title tidySE
#' @description Tidy and normalize a \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'  by adding assays for downstream analyses and visualization. Genes are filtered with
#'  \link[edgeR]{filterByExpr}, TMM normalization factors are calculated with \link[edgeR]{normLibSizes}, and the
#'  expression (\code{E}) matrix from \link[limma]{voomWithQualityWeights} (the quality weights are not stored)
#'  becomes the \code{logCPM} assay, with \code{NA} for the filtered genes. When \code{adjustBatch = TRUE}, a
#'  \code{ComBat} assay holds voom values from \link[sva]{ComBat_seq}-adjusted counts.
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} with a
#'  \code{counts} assay
#' @param group Character specifying the colData column to tidy into the \code{Group} factor
#' @param covariates Character vector of colData columns with numeric covariates
#' @param block Character vector of colData columns with blocking factors
#' @param groups Named character vector for group factor recoding
#' @param adjustBatch Logical indicating whether a \code{Batch} colData column should be adjusted for
#'  with \link[sva]{ComBat_seq} in an additional \code{ComBat} assay
#' @return The input \code{se} with an added \code{logCPM} (and \code{ComBat}) assay, a
#'  \code{retained} \code{rowData} column, and \code{Group}/\code{Sample}/\code{lib.size}/\code{norm.factors}
#'  in the \code{colData}
#' @details The input \code{se} must have a \code{counts} assay. Downstream basicBulk functions expect a
#'  \code{rowData} \code{ID} column and, for symbol-based plots and labelled heatmaps, a \code{symbol}
#'  column. The \code{logCPM}/\code{ComBat} assays are \code{NA} for genes that did not pass
#'  \link[edgeR]{filterByExpr}.
#' @import tidySummarizedExperiment
#' @import edgeR
#' @import limma
#' @import SummarizedExperiment
#' @importFrom dplyr select mutate recode_factor
#' @importFrom tibble as_tibble
#' @importFrom purrr pluck set_names
#' @importFrom forcats as_factor
#' @importFrom cli cli_h1 cli_progress_step cli_warn
#' @importFrom magrittr %>% extract2
#' @seealso \link[augere.de]{runVoom}
#' @references
#' Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways: differential expression analysis of RNA-Seq
#'  experiments using Rsubread and the edgeR quasi-likelihood pipeline. \emph{F1000Research}, 5, 1438.
#'  \doi{10.12688/f1000research.8987.2}
#'
#' Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq
#'  data. \emph{Genome Biology}, 11(3), R25. \doi{10.1186/gb-2010-11-3-r25}
#'
#' Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2025). edgeR v4: powerful differential analysis of sequencing data
#'  with expanded functionality and improved support for small counts and larger datasets. \emph{Nucleic Acids
#'  Research}, 53(2), gkaf018. \doi{10.1093/nar/gkaf018}
#'
#' Law CW, Chen Y, Shi W, Smyth GK (2014). voom: precision weights unlock linear model analysis tools for RNA-seq
#'  read counts. \emph{Genome Biology}, 15(2), R29. \doi{10.1186/gb-2014-15-2-r29}
#'
#' Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat ML, Smyth GK, Ritchie ME (2015). Why
#'  weight? Modelling sample and observational level variability improves power in RNA-seq analyses. \emph{Nucleic
#'  Acids Research}, 43(15), e97. \doi{10.1093/nar/gkv412}
#'
#' Zhang Y, Parmigiani G, Johnson WE (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. \emph{NAR
#'  Genomics and Bioinformatics}, 2(3), lqaa078. \doi{10.1093/nargab/lqaa078}
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' tidySE(se0, group = "cell", block = "dex")
tidySE <- function(se,
                   group = "Group",
                   covariates = NULL,
                   block = NULL,
                   groups = se %>%
                     SummarizedExperiment::colData() %>%
                     as.data.frame() %>%
                     purrr::pluck(group) %>%
                     forcats::as_factor() %>%
                     levels() %>%
                     purrr::set_names(),
                   adjustBatch = FALSE){

  cli::cli_h1("tidySE")

  stopifnot("tidySE requires unique gene (row) names" = !anyDuplicated(rownames(se)))

  # Add Sample
  se$Sample <- colnames(se)

  # Remove any colData ID to prevent a collision with the gene ID in tidy conversions
  if(!is.null(colData(se)$ID)){
    colData(se) <- se %>%
      SummarizedExperiment::colData() %>%
      subset(select = -ID)
  }

  if(is.numeric(SummarizedExperiment::colData(se)[[group]])){
    cli::cli_warn("{group} was numeric and will be now be coded as a factor")
  }

  # Recode group
  se$Group <- se %>%
    SummarizedExperiment::colData() %>%
    tibble::as_tibble() %>%
    dplyr::select(exp = !!group) %>%
    dplyr::mutate(exp = dplyr::recode_factor(exp, !!!groups)) %>%
    purrr::pluck("exp")

  cli::cli_progress_step("Preparing design matrix")
  design <- se %>%
    SummarizedExperiment::colData() %>%
    .designMatrix(groups = "Group", covariates = covariates, block = block)

  cli::cli_progress_step("Filtering low-expression genes and computing normalization factors")
  y <- se %>%
    SummarizedExperiment::assay("counts") %>%
    edgeR::DGEList()
  retained <- edgeR::filterByExpr(y, group = se$Group, large.n = 0)
  y <- y[retained, ] %>%
    edgeR::normLibSizes()

  rowData(se)$retained <- retained
  se$lib.size <- y$samples$lib.size
  se$norm.factors <- y$samples$norm.factors

  cli::cli_progress_step("Adding logCPM assay")
  assay(se, "logCPM") <- y %>%
    limma::voomWithQualityWeights(design = design, plot = FALSE) %>%
    magrittr::extract2("E") %>%
    .padAssay(se)

  if(adjustBatch == TRUE){
    stopifnot(!is.null(se$Batch),
              rlang::is_installed("sva"))

    if(!is.null(covariates)){
      cli::cli_warn("Covariates can not be properly included in the batch adjusted design")
    }

    design <- se %>%
      SummarizedExperiment::colData() %>%
      .designMatrix(groups = "Group",
                    covariates = covariates,
                    block = block %>%
                      .[. != "Batch"] %>%
                      {if(identical(., character(0))){
                        NULL
                      }else{
                        .
                      }})

    cli::cli_progress_step("Adding ComBat assay (voom of ComBat-seq adjusted counts)")
    assay(se, "ComBat") <- se %>%
      SummarizedExperiment::assay("counts") %>%
      .[retained, ] %>%
      sva::ComBat_seq(counts = .,
                      batch = se$Batch,
                      group = se$Group) %>%
      edgeR::DGEList() %>%
      edgeR::normLibSizes() %>%
      limma::voomWithQualityWeights(design = design, plot = FALSE) %>%
      magrittr::extract2("E") %>%
      .padAssay(se)
  }

  return(se)
}

#' @title .designMatrix
#' @description Build an additive design matrix from colData columns, following the public
#'  \link[augere.de]{processSimpleDesignMatrix} approach. Groups become a no-intercept factor block,
#'  covariates are coerced numeric, and blocking variables become factors.
#' @param cd A \link[S4Vectors:DataFrame-class]{DataFrame} or data.frame of sample metadata (colData)
#' @param groups String naming the grouping-factor column (NULL gives an intercept model)
#' @param covariates Character vector of numeric covariate columns
#' @param block Character vector of blocking-factor columns
#' @return A design matrix
#' @importFrom stats model.matrix as.formula
#' @importFrom purrr set_names map
#' @importFrom rlang %||%
#' @importFrom magrittr %>% extract2
#' @keywords internal
#'
.designMatrix <- function(cd,
                          groups = NULL,
                          covariates = NULL,
                          block = NULL){

  cd <- cd %>%
    as.data.frame()

  model.data <- c(if(!is.null(groups)){list(group. = cd %>% magrittr::extract2(groups) %>% factor())},
                  (covariates %||% character(0)) %>%
                    purrr::set_names(paste0("covariate.", make.names(., unique = TRUE))) %>%
                    purrr::map(~ cd %>% magrittr::extract2(.x) %>% as.numeric()),
                  (block %||% character(0)) %>%
                    purrr::set_names(paste0("block.", make.names(., unique = TRUE), ".")) %>%
                    purrr::map(~ cd %>% magrittr::extract2(.x) %>% factor()))

  if(length(model.data) == 0){
    return(stats::model.matrix(~ 1, data = cd))
  }

  model.data %>%
    names() %>%
    paste(collapse = " + ") %>%
    paste(if(!is.null(groups)){"~ 0 +"}else{"~"}, .) %>%
    stats::as.formula() %>%
    stats::model.matrix(data = model.data)
}

#' @title .padAssay
#' @description Expand a matrix of values for a subset of genes to the full gene set of \code{se},
#'  filling absent genes (e.g. those filtered out before voom) with \code{NA}. Values are inserted by
#'  their original row identities, not by position.
#' @param mat A matrix of values whose rownames are a subset of the genes of \code{se}
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @return A matrix with the dimensions and dimnames of \code{se}
#' @keywords internal
#'
.padAssay <- function(mat,
                      se){
  full <- matrix(NA_real_,
                 nrow = nrow(se),
                 ncol = ncol(se),
                 dimnames = list(rownames(se), colnames(se)))
  full[rownames(mat), colnames(mat)] <- mat
  return(full)
}
