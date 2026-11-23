#' @title printSigTibble
#' @description Print the sigTibble
#' @param sigTibble A tibble of tidied significant DGE results
#' @param rearrange Logical indicating whether to sort by contrast and FDR
#' @param BiocStyle Logical indicating whether the table is rendered in a BiocStyle vignette, where the per-column filter boxes misalign, so they are dropped in favour of the global search
#' @return A datatable
#' @import DT
#' @importFrom dplyr select rename_with case_when arrange
#' @importFrom tidyselect any_of
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayDGE)
#' airwayDGE %>%
#'     makeSigTibble() %>%
#'     printSigTibble()
printSigTibble <- function(sigTibble,
                           rearrange = TRUE,
                           BiocStyle = FALSE){

  sigTibble %>%
    dplyr::select(tidyselect::any_of(c("contrast",
                                       "dupe_count",
                                       "symbol",
                                       "desc",
                                       "LogFC",
                                       "PValue",
                                       "FDR"))) %>%
    dplyr::rename_with(~ dplyr::case_when(
      . ==  "contrast" ~ "Contrast",
      . ==  "dupe_count" ~ "Hits",
      . ==  "symbol" ~ "Symbol",
      . == "desc" ~ "Description",
      TRUE ~ .)) %>%
    {if(rearrange == TRUE){
      dplyr::arrange(., Contrast, FDR)
    }else if(rearrange == FALSE){
      .
    }} %>%
    DT::datatable(filter = ifelse(BiocStyle == FALSE, "top", "none"),
                  extensions = "Buttons",
                  rownames = FALSE,
                  options = list(dom = 'Bfrtip',
                                 buttons = c('csv', 'excel'))) %>%
    DT::formatSignif(columns = c("LogFC", "PValue", "FDR"), digits = 2)
}

#' @title printCorTibble
#' @description Print the corTibble
#' @param corTibble A tibble of tidied correlations
#' @param BiocStyle Logical indicating whether the table is rendered in a BiocStyle vignette, where the per-column filter boxes misalign, so they are dropped in favour of the global search
#' @return A datatable
#' @import DT
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' data(airwayDGE)
#' corRes <- tidyr::crossing(x = names(airwayDGE), y = names(airwayDGE)) %>%
#'     dplyr::filter(x != y) %>%
#'     dplyr::mutate(name = paste(stringr::str_remove(x, " vs .*"), "vs",
#'                                stringr::str_remove(y, " vs .*"))) %>%
#'     prepareCor(airwayDGE, cd = SummarizedExperiment::colData(se),
#'                format = "long")
#' printCorTibble(corRes)
printCorTibble <- function(corTibble,
                           BiocStyle = FALSE){
  corTibble %>%
    DT::datatable(filter = ifelse(BiocStyle == FALSE, "top", "none"),
                  extensions = "Buttons",
                  options = list(dom = 'Blfrtip',
                                 buttons = c('csv', 'excel'))) %>%
    DT::formatSignif(columns = "Correlation",  digits = 2)
}

#' printPlots
#' @title Print ggplots
#' @description Create tabs of ggplots in R Markdown
#' @param .x plot
#' @param .y title
#' @param depth Number of # to use
#' @importFrom knitr knit_print
#' @export
#'
#' @return `NULL`; called for its side effect of printing a markdown header followed by the plot.
#' @examples
#' data(airwayGSEA)
#' airwayGSEA %>%
#'     GSEAplotList() %>%
#'     purrr::iwalk(printPlots)
printPlots <- function(.x, .y, depth = 2){
  cat("\n")
  cat(paste(rep("#", depth), collapse = ""), .y)
  cat("\n")
  knitr::knit_print(.x)
  cat("\n")
}

#' printHeatmaps
#' @title Print heatmaps
#' @inheritParams printPlots
#' @description Create tabs of heatmaps in R Markdown
#' @export
#'
#' @return `NULL`; called for its side effect of printing a markdown header followed by the heatmap.
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' data(airwayDGE)
#' list(DEGs = basicHeatmap(se, sigTibble = makeSigTibble(airwayDGE))) %>%
#'     purrr::iwalk(printHeatmaps)
printHeatmaps <- function(.x, .y, depth = 2){
  cat("\n")
  cat(paste(rep("#", depth), collapse = ""), .y)
  cat("\n")
  print(.x)
  cat("\n")
}
