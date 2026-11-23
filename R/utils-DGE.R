#' @title tidyVoom
#' @description Tidy a single contrast's DGE stats into a tibble for downstream functions
#' @param contrast A \link[S4Vectors:DataFrame-class]{DataFrame} of voom results for one contrast from \link[augere.de]{runVoom}
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @return A tibble
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr left_join select everything any_of mutate
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' se %>%
#'     augere.de::runVoom(groups = "cell",
#'                        block = "dex",
#'                        comparisons = list(c("N61311", "N052611")),
#'                        assay = "counts",
#'                        output.dir = tempfile(),
#'                        save.results = FALSE) %>%
#'     magrittr::extract2("results") %>%
#'     magrittr::extract2(1) %>%
#'     tidyVoom(se)
tidyVoom <- function(contrast, se){
  contrast %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    stats::na.omit() %>%
    dplyr::left_join(se %>%
                       SummarizedExperiment::rowData() %>%
                       tibble::as_tibble(),
                     by = "ID") %>%
    tibble::as_tibble() %>%
    {if("desc" %in% names(.)){
      dplyr::mutate(., desc = stringr::str_remove(desc, " \\[.*?\\]"))
    }else{
      .
    }} %>%
    dplyr::select(dplyr::any_of(c("ID", "symbol", "desc")),
                  dplyr::everything())
}

#' @title makeDGEdata
#' @description Create a named list of DGE stats tibbles for downstream functions
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param voomList A list of \link[S4Vectors:DataFrame-class]{DataFrame}s of voom results from \link[augere.de]{runVoom}
#' @return A list of tibbles
#' @importFrom purrr map set_names
#' @importFrom stringr str_remove str_remove_all str_replace str_trim
#' @importFrom magrittr %>%
#' @references
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses
#'  for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research}, 43(7), e47. \doi{10.1093/nar/gkv007}
#'
#' Law CW, Chen Y, Shi W, Smyth GK (2014). voom: precision weights unlock linear model analysis tools for RNA-seq
#'  read counts. \emph{Genome Biology}, 15(2), R29. \doi{10.1186/gb-2014-15-2-r29}
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' voomList <- se %>%
#'     augere.de::runVoom(groups = "cell",
#'                        block = "dex",
#'                        comparisons = list(c("N61311", "N052611")),
#'                        assay = "counts",
#'                        output.dir = tempfile(),
#'                        save.results = FALSE) %>%
#'     magrittr::extract2("results")
#' makeDGEdata(se, voomList)
makeDGEdata <- function(se,
                        voomList){

  voomList %>%
    purrr::map(~ basicBulk::tidyVoom(.x, se)) %>%
    purrr::set_names(names(.) %>%
                       stringr::str_remove("Increase in ") %>%
                       stringr::str_replace("` over `", "` vs `") %>%
                       stringr::str_remove_all("`") %>%
                       stringr::str_trim())
}

#' @title makeSigTibble
#' @description Make the sigTibble
#' @param DGEdata A list of tibbles with tidied voom results
#' @param FDRcutoff Numeric for the FDR cut-off for DEGs, default is 0.05
#' @param logFCcutoff Numeric for the absolute Log2FC cut-off for DEGs, default is 1
#' @return A tibble of tidied significant DGE results, with \code{contrast} a factor whose levels are all the
#'  \code{DGEdata} names
#' @importFrom dplyr filter mutate case_when
#' @importFrom purrr map list_rbind
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayDGE)
#' makeSigTibble(airwayDGE)
makeSigTibble <- function(DGEdata,
                          FDRcutoff = 0.05,
                          logFCcutoff = 1){
  DGEdata %>%
    purrr::map(~ dplyr::filter(., FDR < FDRcutoff & abs(LogFC) > logFCcutoff)) %>%
    purrr::list_rbind(names_to = "contrast") %>%
    dplyr::mutate(contrast = factor(contrast, levels = sort(unique(names(DGEdata)), method = "radix")),
                  direction = dplyr::case_when(LogFC > 0 ~ "Up",
                                               LogFC < 0 ~ "Down"))
}

#' @title makeSigCountTibble
#' @description Prepare the sigTibble for \link{DEGbar}
#' @param sigTibble A tibble from \link{makeSigTibble}. Every \code{contrast} level gets a row, with zero counts for
#'  contrasts without DEGs; to count only some contrasts, subset the \code{DGEdata} before \link{makeSigTibble}
#' @return A tibble
#' @importFrom janitor tabyl
#' @importFrom dplyr mutate
#' @importFrom forcats fct_infreq fct_rev
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayDGE)
#' airwayDGE %>%
#'     makeSigTibble() %>%
#'     makeSigCountTibble()
makeSigCountTibble <- function(sigTibble){
  if(nrow(sigTibble) == 0){
    return(tibble::tibble(contrast = factor(levels(sigTibble$contrast),
                                            levels = rev(levels(sigTibble$contrast))),
                          Up = 0,
                          Down = 0,
                          total = 0))
  }

  sigTibble %>%
    dplyr::mutate(direction = factor(direction, levels = c("Up", "Down"))) %>%
    janitor::tabyl(contrast, direction) %>%
    dplyr::mutate(total = Up + Down,
                  contrast = contrast %>%
                    forcats::fct_infreq(w = total) %>%
                    forcats::fct_rev())
}

#' @title joinFeatures
#' @description Join features in the tidy transcriptomics style
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param features A character vector of features
#' @param assay Assay from se to use
#' @param symbol Logical indicating whether the feature character vector contains gene symbols,
#'  FALSE indicates it contains Ensembl IDs
#' @param score Logical indicating whether a gene set score should be returned
#' @param long Logical indicating whether long format should be returned (FALSE for wide matrix style)
#' @return A tibble
#' @importFrom tibble tibble rownames_to_column as_tibble
#' @importFrom dplyr summarise across as_tibble left_join
#' @importFrom tidyselect any_of
#' @importFrom tidyr pivot_longer
#' @import SummarizedExperiment
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' joinFeatures(se, c("RARRES1", "TMEM176B", "GBP3"), symbol = TRUE)
joinFeatures <- function(se,
                         features,
                         symbol = FALSE,
                         score = FALSE,
                         long = TRUE,
                         assay = "logCPM"){

  type <- dplyr::case_when(symbol == FALSE ~ "ID",
                           TRUE ~ "symbol")

  featureTibble <- tibble::tibble(!!type := features) %>%
    dplyr::left_join(se %>%
                       SummarizedExperiment::assay(assay)  %>%
                       as.data.frame() %>%
                       tibble::rownames_to_column("ID") %>%
                       dplyr::left_join(se %>%
                                          rowData() %>%
                                          tibble::as_tibble() %>%
                                          dplyr::select(ID, symbol),
                                        by = "ID"),
                     by = type)

  if(score == TRUE){
    featureTibble <- featureTibble %>%
      dplyr::summarise(dplyr::across(.cols = !tidyselect::any_of(c("ID", "symbol")),
                                     .fns = ~ mean(.x, na.rm = TRUE)))
  }

  if(long == TRUE){
    featureTibble <- featureTibble %>%
      tidyr::pivot_longer(cols = !tidyselect::any_of(c("ID", "symbol")),
                          names_to = "Sample",
                          values_to = assay) %>%
      dplyr::left_join(se %>%
                         SummarizedExperiment::colData() %>%
                         tibble::as_tibble(),
                       by = "Sample")
  }

  return(featureTibble)

}
