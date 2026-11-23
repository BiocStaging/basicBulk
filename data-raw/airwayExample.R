#' @title Create example data
#' @description Build the example objects shipped with basicBulk from the
#'  \link[airway]{airway} dataset via the package's own workflow
#'  (`tidySE()` -> `augere.de::runVoom()` -> `makeDGEdata()`, and
#'  `augere.gsea::runContrast()` -> `tidyGSEAlist()`). Requires augere.de and
#'  augere.gsea; the objects are saved with `usethis::use_data()`. The DGEdata
#'  is subset to a compact gene set after testing so the shipped `.rda` files
#'  stay small. The SE itself is not shipped: examples rebuild it from
#'  `augere.de::loadExampleDataset()` with `tidySE()`.
library(basicBulk)
library(magrittr)

organism <- "Homo sapiens"

se0 <- augere.de::loadExampleDataset()
SummarizedExperiment::rowData(se0)$ID <- rownames(se0)

airwaySE <- se0 %>%
    basicBulk::tidySE(group = "cell", block = "dex")

comparisonList <- list(c("N61311", "N052611"),
                       c("N061011", "N052611"),
                       c("N080611", "N052611"))

voomList <- airwaySE %>%
    augere.de::runVoom(groups = "cell",
                       block = "dex",
                       comparisons = comparisonList,
                       assay = "counts",
                       output.dir = tempfile(),
                       save.results = FALSE) %>%
    magrittr::extract2("results")

airwayDGE <- airwaySE %>%
    basicBulk::makeDGEdata(voomList)

gene.sets <- basicBulk::getReactome(organism)

set.seed(100)
airwayGSEA <- comparisonList %>%
    purrr::set_names(sapply(., `[`, 1)) %>%
    purrr::map(~ airwaySE %>%
                   augere.gsea::runContrast(sets = gene.sets,
                                            groups = "cell",
                                            block = "dex",
                                            comparison = .x,
                                            methods = "camera",
                                            assay = "counts",
                                            output.dir = tempfile(),
                                            save.results = FALSE)) %>%
    basicBulk::tidyGSEAlist()

## Subset the DGEdata to a compact gene set: the most variable retained genes
## plus a few named genes used in examples, so the shipped data stays small.
named <- c("FKBP5", "DUSP1", "PER1", "RARRES1", "TMEM176B", "GBP3")

retained <- rownames(airwaySE)[SummarizedExperiment::rowData(airwaySE)$retained]
vars <- apply(SummarizedExperiment::assay(airwaySE, "logCPM")[retained, ], 1, var)
topVar <- names(sort(vars, decreasing = TRUE))[seq_len(min(3000, length(vars)))]

namedIDs <- rownames(airwaySE)[SummarizedExperiment::rowData(airwaySE)$symbol %in% named]
keepIDs <- union(topVar, namedIDs)

airwayDGE <- airwayDGE %>%
    purrr::map(~ dplyr::filter(.x, ID %in% keepIDs))

usethis::use_data(airwayDGE, airwayGSEA, overwrite = TRUE)
