#' @title runPCA
#' @description Wrapper for generating a PCA with PCAtools via \link[PCAtools]{pca}
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param assay Assay from se to use
#' @param removeVar Argument to remove lower percent of variables for the \link[PCAtools]{pca}
#' @import PCAtools
#' @import SummarizedExperiment
#' @importFrom stats na.omit
#' @importFrom magrittr %>%
#' @export
#'
#' @return A \link[PCAtools]{pca} object.
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' runPCA(se)
runPCA <- function(se,
                   assay = "logCPM",
                   removeVar = 0.1){
  se %>%
    SummarizedExperiment::assay(assay) %>%
    stats::na.omit() %>%
    PCAtools::pca(metadata = SummarizedExperiment::colData(se),
                  removeVar = removeVar)
}

#' @title PCAplot
#' @description Wrapper for generating a PCA biplot with PCAtools via \link[PCAtools]{pca} and \link[PCAtools]{biplot}
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param PCA A pca object from \link[PCAtools]{pca}
#' @param assay Assay from se to use
#' @param colby colData variable to color points by
#' @param lab colData variable to label points by, or \code{NULL} for no labels
#' @param alpha Numeric specifying alpha for points
#' @inheritParams PCAtools::biplot
#' @inheritDotParams PCAtools::biplot
#' @inheritParams ggsci::scale_color_d3
#' @return A ggplot2
#' @import PCAtools
#' @importFrom ggsci scale_color_d3
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' PCAplot(se)
PCAplot <- function(se,
                    PCA = NULL,
                    assay = "logCPM",
                    colby = "Group",
                    lab = "Sample",
                    pointSize = 6,
                    labSize = 4,
                    alpha = 1,
                    legendPosition = "bottom",
                    legendLabSize = 16,
                    legendIconSize = 6,
                    palette = "category20",
                    x = "PC1",
                    y = "PC2",
                    ...){

  if(is.null(PCA)){
    PCA <- se %>%
      basicBulk::runPCA(assay = assay)
  }

  PCA %>%
    PCAtools::biplot(x = x,
                     y = y,
                     colby = colby,
                     lab = if(!is.null(lab)){as.character(PCA$metadata[[lab]])}else{NULL},
                     pointSize = pointSize,
                     labSize = labSize,
                     legendPosition = legendPosition,
                     legendLabSize = legendLabSize,
                     legendIconSize = legendIconSize,
                     xlab = paste0(x, " (", round(PCA$variance[x], digits = 2), "%)"),
                     ylab = paste0(y, " (", round(PCA$variance[y], digits = 2), "%)"),
                     ...) +
    ggsci::scale_color_d3(palette, alpha = alpha) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = legendPosition,
                   legend.box = "vertical",
                   legend.margin = ggplot2::margin()) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = legendIconSize),
                                                  order = 0),
                    color = ggplot2::guide_legend(override.aes = list(size = legendIconSize),
                                                  order = 1))
}

#' @title plotPCAloadings
#' @description Plot of results from PCA loadings from \link[PCAtools]{getLoadings}
#' @inheritParams PCAplot
#' @import PCAtools
#' @import ggpubr
#' @import SummarizedExperiment
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr select left_join mutate case_when group_by arrange desc slice ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom MetBrewer met.brewer
#' @export
#'
#' @return A ggplot of the PC1 and PC2 loadings.
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' plotPCAloadings(runPCA(se), se)
plotPCAloadings <- function(PCA,
                            se){

  PCA %>%
    PCAtools::getLoadings() %>%
    tibble::rownames_to_column("ID") %>%
    tibble::as_tibble() %>%
    dplyr::select(ID, PC1, PC2) %>%
    tidyr::pivot_longer(!ID,
                        names_to = "PC",
                        values_to = "Loading") %>%
    dplyr::left_join(se %>%
                       rowData() %>%
                       tibble::as_tibble() %>%
                       dplyr::select(ID, Symbol = symbol),
                     by = "ID") %>%
    dplyr::mutate(Direction = dplyr::case_when(Loading > 0 ~ "Up",
                                               Loading < 0 ~ "Down")) %>%
    dplyr::group_by(PC, Direction) %>%
    dplyr::arrange(dplyr::desc(abs(Loading))) %>%
    dplyr::slice(1:5) %>%
    dplyr::ungroup() %>%
    ggpubr::ggdotchart(x = "Symbol",
                       y = "Loading",
                       color = "Direction",
                       sorting = "descending",
                       add = "segments",
                       add.params = list(color = "grey40", size = 1),
                       dot.size = 2.5,
                       rotate = TRUE,
                       ggtheme = ggplot2::theme_bw()) +
    ggplot2::facet_wrap(~PC, ncol = 1, scales = "free") +
    ggplot2::geom_hline(yintercept = 0, linetype = 1, color = "grey40") +
    ggpubr::rotate_x_text(angle = 45) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
    ggplot2::scale_color_manual(values = MetBrewer::met.brewer("Cassatt1")[c(1, 8)] %>%
                                  as.character() %>%
                                  purrr::set_names(c("Up", "Down")),
                                guide = "none") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
}

#' @title PCAgrid
#' @description A grid with a PCA biplot and top loadings for PCs 1 and 2
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param assay Assay from se to use
#' @param removeVar Argument to remove lower percent of variables for \link[PCAtools]{pca}
#' @param pal Named character vector of colors for the PCA biplot
#' @param ... Additional arguments for PCAplot()
#' @import PCAtools
#' @import SummarizedExperiment
#' @import ggplot2
#' @importFrom purrr imap
#' @importFrom cowplot plot_grid
#' @importFrom magrittr %>%
#' @references
#' Blighe K, Andrews J, Lun A (2026). PCAtools: Everything Principal Components Analysis. R package.
#'  \doi{10.18129/B9.bioc.PCAtools}
#' @export
#'
#' @return A ggplot (a cowplot grid of the PCA biplot and loadings).
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' PCAgrid(se, colby = "dex", shape = "cell")
PCAgrid <- function(se,
                    assay = "logCPM",
                    removeVar = 0.1,
                    pal = NULL,
                    ...){

  PCA <- se %>%
    basicBulk::runPCA(assay = assay,
                    removeVar = removeVar)

  list("PCA Biplot" = basicBulk::PCAplot(PCA = PCA,
                                       ...) %>%
         {if(!is.null(pal)){
           . + ggplot2::scale_color_manual(values = pal)
         }else{
           .
         }},
       "PCA Loadings" = basicBulk::plotPCAloadings(PCA,
                                                 se)) %>%
    purrr::imap(~ .x +
                  ggplot2::labs(title = .y)) %>%
    cowplot::plot_grid(plotlist = .,
                       rel_widths = c(2.5,1))
}
