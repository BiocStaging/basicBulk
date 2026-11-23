#' @title basicMA
#' @description Create a grid of MA plots to examinine the LogFC (M) and average logCPM (A)
#' @inheritParams basicVolcano
#' @param pointShape Numeric  (or ".") specifying shape of points
#' @param pointSize Numeric specifying size of points
#' @return A ggplot2 wrap
#' @importFrom dplyr bind_rows mutate case_when
#' @importFrom forcats as_factor
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayDGE)
#' basicMA(airwayDGE)
basicMA <- function(DGEdata,
                    nrow = NULL,
                    ncol = NULL,
                    FDRcutoff = 0.05,
                    logFCcutoff = 1,
                    pointShape = 19,
                    pointSize = 0.1){
  DGEdata %>%
    dplyr::bind_rows(.id = "contrast") %>%
    dplyr::mutate(contrast = forcats::as_factor(contrast),
                  color = dplyr::case_when(abs(LogFC) > logFCcutoff & FDR < FDRcutoff & sign(LogFC) == 1 ~ "red",
                                           abs(LogFC) > logFCcutoff & FDR < FDRcutoff & sign(LogFC) == -1 ~ "blue",
                                           .default = "grey60")) %>%
    ggplot2::ggplot(ggplot2::aes(x = AveExpr,
                                 y = LogFC,
                                 color = I(color))) +
    ggplot2::geom_point(shape = pointShape,
                        alpha = 0.8,
                        size = pointSize) +
    ggplot2::facet_wrap(~contrast, nrow = nrow, ncol = ncol) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = "Average logCPM")
}

#' @title DEGbar
#' @description Bar plot of total DEGs
#' @param sigTabyl A tibble from \link{makeSigCountTibble} or a sigTibble run
#'  through \code{janitor::tabyl(contrast, direction)} and mutated with other
#'  variables of interest
#' @param x Character specifying the x-axis variable
#' @param addCounts Logical indicating whether to add text with totals
#' @param flip Logical indicating whether to flip the coordinates
#' @param linewidth Numeric specifying linewidth for black outline
#' @return A ggplot2
#' @import ggplot2
#' @import shadowtext
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate case_when arrange case_when
#' @importFrom forcats fct_rev
#' @importFrom ggpubr rotate_x_text
#' @importFrom rlang sym
#' @importFrom purrr pluck
#' @importFrom magrittr %>%
#' @importFrom MetBrewer met.brewer
#' @export
#'
#' @examples
#' data(airwayDGE)
#' airwayDGE %>%
#'     makeSigTibble() %>%
#'     makeSigCountTibble() %>%
#'     DEGbar()
DEGbar <- function(sigTabyl,
                   x = "contrast",
                   addCounts = TRUE,
                   flip = TRUE,
                   linewidth = 0.1){
  p1 <- sigTabyl %>%
    tidyr::pivot_longer(c(Up, Down),
                        names_to = "Direction",
                        values_to = "DEGs") %>%
    dplyr::mutate(DEGs = dplyr::case_when(Direction == "Down" ~ -DEGs,
                                          TRUE ~ DEGs),
                  Direction = Direction %>%
                    as.factor() %>%
                    forcats::fct_rev()) %>%
    ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x),
                                 y = DEGs,
                                 fill = Direction)) +
    ggplot2::geom_bar(stat = "identity",
                      position = "stack",
                      color = "black",
                      linewidth = linewidth) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 12),
                   legend.position = "bottom",
                   axis.title.x = ggplot2::element_blank()) +
    ggplot2::labs(y = "Number of DEGs") +
    ggpubr::rotate_x_text(45) +
    ggplot2::scale_y_continuous(labels = abs) +
    ggplot2::scale_fill_manual(values = MetBrewer::met.brewer("Cassatt1")[c(1, 8)] %>%
                                 as.character() %>%
                                 purrr::set_names(c("Up", "Down")))

  if(addCounts == TRUE){
    p1 <- p1 +
      ggplot2::aes(x = factor(!!rlang::sym(x),
                              levels = p1$data %>%
                                dplyr::arrange(total) %>%
                                purrr::pluck(x) %>%
                                unique())) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      shadowtext::geom_shadowtext(ggplot2::aes(label = total),
                                  color = dplyr::case_when(p1$data$Direction == "Down" ~ "Transparent",
                                                           .default = "White"),
                                  y = 0,
                                  size = 4,
                                  angle = dplyr::case_when(flip == TRUE ~ 0,
                                                           flip == FALSE ~ 90))
  }

 if(flip == TRUE){
   p1 <- p1 +
     ggplot2::coord_flip() +
     ggplot2::aes(x = !!rlang::sym(x) %>%
                    forcats::fct_rev())
 }

  p1
}

#' @title prepareCor
#' @description Calculate pairwise logFC correlations with \link[gg4way]{getCor}, returned in long or wide
#'  (\link[corrr]{corrr}) format
#' @param corTibble A dataframe used to specify correlations through the following columns: x, y, name
#' @param DGEdata A list of DGE results
#' @param cd colData from a \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param format Character specifying the format of the returned data frame, "wide" (the default) or "long"
#' @inheritParams purrr::pmap
#' @return A dataframe
#' @import gg4way
#' @import corrr
#' @importFrom purrr pmap set_names
#' @importFrom tibble enframe
#' @importFrom tidyr unnest separate_wider_delim pivot_wider
#' @importFrom dplyr mutate select across rename
#' @importFrom magrittr %>%
#' @references
#' Laufer BI, Friedman BA (2026). gg4way: 4way Plots of Differential Expression. R package.
#'  \doi{10.18129/B9.bioc.gg4way}
#'
#' Kuhn M, Jackson S, Cimentada J (2025). corrr: Correlations in R. R package.
#'  \doi{10.32614/CRAN.package.corrr}
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
#'     prepareCor(airwayDGE, cd = SummarizedExperiment::colData(se))
prepareCor <- function(corTibble,
                       DGEdata,
                       cd,
                       format = c("wide", "long"),
                       .progress = TRUE){

  format <- match.arg(format)

  corRes <- corTibble %>%
    purrr::pmap(~ with(list(...),
                       gg4way::gg4way(DGEdata,
                                      x = x,
                                      y = y,
                                      logFC = "LogFC",
                                      FDR = "FDR") %>%
                         gg4way::getCor()),
                .progress = .progress) %>%
    purrr::set_names(corTibble$name) %>%
    tibble::enframe() %>%
    tidyr::unnest(value) %>%
    tidyr::separate_wider_delim(name, delim = " vs ", names = c("x", "y")) %>%
    dplyr::mutate(Correlation = value %>%
                    as.numeric() %>%
                    round(digits = 2)) %>%
    dplyr::select(-value)

  corRes <- corRes %>%
    dplyr::mutate(dplyr::across(c(x, y), ~ factor(.x, levels = levels(cd$Group))))

  if(format == "long"){
    return(corRes)
  }else if(format == "wide"){
    corRes <- corRes %>%
      tidyr::pivot_wider(names_from = y,
                         values_from = Correlation) %>%
      dplyr::rename(term = x) %>%
      dplyr::select(term, dplyr::any_of(levels(cd$Group))) %>%
      dplyr::arrange(term)

    class(corRes) <- c("cor_df", class(corRes))

    return(corRes)
  }
}

#' @title corTriplot
#' @description Plot the lower triangle of logFC correlations from \link{prepareCor}
#' @param corRes cor_df from \link{prepareCor}
#' @return A ggplot
#' @import corrr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @references
#' Laufer BI, Friedman BA (2026). gg4way: 4way Plots of Differential Expression. R package.
#'  \doi{10.18129/B9.bioc.gg4way}
#'
#' Kuhn M, Jackson S, Cimentada J (2025). corrr: Correlations in R. R package.
#'  \doi{10.32614/CRAN.package.corrr}
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
#'     prepareCor(airwayDGE, cd = SummarizedExperiment::colData(se))
#' corTriplot(corRes)
corTriplot <- function(corRes){

  corRes %>%
    corrr::shave() %>%
    corrr::rplot(shape = 21) +
    ggplot2::aes(fill = r,
                 color = "black") +
    ggplot2::scale_color_manual(values = c("black" = "black"), guide = "none") +
    ggplot2::scale_fill_gradient2(guide = ggplot2::guide_colourbar(title = "r",
                                                                   theme = ggplot2::theme(legend.title = ggplot2::element_text(hjust = 0.1),
                                                                                          legend.title.position = "top",
                                                                                          legend.key.width = grid::unit(0.03, "npc"),
                                                                                          legend.key.height = grid::unit(0.2, "npc"),
                                                                                          legend.frame = ggplot2::element_rect(colour = "black",
                                                                                                                               linewidth = 0.1),
                                                                                          legend.ticks = ggplot2::element_line(colour = "black",
                                                                                                                               linewidth = 0.1))),
                                  low = "mediumblue",
                                  mid = "white",
                                  high = "firebrick",
                                  midpoint = 0) +
    ggpubr::rotate_x_text(45)
}

#' @title barStyle
#' @description Makes plots less baring
#' @return An improved ggplot2
#' @import ggplot2
#' @export
#'
#' @examples
#' data(airwayDGE)
#' airwayDGE %>%
#'     makeSigTibble() %>%
#'     makeSigCountTibble() %>%
#'     ggplot2::ggplot(ggplot2::aes(x = contrast, y = total)) +
#'     barStyle()
barStyle <- function(){
  list(ggplot2::geom_bar(stat = "identity",
                         position = ggplot2::position_dodge(),
                         color = "black",
                         linewidth = 0.1),
       ggplot2::coord_flip(),
       ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)),
       ggplot2::theme_bw(),
       ggplot2::theme(legend.position = "none",
                      axis.title.y = ggplot2::element_blank()))
}
