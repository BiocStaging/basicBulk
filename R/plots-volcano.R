#' @title basicVolcano
#' @description A wrap of basic volcanos with DEG counts
#' @param DGEdata A list of DGE results
#' @param FDRcutoff Numeric for the FDR cut-off for DEGs, default is 0.05
#' @param logFCcutoff Numeric for the absolute Log2FC cut-off for DEGs, default is 1
#' @param lab TRUE adds labels for all significant genes, while a character vector of gene symbols will only label that selection
#' @param pointSize Numeric specifying the size of the points
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param labSize Numeric for font size of labels (in pt)
#' @param labFace Character specifying font face of gene labels (i.e bold, italic, bolditalic)
#' @param style Character specifying either a "basic" or "enhanced" style.
#' @return A ggplot2 wrap
#' @importFrom dplyr arrange select mutate recode_values bind_rows case_when
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayDGE)
#' basicVolcano(airwayDGE)
basicVolcano <- function(DGEdata,
                         FDRcutoff = 0.05,
                         logFCcutoff = 1,
                         lab = NULL,
                         pointSize = 0.1,
                         nrow = NULL,
                         ncol = NULL,
                         labSize = 9,
                         labFace = "italic",
                         style = c("basic", "enhanced")){

  style <- match.arg(style)

  counts <- DGEdata %>%
    basicBulk::makeSigTibble(FDRcutoff = FDRcutoff,
                           logFCcutoff = logFCcutoff) %>%
    basicBulk::makeSigCountTibble() %>%
    dplyr::arrange(contrast) %>%
    dplyr::select(-total) %>%
    tidyr::pivot_longer(cols = c(Down, Up),
                        names_to = "direction",
                        values_to = "count") %>%
    dplyr::mutate(x = dplyr::recode_values(direction,
                                           "Up" ~ Inf,
                                           "Down" ~ -Inf),
                  y = Inf,
                  color = dplyr::case_when(style == "basic" & direction == "Up" ~ "firebrick",
                                           style == "basic" & direction == "Down" ~ "mediumblue",
                                           style == "enhanced" ~ "firebrick"),
                  hjust = dplyr::recode_values(direction,
                                               "Up" ~ 1,
                                               "Down" ~ 0),
                  vjust = dplyr::recode_values(direction,
                                               "Up" ~ 1,
                                               "Down" ~ 1))

  plotData <- DGEdata %>%
    dplyr::bind_rows(.id = "contrast") %>%
    dplyr::mutate(contrast = forcats::as_factor(contrast)) %>%
    {if(style == "basic"){
      dplyr::mutate(., color = dplyr::case_when(abs(LogFC) > logFCcutoff & FDR < FDRcutoff & sign(LogFC) == 1 ~ "Up",
                                                abs(LogFC) > logFCcutoff & FDR < FDRcutoff & sign(LogFC) == -1 ~ "Down",
                                                .default = "NS") %>%
                      factor(levels = c("Up", "Down", "NS")))
    }else if(style == "enhanced"){
      dplyr::mutate(., color = dplyr::case_when(abs(LogFC) > logFCcutoff & FDR < FDRcutoff ~ "FDR & LogFC",
                                                abs(LogFC) <= logFCcutoff & FDR < FDRcutoff ~ "FDR",
                                                abs(LogFC) > logFCcutoff & FDR >= FDRcutoff ~ "LogFC",
                                                .default = "NS") %>%
                      factor(levels = c("FDR & LogFC", "FDR", "LogFC", "NS")))}}
  (plotData %>%
      ggplot2::ggplot(ggplot2::aes(x = LogFC,
                                   y = -log10(FDR),
                                   color = color)) +
      ggplot2::geom_point(shape = 19,
                          size = pointSize) +
      ggplot2::facet_wrap(~contrast, nrow = nrow, ncol = ncol) +
      ggplot2::geom_label(data = counts,
                          ggplot2::aes(x = x,
                                       y = y,
                                       color = I(color),
                                       hjust = hjust,
                                       vjust = vjust,
                                       label = count),
                          size = labSize/ggplot2::.pt,
                          fill = NA,
                          show.legend = FALSE,
                          label.padding = ggplot2::unit(0.1, "lines"),
                          linewidth = NA) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = expression("Log"["2"] (FC)),
                    y = expression("-Log"["10"](italic(FDR)))) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())) %>%
    {if(!is.null(lab)){
      . +
        ggrepel::geom_text_repel(data = if(isTRUE(lab)){
          plotData %>%
            dplyr::filter(abs(LogFC) > logFCcutoff & FDR < FDRcutoff)
        }else{
          plotData %>%
            dplyr::filter(symbol %in% lab)
        },
        ggplot2::aes(label = symbol),
        segment.color = "grey30",
        segment.size = 0.5,
        min.segment.length = 0,
        arrow = grid::arrow(length = ggplot2::unit(0.0025, "npc"),
                            ends = "first",
                            type = "closed"),
        size = labSize/ggplot2::.pt,
        fontface = labFace,
        na.rm = TRUE,
        direction = "both")
    }else .} %>%
    {if(style == "enhanced"){
      . +
        ggplot2::geom_hline(yintercept = -log10(FDRcutoff),
                            color = "grey30",
                            linetype = "dashed") +
        ggplot2::geom_vline(xintercept = c(-logFCcutoff, logFCcutoff),
                            color = "grey30",
                            linetype = "dashed") +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::scale_color_manual(values = c("NS" = "grey80",
                                               "LogFC" = "forestgreen",
                                               "FDR" = "mediumblue",
                                               "FDR & LogFC" = "firebrick"),
                                    labels = c("NS" = "Not Significant",
                                               "LogFC" = expression(Log[2] ~ "Fold Change"),
                                               "FDR" = "FDR",
                                               "FDR & LogFC" = expression("FDR & " ~ Log[2] ~ "Fold Change")),
                                    name = NULL,
                                    guide = ggplot2::guide_legend(override.aes = list(size = 4)))

    }else{
      . +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_color_manual(values = c("Up" = "firebrick",
                                               "Down" = "mediumblue",
                                               "NS" = "grey80"))}}
  }
