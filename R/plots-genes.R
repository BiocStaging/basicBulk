#' @title boxPlot
#' @description Makes a \link[ggpubr]{ggboxplot}
#' @param data The data
#' @param x Column for the x-axis
#' @param y Column for the y-axis
#' @param color Column for the color
#' @param size Numeric specifying size of the dots
#' @param plotly Logical indicating whether to add the text aes for \link[plotly]{ggplotly}
#' @param ... Additional arguments for \link[ggpubr]{ggboxplot}
#' @return A stylish ggppubr boxplot
#' @import ggpubr
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' se %>%
#'     tibble::as_tibble() %>%
#'     dplyr::filter(symbol == "FKBP5") %>%
#'     boxPlot(x = "Group", y = "logCPM", color = "Group")
boxPlot <- function(data,
                    x,
                    y,
                    color,
                    size = 1.25,
                    plotly = FALSE,
                    ...){

  p1 <- data %>%
    ggpubr::ggboxplot(x = x,
                      y = y,
                      color = color,
                      outlier.shape = NA,
                      add = "jitter",
                      linewidth = 0.3,
                      width = 0.4,
                      alpha = 0.2,
                      add.params = list(size = size, alpha = 0.4),
                      ...) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right",
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  if(plotly == TRUE){
    p1 <- p1 +
      ggplot2::aes(text = Sample)
  }else if(plotly == FALSE){
    p1 <- p1 +
      ggsci::scale_color_d3("category20") +
      ggsci::scale_fill_d3("category20")
  }

  return(p1)
}

#' @title dotPlot
#' @description Makes a \link[ggpubr]{ggstripchart}
#' @inheritParams ggpubr::ggstripchart
#' @inheritDotParams ggpubr::ggstripchart -position -fill -ggtheme
#' @param plotly Logical indicating whether to add the text aes for \link[plotly]{ggplotly}
#' @param alpha Numeric specifying the alpha
#' @return A stylish ggppubr strip chart
#' @import ggpubr
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' se %>%
#'     tibble::as_tibble() %>%
#'     dplyr::filter(symbol == "FKBP5") %>%
#'     dotPlot(x = "Group", y = "logCPM", color = "Group")
dotPlot <- function(data,
                    x,
                    y,
                    color,
                    shape = 21,
                    size = 1.25,
                    alpha = 0.6,
                    jitter = 0,
                    plotly = FALSE,
                    ...){

  p1 <- data %>%
    ggpubr::ggstripchart(x = x,
                         y = y,
                         fill = color,
                         color = ifelse(shape == 21, "black", color),
                         shape = shape,
                         size = size,
                         alpha = alpha,
                         position = ggplot2::position_jitter(width = jitter,
                                                             height = 0,
                                                             seed = 123),
                         ggtheme = ggplot2::theme_bw(),
                         ...) +
    ggplot2::theme(legend.position = "right",
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  if(plotly == TRUE){
    p1 <- p1 +
      ggplot2::aes(text = Sample)
  }else if(plotly == FALSE){
    if(shape == 21){
      p1 <- p1 +
        ggsci::scale_fill_d3("category20")
    }else if(shape != 21){
      p1 <- p1 +
        ggsci::scale_color_d3("category20")
    }
  }

  return(p1)
}

#' @title genePlot
#' @description Faceted plots for genes of interest
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object, or a tibble of
#'  gene set scores from \link{joinFeatures} with \code{score = TRUE} stacked with the gene set names in a
#'  \code{symbol} column
#' @param genes A character vector specifying the genes (or gene sets) to plot
#' @param DGEdata A list of DGE results, which requires a
#'  \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} \code{se}
#' @param assay Assay from se to use
#' @param addStats Logical indicating whether to add significance annotations from \code{DGEdata}
#' @param ref.group Character specifying reference group for statistics
#' @param logFCcutoff Numeric for the absolute Log2FC cut-off for DEGs, default is 1
#' @param plotType A character specifying a "box" or "dot" plot
#' @param size Numeric specifying size of the dots
#' @param jitter Numeric specifying jitter of the dots (for dot plot only)
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param nrow_legend Number of rows for the legend
#' @param tip.length Numeric specifying length of bracket tips
#' @param ggprism Logical indicating whether a prism theme should be applied
#' @param nudge_y Numeric specifying the multiplier to increase spacing between statistics results and max y value;
#'  the default is 0.1 with a \code{ref.group} (its asterisks start at the max y value) and 0 without (the brackets
#'  are already spaced above it)
#' @inheritDotParams ggprism::add_pvalue
#' @return A ggplot
#' @import ggplot2
#' @import ggprism
#' @import rstatix
#' @importFrom tibble as_tibble
#' @importFrom purrr imap list_rbind pmap
#' @importFrom dplyr filter mutate select distinct group_by if_else inner_join
#' @importFrom tidyr separate_wider_delim
#' @importFrom stats reformulate
#' @importFrom rlang sym
#' @importFrom cli cli_warn
#' @importFrom magrittr %>% %T>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' genePlot(se, genes = c("RARRES1", "TMEM176B", "GBP3"))
#'
#' ## Gene set scores
#' gene.sets <- getReactome("Homo sapiens")
#' gene.sets[c("REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL",
#'             "REACTOME_MUSCLE_CONTRACTION")] %>%
#'     purrr::map(~ joinFeatures(se, features = .x, score = TRUE)) %>%
#'     purrr::list_rbind(names_to = "symbol") %>%
#'     genePlot(genes = unique(.$symbol))
genePlot <- function(se,
                     genes,
                     DGEdata = NULL,
                     assay = "logCPM",
                     addStats = !is.null(DGEdata),
                     ref.group = NULL,
                     logFCcutoff = 1,
                     plotType = c("dot", "box"),
                     size = 2,
                     jitter = 0.1,
                     nrow = NULL,
                     ncol = NULL,
                     nrow_legend = NULL,
                     tip.length = 0,
                     ggprism = FALSE,
                     nudge_y = ifelse(is.null(ref.group), 0, 0.1),
                     ...){

  plotType <- match.arg(plotType)

  genes <- genes %>%
    forcats::as_factor()

  geneData <- se %>%
    tibble::as_tibble() %T>%
    {stopifnot(assay %in% colnames(.))} %>%
    dplyr::filter(symbol %in% !!genes,
                  !is.na(!!rlang::sym(assay))) %>%
    dplyr::mutate(symbol = factor(symbol, levels = genes))

  stopifnot("None of the genes have values to plot" = nrow(geneData) > 0)

  if(any(genes %ni% geneData$symbol)){
    cli::cli_warn("No {assay} values to plot for {setdiff(genes, geneData$symbol)}")
  }

  if(plotType == "box"){
    genePlots <- geneData %>%
      basicBulk::boxPlot(x = "Group",
                       y = assay,
                       color = "Group",
                       size = size)
  }else if(plotType == "dot"){
    genePlots <- geneData %>%
      basicBulk::dotPlot(x = "Group",
                       y = assay,
                       color = "Group",
                       size = size,
                       jitter = jitter)
  }

  genePlots <- genePlots +
    ggplot2::facet_wrap(~symbol, scales = "free", nrow = nrow, ncol = ncol) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = nrow_legend, theme = ggplot2::theme(legend.byrow = TRUE)),
                    fill = ggplot2::guide_legend(nrow = nrow_legend, theme = ggplot2::theme(legend.byrow = TRUE)))

  if(ggprism == TRUE){

    genePlots <- genePlots +
      ggprism::theme_prism(base_line_size = 0.5, base_size = 9, axis_text_angle = 45) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 9, face = "bold"),
                     strip.text = ggplot2::element_text(size = 9, face = "bold.italic"))
  }

  genePlots <- genePlots +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x =  ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  if(addStats == TRUE){

    stats <- DGEdata %>%
      purrr::imap(~ .x %>%
                    dplyr::filter(symbol %in% geneData$symbol) %>%
                    dplyr::mutate(contrast = .y)) %>%
      purrr::list_rbind() %>%
      dplyr::select(symbol, contrast, p.adj = FDR, LogFC) %>%
      tidyr::separate_wider_delim(contrast,
                                  delim = " vs ",
                                  names = c("group1", "group2"),
                                  too_few = "align_start") %>%
      dplyr::filter(p.adj < 0.05 & abs(LogFC) > logFCcutoff,
                    group1 %in% geneData$Group,
                    group2 %in% geneData$Group) %>%
      {if(is.null(ref.group)){
        .
      }else{
        dplyr::filter(., group1 == ref.group | group2 == ref.group) %>%
          dplyr::mutate(group2 = dplyr::if_else(group1 == ref.group, group2, group1),
                        group1 = ref.group)
      }} %>%
      dplyr::mutate(p.adj.signif = "*")

    signif <- geneData %>%
      dplyr::group_by(symbol) %>%
      {if(nrow(stats) == 0){
        stats
      }else{
        rstatix::get_y_position(.,
                                formula = stats::reformulate("Group", assay),
                                comparisons = stats %>%
                                  dplyr::distinct(group1, group2) %>%
                                  purrr::pmap(c),
                                step.increase = dplyr::case_when(is.null(ref.group) ~ 0.12,
                                                                 TRUE ~ 0),
                                scales = "free") %>%
          dplyr::inner_join(stats, by = c("symbol", "group1", "group2")) %>%
          dplyr::mutate(symbol = forcats::as_factor(symbol),
                        y.position = y.position + abs(y.position) * nudge_y)
      }}

    genePlots <- genePlots %>%
      {if(nrow(signif) == 0){
        .
      }else{
        . +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
          ggprism::add_pvalue(signif, label = "p.adj.signif", label.size = 5, tip.length = tip.length,
                              bracket.size = dplyr::case_when(ggprism == TRUE ~ 0.6,
                                                              TRUE ~ 0.3),
                              remove.bracket = dplyr::case_when(is.null(ref.group) ~ FALSE,
                                                                TRUE ~ TRUE),
                              ...)
      }} +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 9, face = "italic"))
  }

  return(genePlots)
}

#' @title geneHistogram
#' @description Histograms of mean gene expression for each group
#' @inheritParams genePlot
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param geneSymbol An optional character vector specifying the gene to draw a mean expression line for
#' @param palette Character specifying palette name for \link[ggsci]{scale_color_d3}
#' @return A ggplot
#' @import ggplot2
#' @importFrom SummarizedExperiment assays assay
#' @importFrom dplyr group_by summarise ungroup filter
#' @importFrom rlang sym
#' @importFrom ggsci scale_fill_d3
#' @importFrom ggpubr rotate_x_text
#' @importFrom cli cli_warn
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' geneHistogram(se, geneSymbol = "FKBP5")
geneHistogram <- function(se,
                          geneSymbol = NULL,
                          assay = "logCPM",
                          palette = "category20"){

  stopifnot(assay %in% names(SummarizedExperiment::assays(se)))

  p1 <- se %>%
    dplyr::group_by(Group, .feature) %>%
    dplyr::summarise(mean = mean(!!rlang::sym(assay), na.rm = TRUE)) %>%
    dplyr::filter(!is.na(mean)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(ggplot2::aes(x = mean,
                                 fill = Group)) +
    ggplot2::geom_bar(color = "black",
                      linewidth = 0.05) +
    ggplot2::scale_x_binned(breaks = seq(floor(min(SummarizedExperiment::assay(se, assay), na.rm = TRUE)),
                                         ceiling(max(SummarizedExperiment::assay(se, assay), na.rm = TRUE)),
                                         by = 1)) +
    ggplot2::facet_wrap(~Group) +
    ggsci::scale_fill_d3(palette = palette,
                         alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_y_continuous(expand = c(0,0,0.05,0)) +
    ggpubr::rotate_x_text() +
    ggplot2::labs(x = paste("Mean", assay),
                  y = "Count")

  if(!is.null(geneSymbol)){

    stats <- se %>%
      dplyr::filter(symbol == geneSymbol) %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(mean = mean(!!rlang::sym(assay), na.rm = TRUE)) %>%
      dplyr::filter(!is.na(mean)) %>%
      dplyr::ungroup()

    if(nrow(stats) == 0){
      cli::cli_warn("No {assay} values to plot for {geneSymbol}")
    }

    p1 <- p1 +
      ggplot2::geom_vline(data = stats,
                          ggplot2::aes(xintercept = mean,
                                       color = "Group"),
                          color = "firebrick",
                          linetype = "dashed")
  }

  return(p1)

}
