#' @title basicHeatmap
#' @description Create a basic heatmap of DEGs or pathways.
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#'  with "Group" in the colData.
#' @param sigTibble A tibble of significant genes.
#' @param assay Character specifying the name of the assay from the se.
#' @param contrast Character specifying the name of a contrast (\code{"X vs Y"}) to subset by if not NULL: keeps the
#'  samples of its two groups and, for a sigTibble of DEGs, only that contrast's genes.
#' @param pathways Logical indicating whether sigTibble is a pathwayTibble.
#' @param geneSet A gene set of Ensembl IDs or symbols to use instead of the sigTibble; with
#'  \code{cluster_rows = FALSE} the rows keep its order.
#' @param draw Logical indicating whether to call \link[ComplexHeatmap]{draw}.
#' @param heatmap_legend_side Character specifying the side for the heatmap legend.
#' @param annotation_legend_side Character specifying the side for the annotation legend.
#' @param ...  Additional arguments passed onto \link[ComplexHeatmap]{Heatmap}.
#' @inheritParams scaleRows
#' @inheritParams colorRamp3
#' @inheritParams ComplexHeatmap::Heatmap
#' @inheritParams ComplexHeatmap::`draw,HeatmapList-method`
#' @return A \link[ComplexHeatmap:HeatmapList-class]{HeatmapList} object.
#' @import ComplexHeatmap
#' @importFrom SummarizedExperiment assays assay colData rowData
#' @importFrom ggsci pal_d3
#' @importFrom purrr set_names
#' @importFrom dplyr filter pull left_join select arrange
#' @importFrom stringr str_split
#' @importFrom tibble tibble
#' @importFrom grid gpar unit
#' @importFrom magrittr extract %>%
#' @references
#' Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic
#'  data. \emph{Bioinformatics}. \doi{10.1093/bioinformatics/btw313}
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' data(airwayDGE)
#' basicHeatmap(se, sigTibble = makeSigTibble(airwayDGE))
basicHeatmap <- function(se,
                         sigTibble = NULL,
                         geneSet = NULL,
                         assay = "logCPM",
                         contrast = NULL,
                         pathways = FALSE,
                         cluster_columns = FALSE,
                         cluster_rows = TRUE,
                         column_title = NULL,
                         show_row_names = FALSE,
                         show_column_names = FALSE,
                         heatmap_legend_side = "bottom",
                         annotation_legend_side = "right",
                         verbose = FALSE,
                         legendLimits = NULL,
                         top_annotation = NULL,
                         right_annotation = NULL,
                         row_split = NULL,
                         draw = TRUE,
                         ...){

  stopifnot(assay %in% names(SummarizedExperiment::assays(se)))
  stopifnot(is.logical(pathways))

  if(is.null(sigTibble) & is.null(geneSet)){
    stop("No genes specified")
  }

  if(!is.null(sigTibble) & !is.null(geneSet)){
    stop("Only one of sigTibble or geneSet can be specified")
  }

  if(is.null(top_annotation)){
    col <- list("Group" = ggsci::pal_d3()(n = length(levels(se$Group))) %>%
                  purrr::set_names(se$Group %>% levels()))
  }

  if(!is.null(contrast)){
    if(!is.null(sigTibble) & pathways == FALSE){
      sigTibble <- sigTibble %>%
        dplyr::filter(contrast == !!contrast)
    }

    samples <- colnames(se)
    se <- se %>%
      .[, .$Group %in% unlist(stringr::str_split(contrast, " vs "))]
  }

  if(!is.null(sigTibble)){
    genes <- sigTibble %>%
      dplyr::pull(ID) %>%
      unique() %>%
      .[. %in% rownames(SummarizedExperiment::assay(se, assay))]
  }else if(!is.null(geneSet)){
    geneSet <- geneSet[!is.na(geneSet)]
    if(all(stringr::str_detect(geneSet, "^ENS"))){
      genes <- geneSet %>%
        .[. %in% rownames(SummarizedExperiment::assay(se, assay))]
    }else{
      genes <- se %>%
        SummarizedExperiment::rowData() %>%
        as.data.frame() %>%
        dplyr::filter(symbol %in% geneSet) %>%
        {if(isFALSE(cluster_rows)){
          dplyr::arrange(., match(symbol, geneSet))
        }else{
          .
        }} %>%
        rownames()
    }
  }

  if(length(genes) < 2){
    return(NULL)
  }

  heatMatrix <- genes %>%
    SummarizedExperiment::assay(se, assay)[.,] %>%
    as.matrix() %>%
    basicBulk::scaleRows(verbose = verbose) %>%
    .[rowSums(is.na(.)) < ncol(.), , drop = FALSE]

  if(nrow(heatMatrix) < 2){
    return(NULL)
  }

  if(is.null(top_annotation)){
    top_annotation <- ComplexHeatmap::HeatmapAnnotation("Group" = se$Group,
                                                        col = col,
                                                        annotation_legend_param = list(title = NULL,
                                                                                       labels_gp = grid::gpar(fontsize = 9)),
                                                        border = TRUE,
                                                        show_annotation_name = FALSE)

  }

  if(!is.null(contrast) & stats::nobs(top_annotation) %ni% ncol(se)){
    top_annotation <- top_annotation %>%
      subset(., samples %in% colnames(se))
  }

  if(pathways == TRUE){
    if(nlevels(sigTibble$Pathway) > 1){
      stopifnot(rlang::is_installed("pals"))

      geneGroups <- heatMatrix[rownames(heatMatrix) %in% genes,] %>%
        rownames() %>%
        tibble::tibble(ID = .) %>%
        dplyr::left_join(sigTibble)

      right_annotation <- ComplexHeatmap::rowAnnotation("Pathway" = geneGroups %>%
                                                          dplyr::pull(Pathway),
                                                        col = list("Pathway" = pals::glasbey(n = nlevels(sigTibble$Pathway)) %>%
                                                                     purrr::set_names(levels(sigTibble$Pathway)) %>%
                                                                     magrittr::extract(levels(geneGroups$Pathway))),
                                                        border = TRUE,
                                                        show_annotation_name = FALSE,
                                                        annotation_legend_param = list(border = TRUE,
                                                                                       labels_gp = grid::gpar(fontsize = 10)))

      row_split <- geneGroups$Pathway
      }}

  heatmap <- heatMatrix %>%
    ComplexHeatmap::Heatmap(top_annotation = top_annotation,
                            right_annotation = right_annotation,
                            row_labels = if(show_row_names == TRUE){
                              tibble::tibble("ID" = rownames(heatMatrix)) %>%
                                dplyr::left_join(se %>%
                                                   SummarizedExperiment::rowData() %>%
                                                   tibble::as_tibble()) %>%
                                dplyr::pull(symbol)
                            }else{rownames(heatMatrix)},
                            column_order = if(cluster_columns == FALSE){
                              se %>%
                                SummarizedExperiment::colData() %>%
                                as.data.frame() %>%
                                tibble::rownames_to_column() %>%
                                tibble::as_tibble() %>%
                                dplyr::select(rowname,
                                              Group) %>%
                                dplyr::arrange(Group) %>%
                                dplyr::pull(rowname)
                            }else{NULL},
                            col = heatMatrix %>%
                              basicBulk::colorRamp3(legendLimits = legendLimits),
                            name = "Z-score",
                            column_split = if(cluster_columns == FALSE){se$Group}else{NULL},
                            row_split = row_split,
                            row_title = NULL,
                            column_title = column_title,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            show_row_names = show_row_names,
                            show_column_names = show_column_names,
                            row_names_gp = grid::gpar(fontsize = 10, fontface = "italic"),
                            show_row_dend = FALSE,
                            border = TRUE,
                            heatmap_legend_param = list(border = TRUE,
                                                        labels = if(is.null(legendLimits)){
                                                          c(floor(min(heatMatrix)),
                                                            0,
                                                            ceiling(max(heatMatrix)))
                                                        }else{
                                                          c(legendLimits[1], 0, legendLimits[2])},
                                                        at = if(is.null(legendLimits)){
                                                          c(floor(min(heatMatrix)),
                                                               0,
                                                               ceiling(max(heatMatrix)))
                                                          }else{
                                                            c(legendLimits[1], 0, legendLimits[2])},
                                                        direction = "horizontal",
                                                        title_position = "topcenter",
                                                        legend_width = if(pathways == FALSE){
                                                          grid::unit(5, "cm")
                                                        }else{
                                                          grid::unit(4, "cm")}),
                            raster_quality = 5,
                            ...)

  if(draw == TRUE){
    heatmap %>%
      ComplexHeatmap::draw(heatmap_legend_side = heatmap_legend_side,
                           annotation_legend_side = annotation_legend_side,
                           ht_gap = grid::unit(0.5, "cm"))
  }else if(draw == FALSE){
    heatmap
  }

}

#' @title scaleRows
#' @description Scale a matrix by rows to generate Z-scores for a heatmap
#' @param mat A Matrix
#' @param verbose Logical indicating whether to print range of Z-scores
#' @export
#'
#' @return A matrix of row z-scores (the rows of `mat` centred and scaled).
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' SummarizedExperiment::assay(se, "logCPM")[getHVGs(se, n = 20, symbols = FALSE), ] %>%
#'     scaleRows()
scaleRows <- function(mat,
                      verbose = TRUE){
  res <- t(scale(t(mat)))

  if(verbose == TRUE){
    res %>%
      range(finite = TRUE) %>%
      round(digits = 2) %>%
      print()
  }

  return(res)
}

#' @title colorRamp3
#' @description Generate a function for a color palette for a heatmap
#' @param mat A Matrix
#' @param heatColors Character vector of colors that has an odd number of colors with a midpoint; the default is
#'  MetBrewer's Troy palette (10 colors, blue to red) around a white midpoint
#' @param legendLimits Manually override the legend limits
#' @return A function
#' @importFrom circlize colorRamp2
#' @importFrom MetBrewer met.brewer
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' SummarizedExperiment::assay(se, "logCPM")[getHVGs(se, n = 20, symbols = FALSE), ] %>%
#'     scaleRows() %>%
#'     colorRamp3()
colorRamp3 <- function(mat,
                       heatColors = MetBrewer::met.brewer("Troy", n = 10, type = "continuous") %>%
                         as.character() %>%
                         rev() %>%
                         append("#f6f2ee", after = 5),
                       legendLimits = NULL){

  stopifnot(length(heatColors) %% 2 == 1)
  palLength <- (length(heatColors)/2) + 0.5
  mat <- mat[is.finite(mat)]

  if(is.null(legendLimits)){
    c(seq(min(mat), 0, length.out = palLength),
      seq(0, max(mat), length.out = palLength)[-1]) %>%
      circlize::colorRamp2(colors = heatColors)
  }else if(length(legendLimits) == 2){
    c(seq(legendLimits[1], 0, length.out = palLength),
      seq(0, legendLimits[2], length.out = palLength)[-1]) %>%
      circlize::colorRamp2(colors = heatColors)
  }else{
    stop("legendLimits not valid")
  }
}

#' @title getHVGs
#' @description Get highly variable genes.
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param assay Character specifying the name of the assay from the se
#' @param n Integer specifying the number of genes.
#' @param symbols Logical indicating to return gene symbols instead of Ensembl IDs.
#' @return A character vector
#' @importFrom SummarizedExperiment assays assay rowData
#' @importFrom rlang is_installed
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' getHVGs(se, n = 20)
getHVGs <- function(se,
                    assay = "logCPM",
                    n = 100,
                    symbols = TRUE){

  stopifnot(rlang::is_installed("MatrixGenerics"))
  stopifnot(assay %in% names(SummarizedExperiment::assays(se)))

  HVGs <- se %>%
    SummarizedExperiment::assay(assay) %>%
    MatrixGenerics::rowVars() %>%
    sort(decreasing = TRUE) %>%
    utils::head(n) %>%
    names()

  if(symbols == TRUE){
    HVGs <- SummarizedExperiment::rowData(se)[HVGs, "symbol"]
  }

  return(HVGs)
}

#' @title extractHeatmap
#' @description Extract the clustered matrix and relevant metadata from a heatmap for figure source data
#' @param se A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object
#' @param heatmap A drawn heatmap from \link{basicHeatmap}
#' @return A tibble
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom dplyr left_join select relocate replace_values
#' @importFrom stringr str_remove
#' @importFrom ComplexHeatmap row_order
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' data(airwayDGE)
#' heatmap <- basicHeatmap(se, sigTibble = makeSigTibble(airwayDGE))
#' extractHeatmap(se, heatmap)
extractHeatmap <- function(se,
                           heatmap){

  heatmap@ht_list[["Z-score"]]@matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::left_join(se %>%
                       SummarizedExperiment::rowData() %>%
                       tibble::as_tibble() %>%
                       dplyr::select(symbol, ID),
                     by = "ID") %>%
    dplyr::relocate(symbol, .after = "ID") %>%
    `colnames<-`(x = ., value = paste(colnames(.) %>%
                                        dplyr::replace_values(from = SummarizedExperiment::colData(se)$Sample,
                                                              to = as.character(SummarizedExperiment::colData(se)$Group)),
                                      colnames(.)) %>%
                   stringr::str_remove(" ID| symbol")) %>%
    .[unlist(ComplexHeatmap::row_order(heatmap), use.names = FALSE), ] %>%
    tibble::as_tibble()
}
