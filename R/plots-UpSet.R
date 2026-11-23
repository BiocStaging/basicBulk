#' @title basicUpset
#' @description A wrapper for \link[ComplexUpset]{upset}
#' @param DGEdata A list of DGE results
#' @param directionSplit Logical indicating whether gene lists should be split by direction, default is TRUE
#' @param FDRcutoff Numeric for the FDR cut-off for DEGs, default is 0.05
#' @param logFCcutoff Numeric for the absolute Log2FC cut-off for DEGs, default is 1
#' @param fontSize Integer specifying font size
#' @param fontSizeOverall Integer specifying font size for overall set size axis labels
#' @inheritParams ComplexUpset::upset
#' @inheritParams ComplexUpset::upset_data
#' @inheritParams ComplexUpset::intersection_size
#' @return A ggplot2
#' @import ComplexUpset
#' @importFrom purrr set_names imap list_rbind map pluck
#' @importFrom stringr str_replace_all
#' @importFrom tibble lst
#' @importFrom dplyr filter mutate case_when
#' @importFrom forcats fct_rev
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom MetBrewer met.brewer
#' @references
#' Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H (2014). UpSet: visualization of intersecting sets.
#'  \emph{IEEE Transactions on Visualization and Computer Graphics}, 20(12), 1983-1992. \doi{10.1109/TVCG.2014.2346248}
#'
#' Krassowski M (2020). ComplexUpset. \doi{10.5281/zenodo.3700590}
#' @export
#'
#' @examples
#' data(airwayDGE)
#' basicUpset(airwayDGE)
basicUpset <- function(DGEdata,
                       directionSplit = TRUE,
                       FDRcutoff = 0.05,
                       logFCcutoff = 1,
                       n_intersections = 40,
                       width_ratio = 0.1,
                       height_ratio = 0.95,
                       counts = TRUE,
                       sort_sets = "descending",
                       fontSize = 12,
                       fontSizeOverall = 9){

  DGEdata <- DGEdata %>%
    purrr::set_names(names(.) %>%
                       stringr::str_replace_all("_", "-"))

  if(directionSplit == TRUE){

    geneList <- DGEdata %>%
      purrr::imap(~ .x %>%
                    dplyr::filter(FDR < FDRcutoff & abs(LogFC) > logFCcutoff) %>%
                    dplyr::mutate(Direction = dplyr::case_when(LogFC > 0 ~ glue::glue("{.y} Up"),
                                                               LogFC < 0 ~ glue::glue("{.y} Down")) %>%
                                    as.factor() %>%
                                    forcats::fct_rev())) %>%
      purrr::list_rbind() %>%
      split(.$Direction) %>%
      base::Filter(function(x) dim(x)[1] > 0, .) %>%
      purrr::map(purrr::pluck("ID"))

    topSets <- geneList %>%
      rev() %>%
      basicBulk::fromList2() %>%
      ComplexUpset::upset_data(.,
                               names(.),
                               n_intersections = n_intersections,
                               sort_sets = sort_sets) %>%
      magrittr::extract2("sets_ordering_in_ids") %>%
      stringr::str_replace_all("_", "-")

    queryList <- geneList %>%
      names() %>%
      tibble::tibble(set = .) %>%
      dplyr::filter(set %in% topSets) %>%
      dplyr::mutate(fill = dplyr::case_when(stringr::str_detect(set, "Up") ~ MetBrewer::met.brewer("Cassatt1")[1] %>% as.character(),
                                            stringr::str_detect(set, "Down") ~ MetBrewer::met.brewer("Cassatt1")[8] %>% as.character())) %>%
      purrr::pmap(~ with(list(...),
                         tibble::lst(set = {{set}},
                                     intersect = NULL,
                                     group_by_group = NULL,
                                     only_components = NULL,
                                     fill = {{fill}})))

    matrix_aes <- ComplexUpset::intersection_matrix(geom = ggplot2::geom_point(shape = "circle filled",
                                                                               size = 3.5,
                                                                               stroke = 0.45))

  }else if(directionSplit == FALSE){

    geneList <- DGEdata %>%
      purrr::map(~ .x %>%
                   dplyr::filter(FDR < FDRcutoff & abs(LogFC) > logFCcutoff) %>%
                   purrr::pluck("ID"))

    queryList <- NULL
    matrix_aes <- ComplexUpset::intersection_matrix()
  }

  geneList %>%
    rev() %>%
    basicBulk::fromList2() %>%
    ComplexUpset::upset(.,
                        names(.),
                        n_intersections = n_intersections,
                        width_ratio = width_ratio,
                        height_ratio = height_ratio,
                        base_annotations = list("Intersection size"= ComplexUpset::intersection_size(counts = counts,
                                                                                                     text = list(angle = 90,
                                                                                                                 vjust = 0.5,
                                                                                                                 hjust = -0.001)) +
                                                  ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                                                                 axis.text = ggplot2::element_text(size = fontSize),
                                                                 axis.title = ggplot2::element_text(size = fontSize)) +
                                                  ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0))),
                        sort_sets = sort_sets,
                        set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(color = "black"),
                                                                 position = "right"),
                        queries = queryList,
                        themes = ComplexUpset::upset_modify_themes(
                          list("intersections_matrix" = ggplot2::theme(axis.text = ggplot2::element_text(size = fontSize),
                                                                       axis.title = ggplot2::element_blank()),
                               "overall_sizes" = ggplot2::theme(axis.text = ggplot2::element_text(size = fontSizeOverall, angle = 45, hjust = 1),
                                                                axis.title = ggplot2::element_text(size = fontSize)))),
                        matrix = matrix_aes)
}

#' @title fromList2
#' @description Convert a named list of sets to a data.frame of set membership with the elements as rownames
#' @param input A named list of character vectors
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @export
#'
#' @return A data.frame of 0/1 set membership, with one row per element (rownames) and one column per set.
#' @examples
#' data(airwayDGE)
#' airwayDGE %>%
#'     purrr::map(~ dplyr::filter(.x, FDR < 0.05) %>% dplyr::pull(ID)) %>%
#'     fromList2()
fromList2 <- function(input){

  elements <- input %>%
    unlist() %>%
    unique()

  input %>%
    purrr::map(~ as.integer(elements %in% .x)) %>%
    data.frame(row.names = elements,
               check.names = FALSE)
}

#' @title extractUpSet
#' @description Extract gene IDs from an intersection in an UpSet plot
#' @param DGEcontrasts A character vector of the UpSet sets to intersect: \code{DGEdata} names followed by
#'   \code{" Up"} or \code{" Down"}
#' @param DGEdata A list of DGE results
#' @param FDRcutoff Numeric for the FDR cut-off for DEGs, default is 0.05
#' @param logFCcutoff Numeric for the absolute Log2FC cut-off for DEGs, default is 1
#' @importFrom purrr imap list_rbind map pluck
#' @importFrom dplyr filter mutate case_when if_all all_of
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @export
#'
#' @return A character vector of the gene IDs in the requested overlap.
#' @examples
#' data(airwayDGE)
#' names(airwayDGE) %>%
#'     paste("Up") %>%
#'     extractUpSet(airwayDGE)
extractUpSet <- function(DGEcontrasts,
                         DGEdata,
                         FDRcutoff = 0.05,
                         logFCcutoff = 1){

  DGEdata %>%
    purrr::imap(~ .x %>%
                  dplyr::filter(FDR < FDRcutoff & abs(LogFC) > logFCcutoff) %>%
                  dplyr::mutate(Direction = dplyr::case_when(LogFC > 0 ~ glue::glue("{.y} Up"),
                                                             LogFC < 0 ~ glue::glue("{.y} Down")))) %>%
    purrr::list_rbind() %>%
    dplyr::mutate(Direction = factor(Direction,
                                     levels = unique(paste(rep(names(DGEdata), each = 2), c("Up", "Down"))))) %>%
    split(.$Direction) %>%
    purrr::map(purrr::pluck("ID")) %>%
    basicBulk::fromList2() %>%
    as.data.frame() %>%
    dplyr::filter(dplyr::if_all(.cols = dplyr::all_of(DGEcontrasts), .fns = ~ .x == 1)) %>%
    dplyr::filter(rowSums(.) == length(DGEcontrasts)) %>%
    rownames()
}
