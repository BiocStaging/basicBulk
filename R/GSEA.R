#' @title getReactome
#' @description Get the Reactome pathway collection (MSigDB C2:CP:REACTOME) for GSEA as a named list
#'  of Ensembl gene IDs, via \link[msigdbr]{msigdbr}. Set names follow the \code{REACTOME_*}
#'  convention expected by \link{tidyGSEAlist} and \link{makePathwayTibble}, and the Ensembl IDs
#'  match the row names used by \link[augere.gsea]{runContrast}.
#' @param organism Organism of interest c("Homo sapiens", "Mus musculus")
#' @return A named list of character vectors of Ensembl gene IDs
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr filter distinct group_by summarise
#' @importFrom tibble deframe
#' @importFrom magrittr %>%
#' @references
#' Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES,
#'  Mesirov JP (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression
#'  profiles. \emph{Proceedings of the National Academy of Sciences}, 102(43), 15545-15550.
#'  \doi{10.1073/pnas.0506580102}
#'
#' Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P (2015). The Molecular Signatures Database
#'  hallmark gene set collection. \emph{Cell Systems}, 1(6), 417-425. \doi{10.1016/j.cels.2015.12.004}
#'
#' Ragueneau E, Gong C, Sinquin P, Sevilla C, Beavers D, Grentner A, Griss J, Hogue GFJ, Li NT, Matthews L, May B,
#'  Milacic M, Mohammadi H, Petryszak R, Rothfels K, Shamovsky V, Stephan R, Tiwari K, Weiser J, Wright A, Gillespie M,
#'  Wu G, Stein L, Hermjakob H, D'Eustachio P (2025). The Reactome Knowledgebase 2026. \emph{Nucleic Acids Research},
#'  54(D1), D673-D681. \doi{10.1093/nar/gkaf1223}
#' @export
#'
#' @examples
#' gene.sets <- getReactome("Homo sapiens")
getReactome <- function(organism){

  stopifnot(organism %in% c("Homo sapiens", "Mus musculus"))

  msigdbr::msigdbr(db_species = "HS",
                   species = organism,
                   collection = "C2",
                   subcollection = "CP:REACTOME") %>%
    dplyr::filter(!is.na(ensembl_gene) & ensembl_gene != "") %>%
    dplyr::distinct(gs_name, ensembl_gene) %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(genes = list(ensembl_gene), .groups = "drop") %>%
    tibble::deframe()
}

#' @title tidyGSEAlist
#' @description Create a named list of GSEA tibbles for downstream functions
#' @param GSEAlist A named list of \link[augere.gsea]{runContrast} results, one per contrast
#' @param GSEAmethod Character specifying name of GSEA method to extract
#' @return A list of tibbles
#' @importFrom purrr map
#' @importFrom magrittr %>% extract2
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr separate_wider_delim
#' @importFrom dplyr mutate arrange
#' @importFrom stringr str_replace_all regex
#' @references
#' Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation. \emph{Nucleic
#'  Acids Research}, 40(17), e133. \doi{10.1093/nar/gks461}
#' @export
#'
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
#' se <- tidySE(se0, group = "cell", block = "dex")
#' gene.sets <- getReactome("Homo sapiens")
#' list(N61311 = se %>%
#'          augere.gsea::runContrast(sets = gene.sets,
#'                                   groups = "cell",
#'                                   block = "dex",
#'                                   comparison = c("N61311", "N052611"),
#'                                   methods = "camera",
#'                                   assay = "counts",
#'                                   output.dir = tempfile(),
#'                                   save.results = FALSE)) %>%
#'     tidyGSEAlist()
tidyGSEAlist <- function(GSEAlist,
                         GSEAmethod = "camera"){

  GSEAlist %>%
    purrr::map(~ magrittr::extract2(., GSEAmethod) %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column("name") %>%
                 tibble::as_tibble() %>%
                 tidyr::separate_wider_delim(name,
                                             delim = stringr::regex("[^[:alnum:]]+"),
                                             names = c("Database", "Term"),
                                             too_many = "merge",
                                             too_few = "align_start") %>%
                 dplyr::mutate(Term = Term %>%
                                 stringr::str_replace_all("_", " ")) %>%
                 dplyr::arrange(FDR))
}

#' @title GSEAplot
#' @description Plot GSEA results
#' @param GSEA A tibble of GSEA results for one contrast, an element of the \link{tidyGSEAlist} output
#' @param n_terms Numeric specifying number of terms per a facet (Up and Down)
#' @return A ggplot2
#' @import ggplot2
#' @importFrom tidyr drop_na
#' @importFrom dplyr mutate group_by slice ungroup
#' @importFrom stringr str_trunc
#' @importFrom forcats as_factor fct_rev
#' @importFrom magrittr %>%
#' @importFrom MetBrewer met.brewer
#' @export
#'
#' @examples
#' data(airwayGSEA)
#' GSEAplot(airwayGSEA[[1]])
GSEAplot <- function(GSEA,
                     n_terms = 5){

  GSEA %>%
    tidyr::drop_na(Direction) %>%
    dplyr::mutate("-log10.FDR" = -log10(FDR),
                  Term = Term %>%
                    stringr::str_replace("^[[:lower:]]", toupper) %>%
                    forcats::as_factor() %>%
                    forcats::fct_rev(),
                  Direction = Direction %>%
                    factor(levels = c("Up", "Down"))) %>%
    dplyr::group_by(Direction) %>%
    dplyr::slice(1:n_terms) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(ggplot2::aes(x = Term,
                                 y = `-log10.FDR`,
                                 fill = `Direction`,
                                 group = `Database`)) +
    ggplot2::scale_fill_manual(values = MetBrewer::met.brewer("VanGogh2")[c(1, 8)] %>%
                                 as.character() %>%
                                 purrr::set_names(c("Up", "Down"))) +
    basicBulk::barStyle() +
    ggplot2::scale_x_discrete(labels = ~ stringr::str_trunc(.x, 45)) +
    ggplot2::labs(y = expression("-Log"["10"](italic(FDR)))) +
    ggplot2::geom_hline(color = ggplot2::alpha("black", 0.5),
                        linewidth = 0.75,
                        yintercept = -log10(0.05)) +
    ggplot2::facet_grid(Direction ~ ., scales = "free") +
    ggplot2::theme(legend.position = "none")
}

#' @title GSEAplotList
#' @description Create a list of GSEA plots for all comparisons
#' @param GSEAlist A named list of GSEA result tibbles from \link{tidyGSEAlist}
#' @inheritParams GSEAplot
#' @return A list of ggplots
#' @import ggplot2
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayGSEA)
#' GSEAplotList(airwayGSEA)
GSEAplotList <- function(GSEAlist,
                         n_terms = 5){

  GSEAlist %>%
    purrr::map(~ .x %>%
                 basicBulk::GSEAplot(n_terms = n_terms))
}

#' @title printGSEAtibble
#' @description Print a table of significant GSEA results
#' @param GSEAlist A named list of GSEA result tibbles from \link{tidyGSEAlist}
#' @param BiocStyle Logical indicating whether the table is rendered in a BiocStyle vignette, where the per-column filter boxes misalign, so they are dropped in favour of the global search
#' @return A datatable
#' @import DT
#' @importFrom purrr map list_rbind
#' @importFrom dplyr select arrange mutate
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayGSEA)
#' printGSEAtibble(airwayGSEA)
printGSEAtibble <- function(GSEAlist,
                            BiocStyle = FALSE){

  GSEAlist %>%
    purrr::map(~ .x %>%
                 dplyr::filter(FDR < 0.05)) %>%
    purrr::list_rbind(names_to = "Contrast") %>%
    dplyr::select(-Database) %>%
    dplyr::mutate(Direction = Direction %>%
                    as.factor() %>%
                    forcats::fct_rev()) %>%
    dplyr::arrange(Contrast, FDR) %>%
    DT::datatable(filter = ifelse(BiocStyle == FALSE, "top", "none"),
                  extensions = "Buttons",
                  rownames = FALSE,
                  options = list(dom = 'Bfrtip',
                                 buttons = c('csv', 'excel'))) %>%
    DT::formatSignif(columns = c("PValue", "FDR"), digits = 2)
}

#' @title topGSEAplot
#' @description A plot of the top GSEA terms from each contrast
#' @param GSEAlist A list of results from GSEA
#' @param n_terms Numeric indicating number of terms per a direction from each contrast
#' @return A ggplot2
#' @importFrom purrr map list_rbind
#' @importFrom dplyr filter mutate group_by slice ungroup pull arrange case_when
#' @importFrom stats na.omit
#' @importFrom scales label_wrap
#' @importFrom forcats as_factor
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom MetBrewer met.brewer
#' @export
#'
#' @examples
#' data(airwayGSEA)
#' topGSEAplot(airwayGSEA)
topGSEAplot <- function(GSEAlist,
                        n_terms = 2){

  terms <- GSEAlist %>%
    purrr::map(~ .x %>%
                 dplyr::mutate(Direction = factor(Direction, levels = c("Up", "Down"))) %>%
                 dplyr::group_by(Direction) %>%
                 dplyr::slice(1:n_terms) %>%
                 dplyr::ungroup()) %>%
    purrr::list_rbind() %>%
    dplyr::filter(FDR < 0.05) %>%
    dplyr::select(Term) %>%
    dplyr::distinct()

  GSEAlevels <- terms %>%
    dplyr::pull(Term) %>%
    unique() %>%
    forcats::as_factor()

  GSEA <- GSEAlist %>%
    purrr::list_rbind(names_to = "contrast") %>%
    dplyr::mutate(contrast = contrast %>% forcats::as_factor(),
                  "-log10.FDR" = -log10(FDR))

  plotData <- GSEA %>%
    dplyr::right_join(terms, by = "Term") %>%
    stats::na.omit() %>%
    dplyr::arrange(`-log10.FDR`, .by_group = TRUE) %>%
    dplyr::mutate(Term = Term %>%
                    factor(levels = GSEAlevels) %>%
                    droplevels(),
                  alpha = dplyr::case_when(FDR >= 0.05 ~ 0.1,
                                           FDR < 0.05 ~ 1))

  plotData %>%
    ggplot2::ggplot(ggplot2::aes(x = contrast,
                                 y = Term,
                                 size = `-log10.FDR`,
                                 fill = Direction,
                                 alpha = I(alpha))) +
    ggplot2::geom_point(color = "black",
                        pch = 21) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right",
                   text = ggplot2::element_text(size = 10),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1,
                                                       size = 10),
                   axis.text.y = ggplot2::element_text(size = 10)) +
    ggplot2::scale_fill_manual(breaks = c("Up", "Down"),
                               values = MetBrewer::met.brewer("VanGogh2")[c(1, 8)] %>%
                                 as.character() %>%
                                 purrr::set_names(c("Up", "Down"))) +
    ggplot2::labs(size = expression("-log"[10](italic(FDR)))) +
    ggplot2::scale_y_discrete(labels = scales::label_wrap(45), limits = rev) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5)),
                    size = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_alpha(guide = "none")
}

#' @title makePathwayTibble
#' @description Make a tibble of pathways for \link{basicHeatmap}
#' @param pathways A character vector of pathway terms from \link{printGSEAtibble}
#' @param gene.sets Gene sets from \link{getReactome}
#' @return A tibble
#' @importFrom AnnotationDbi unlist2
#' @importFrom tibble enframe
#' @importFrom dplyr mutate
#' @importFrom stringr str_remove str_replace_all
#' @importFrom forcats as_factor
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' data(airwayGSEA)
#' gene.sets <- getReactome("Homo sapiens")
#' airwayGSEA[[1]] %>%
#'     dplyr::slice(1:2) %>%
#'     dplyr::pull(Term) %>%
#'     stringr::str_replace_all(" ", "_") %>%
#'     paste0("REACTOME_", .) %>%
#'     makePathwayTibble(gene.sets)
makePathwayTibble <- function(pathways,
                              gene.sets){

 gene.sets[names(gene.sets) %in% pathways] %>%
    as.list() %>%
    AnnotationDbi::unlist2() %>%
    tibble::enframe(name = "Pathway", value = "ID") %>%
    dplyr::mutate(Pathway = Pathway %>%
                    stringr::str_remove("REACTOME_") %>%
                    stringr::str_replace_all("_", " ") %>%
                    forcats::as_factor())
}
