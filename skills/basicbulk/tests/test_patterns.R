# Checks that the design-pattern recipes in references/design-patterns.md produce the expected contrasts
# (count, names, and direction) on mock SummarizedExperiments. Needs basicBulk and augere.de installed.
# Run with: Rscript skills/basicbulk/tests/test_patterns.R
library(basicBulk)
library(SummarizedExperiment)
library(magrittr)

# Helper: create a mock SE with random counts, a counts assay and rowData ID/symbol
make_mock_se <- function(colData, n_genes = 500) {
  set.seed(42)
  n_samples <- nrow(colData)
  counts <- matrix(rpois(n_genes * n_samples, lambda = 100), nrow = n_genes)
  rownames(counts) <- paste0("ENSG000000", sprintf("%05d", seq_len(n_genes)))
  colnames(counts) <- paste0("S", sprintf("%02d", seq_len(n_samples)))
  cd <- S4Vectors::DataFrame(colData)
  rownames(cd) <- colnames(counts)
  se <- SummarizedExperiment(assays = list(counts = counts), colData = cd)
  rowData(se)$ID <- rownames(se)
  rowData(se)$symbol <- paste0("GENE", seq_len(n_genes))
  se
}

one_factor <- function(groups, n_reps = 3) {
  tibble::tibble(Group = factor(rep(groups, each = n_reps), levels = groups))
}

# =============================================================================
# Pattern 1: Single model, simple comparisons
# =============================================================================
cat("\n=== Pattern 1 ===\n")

se1 <- make_mock_se(one_factor(c("PBS", "DrugA", "DrugB", "DrugC")))

comparisonList <- levels(se1$Group)[c(2:nlevels(se1$Group))] %>%
  purrr::map(~ list(.x, levels(se1$Group)[1]) %>%
               purrr::list_c())

voomList <- se1 %>%
  augere.de::runVoom(groups = "Group",
                     comparisons = comparisonList,
                     assay = "counts",
                     output.dir = tempfile(),
                     save.results = FALSE) %>%
  magrittr::extract2("results")

cat("Pattern 1 contrasts:", length(voomList), "\n")
cat("Names:", paste(names(voomList), collapse = ", "), "\n")
stopifnot(length(voomList) == 3,
          identical(names(se1 %>% basicBulk::makeDGEdata(voomList)),
                    c("DrugA vs PBS", "DrugB vs PBS", "DrugC vs PBS")))

# =============================================================================
# Pattern 2: Compound group names, explicit pairings via designTibble
# =============================================================================
cat("\n=== Pattern 2 ===\n")

se2 <- make_mock_se(one_factor(c("NeuronPBS", "NeuronDrug", "AstroPBS", "AstroDrug")))

designTibble <- tibble::tibble(
  exp = c("NeuronDrug", "AstroDrug"),
  control = c("NeuronPBS", "AstroPBS"),
  name = glue::glue("{exp} vs {control}")
)

voomList2 <- designTibble %>%
  purrr::pmap(~ with(list(...),
                     augere.de::runVoom(se2,
                                        groups = "Group",
                                        comparisons = c(exp, control),
                                        assay = "counts",
                                        output.dir = tempfile(),
                                        save.results = FALSE) %>%
                       magrittr::extract2("results")),
              .progress = TRUE) %>%
  purrr::list_c()

cat("Pattern 2 contrasts:", length(voomList2), "\n")
cat("Names:", paste(names(voomList2), collapse = ", "), "\n")
stopifnot(length(voomList2) == 2,
          identical(names(se2 %>% basicBulk::makeDGEdata(voomList2)),
                    as.character(designTibble$name)))

# =============================================================================
# Pattern 3: Single-variable subset via subset.factor (guarded recipe)
# =============================================================================
cat("\n=== Pattern 3 ===\n")

pattern3 <- function(se) {
  control <- levels(se$Group)[1]

  comparisonList <- levels(se$Group)[c(2:nlevels(se$Group))] %>%
    purrr::map(~ list(.x, control) %>%
                 purrr::list_c())

  unique(se$CellType) %>%
    forcats::as_factor() %>%
    purrr::set_names() %>%
    purrr::map(function(subset_level) {
      subset_groups <- se %>%
        dplyr::filter(CellType == subset_level) %>%
        dplyr::pull(Group) %>%
        forcats::fct_drop() %>%
        levels()

      if (!control %in% subset_groups || length(subset_groups) < 2) {
        message(glue::glue("Skipping {subset_level}: control missing or fewer than 2 groups"))
        return(NULL)
      }

      valid_comparisons <- comparisonList %>%
        purrr::keep(~ all(.x %in% subset_groups))

      if (length(valid_comparisons) == 0) {
        message(glue::glue("Skipping {subset_level}: no valid comparisons"))
        return(NULL)
      }

      results <- augere.de::runVoom(se,
                                    groups = "Group",
                                    comparisons = valid_comparisons,
                                    subset.factor = "CellType",
                                    subset.levels = subset_level,
                                    assay = "counts",
                                    output.dir = tempfile(),
                                    save.results = FALSE) %>%
        magrittr::extract2("results")
      purrr::set_names(results, glue::glue("{subset_level}:{names(results)}"))
    }, .progress = TRUE) %>%
    purrr::compact() %>%
    purrr::list_c()
}

design3 <- tidyr::crossing(
  Group = factor(c("PBS", "DrugA", "DrugB"), levels = c("PBS", "DrugA", "DrugB")),
  CellType = factor(c("Neuron", "Astro", "Microglia"), levels = c("Neuron", "Astro", "Microglia"))
) %>%
  dplyr::slice(rep(seq_len(dplyr::n()), each = 3))

voomList3 <- pattern3(make_mock_se(design3))
cat("Pattern 3 contrasts:", length(voomList3), "\n")
cat("Names:", paste(names(voomList3), collapse = ", "), "\n")
stopifnot(length(voomList3) == 6,  # 3 cell types x 2 comparisons
          all(stringr::str_detect(names(voomList3), "over `PBS`$")))

# Unbalanced: Microglia has no PBS samples, so it is skipped rather than re-baselined
voomList3u <- pattern3(make_mock_se(dplyr::filter(design3, !(CellType == "Microglia" & Group == "PBS"))))
cat("Pattern 3 (unbalanced) contrasts:", length(voomList3u), "\n")
stopifnot(length(voomList3u) == 4,
          all(stringr::str_detect(names(voomList3u), "over `PBS`$")))

# =============================================================================
# Pattern 4: Crossed multi-factor split (guarded recipe)
# =============================================================================
cat("\n=== Pattern 4 ===\n")

pattern4 <- function(se) {
  control <- levels(se$Group)[1]

  comparisonList <- levels(se$Group)[c(2:nlevels(se$Group))] %>%
    purrr::map(~ list(.x, control) %>%
                 purrr::list_c())

  splitVoom <- function(Region, Cohort) {
    splitSe <- se %>%
      dplyr::filter(Region == !!Region & Cohort == !!Cohort) %>%
      dplyr::mutate(Group = forcats::fct_drop(Group))

    split_groups <- levels(splitSe$Group)

    if (!control %in% split_groups || length(split_groups) < 2) {
      message(glue::glue("Skipping {Region}:{Cohort}: control missing or fewer than 2 groups"))
      return(NULL)
    }

    valid_comparisons <- comparisonList %>%
      purrr::keep(~ all(.x %in% split_groups))

    if (length(valid_comparisons) == 0) {
      message(glue::glue("Skipping {Region}:{Cohort}: no valid comparisons"))
      return(NULL)
    }

    voomList <- splitSe %>%
      augere.de::runVoom(groups = "Group",
                         comparisons = valid_comparisons,
                         assay = "counts",
                         output.dir = tempfile(),
                         save.results = FALSE) %>%
      magrittr::extract2("results")

    splitSe %>%
      basicBulk::makeDGEdata(voomList = voomList) %>%
      purrr::set_names(glue::glue("{Region}:{Cohort}:{names(.)}"))
  }

  tidyr::crossing(
    Region = unique(se$Region) %>% forcats::as_factor(),
    Cohort = unique(se$Cohort) %>% forcats::as_factor()
  ) %>%
    purrr::pmap(~ with(list(...),
                       splitVoom(Region = Region, Cohort = Cohort)),
                .progress = TRUE) %>%
    purrr::compact() %>%
    purrr::list_c()
}

design4 <- tidyr::crossing(
  Group = factor(c("Vehicle", "DrugA", "DrugB"), levels = c("Vehicle", "DrugA", "DrugB")),
  Region = factor(c("Cortex", "Cerebellum"), levels = c("Cortex", "Cerebellum")),
  Cohort = factor(c("CohortA", "CohortB"), levels = c("CohortA", "CohortB"))
) %>%
  dplyr::slice(rep(seq_len(dplyr::n()), each = 3))

DGEdata4 <- pattern4(make_mock_se(design4))
cat("Pattern 4 contrasts:", length(DGEdata4), "\n")
cat("Names:", paste(names(DGEdata4), collapse = ", "), "\n")
stopifnot(length(DGEdata4) == 8,  # 2 regions x 2 cohorts x 2 comparisons
          all(stringr::str_detect(names(DGEdata4), " vs Vehicle$")))

# Unbalanced: Cerebellum:CohortB has no Vehicle samples, so it is skipped rather than re-baselined
DGEdata4u <- pattern4(make_mock_se(dplyr::filter(design4, !(Region == "Cerebellum" & Cohort == "CohortB" & Group == "Vehicle"))))
cat("Pattern 4 (unbalanced) contrasts:", length(DGEdata4u), "\n")
stopifnot(length(DGEdata4u) == 6,
          all(stringr::str_detect(names(DGEdata4u), " vs Vehicle$")))

# =============================================================================
cat("\n=== ALL PATTERNS PASS ===\n")
