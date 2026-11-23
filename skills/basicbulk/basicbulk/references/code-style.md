# Code Style Reference

Generated code must match the author's tidyverse + magrittr style. These are not preferences; they are requirements for the skill's output.

The core style rules stay inline in `SKILL.md`. Use this reference for examples, edge cases, and less common idioms.

## Pipe: Magrittr Only

- Always `%>%`, never base R `|>`.
- Use `%T>%` (tee pipe) for side effects that should not break the chain, especially `ggplot2::ggsave()`:

```r
(plot_expression) %T>%
  ggplot2::ggsave("plots/DGE/volcano.png", plot = ., width = 6, height = 4)
```

- Use `%$%` (exposition pipe) to expose columns when needed:

```r
designTibble %$% list(exp, control) %>% purrr::list_transpose()
```

- Use `magrittr::extract2()` instead of `[[` in pipes:

```r
voomList %>% magrittr::extract2("results")
DGEdata %>% magrittr::extract2(2)
```

- Use `magrittr::extract()` instead of `[` in pipes.

## Namespace Style

- Fully qualify package functions: `dplyr::filter()`, `dplyr::mutate()`, `ggplot2::ggplot()`, `SummarizedExperiment::colData()`, `basicBulk::basicVolcano()`, etc.
- The only exceptions are:
  - `factor()` is acceptable in simple setup code when explicitly setting analysis group levels. In tidyverse pipelines and other factor manipulation, prefer `forcats::*`.
  - Base helpers such as `grepl()` are acceptable in shell/Rscript inspection snippets. Generated analysis code should use `stringr::*` for string operations.
- Data columns inside `ggplot2::aes()` are not function calls and do not need namespaces.
- Use `purrr::set_names()`, not `stats::setNames()` or `rlang::set_names()`.

## `purrr` Idioms

- Use `purrr::pmap(~ with(list(...), ...))` for iterating over tibble rows with named access:

```r
designTibble %>%
  purrr::pmap(~ with(list(...),
                     augere.de::runVoom(se, groups = "Group", comparisons = c(exp, control),
                                        assay = "counts", ...)))
```

- Use `purrr::imap(~ ...)` for named iteration where `.y` is the name.
- Use `purrr::iwalk(~ ...)` for side-effect iteration, including printing and saving.
- Use `purrr::set_names()` to name list elements inline:

```r
comparisonList %>%
  purrr::set_names(purrr::map_chr(., magrittr::extract, 1)) %>%
  purrr::imap(~ ...)
```

- Use `purrr::pluck()` instead of `$` for named extraction inside a pipe (for `[[`, use `magrittr::extract2()` as above):

```r
se %>%
  SummarizedExperiment::colData() %>%
  purrr::pluck("Group")
```

- Use `.progress = TRUE` on long-running `purrr::map()` and `purrr::pmap()` calls.

## Tidy Evaluation

- Use `!!` (bang-bang) for injecting variables into dplyr expressions:

```r
dplyr::filter(Group == !!control)
dplyr::select(exp = !!group)
```

- Use `!!!` (bang-bang-bang) for splicing named vectors:

```r
dplyr::recode_factor(exp, !!!groups)
```

- Use `rlang::sym()` for programmatic column access:

```r
ggplot2::aes(x = !!rlang::sym(groupVar))
```

- Use `.data[[var]]` in newer code, especially ggplot2 contexts:

```r
ggplot2::aes(fill = .data[[groupVar]])
```

## Pipe Conditionals

Write a conditional step inside a pipe chain as `{if(cond){ ... }else{ . }}`:

```r
plot %>%
  {if(!is.null(pal)){
    . + ggplot2::scale_color_manual(values = pal)
  }else{
    .
  }}
```

## Factor Handling

- Use `factor()` only in simple setup code when explicitly setting analysis group levels, for example `factor(levels = c("Control", "Treatment"))`.
- Use `forcats::as_factor()` to preserve order in tidyverse pipelines, not `as.factor()` or `factor()`.
- Use `forcats::fct_rev()`, `forcats::fct_drop()`, and `forcats::fct_relevel()` for manipulation.
- Relevel the control group so control is `levels(...)[1]`.

## String Operations

- Use `stringr::` functions in generated analysis code, never base `gsub()`/`sub()`/`grep()`/`grepl()`.
- Base helpers such as `grepl()` are acceptable in shell/Rscript inspection snippets.
- Use `glue::glue()` for string interpolation, especially file paths and dynamic names:

```r
glue::glue("plots/DGE/{.y}_volcano.png")
glue::glue("{exp} vs {control}")
```

- Use `glue::glue_collapse()` for joining vectors into strings.

## Tibble Construction

- Use `tibble::tribble()` for inline data with named columns.
- Use `tidyr::crossing()` for parameter grids:

```r
tidyr::crossing(Region = levels(se$Region),
                Study = levels(se$Study))
```

## ggplot2 Patterns

- Build plots in a `(expression) %T>% ggplot2::ggsave(...)` pattern for simultaneous display and save.
- Wrap ggplot chains that end with `+` in `(...)` before `%T>%`. Without parens, `%T>%` binds to the last `+` operand (e.g., `facet_nested()`), not the full plot, causing `ggsave()` to fail. This applies to any chain ending in `+ ggh4x::facet_nested(...)`, `+ ggplot2::theme(...)`, `+ ggplot2::labs(...)`, etc. (`ggh4x` is a Suggests-level package of basicBulk, installed with `dependencies = TRUE`.)
- Use `ggsci::scale_fill_d3("category20")` and `ggsci::scale_color_d3()` as default palettes.
- Use `ggpubr::rotate_x_text()` for angled x-axis labels.
- Use `ggplot2::theme_bw()` as the base theme.

## Data Access Patterns

- Define `plotAssay` once as `"logCPM"` (or `"ComBat"` with `adjustBatch = TRUE`) after `basicBulk::tidySE()`, then pass `assay = plotAssay` to the basicBulk functions that read expression values (`PCAgrid()`, `basicHeatmap()`, `genePlot()`, `geneHistogram()`, `getHVGs()`, `joinFeatures()`). These default to `assay = "logCPM"`, so leaving out `assay = plotAssay` silently drops the ComBat adjustment when `adjustBatch = TRUE`. Statistics use the same `se` with `assay = "counts"` in the augere calls (`augere.de::runVoom()`, `augere.gsea::runContrast()`).
- Use `SummarizedExperiment::colData()`, `SummarizedExperiment::rowData()`, and `SummarizedExperiment::assay()` fully qualified.
- Use `tibble::as_tibble()` to convert S4 objects to tibbles for dplyr operations.

## Custom Operators

`%ni%` means "not in" and is defined as:

```r
`%ni%` <- purrr::negate(`%in%`)
```
