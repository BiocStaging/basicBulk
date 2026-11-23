# basicBulk

Tidyverse wrappers, tables and plots for (pseudo)bulk RNA-seq analysis of a `SummarizedExperiment`. Differential
expression and gene set testing come from `augere.de::runVoom()` and `augere.gsea::runContrast()` (Imports, declared
with `@importFrom` in `R/basicBulk-package.R` because basicBulk's own code never calls them).

## Commands

From the package root:

- Document (regenerates `NAMESPACE` and `man/`): `Rscript -e 'roxygen2::roxygenise()'`
- Install: `R CMD INSTALL .`
- Render the vignette: `Rscript -e 'rmarkdown::render("vignettes/basicBulk.Rmd")'` (msigdbr's first MSigDB download
  needs network access).
- Rebuild the example data: install this tree first (the script uses the installed basicBulk), then
  `Rscript data-raw/airwayExample.R`, and remove the `LazyData: true` line its `usethis::use_data()` adds to
  DESCRIPTION.
- Check: `R CMD build . && R CMD check basicBulk_<version>.tar.gz`, then `Rscript -e 'BiocCheck::BiocCheck(".")'` (it
  errors unless the package folder is named `basicBulk`).
- Skill recipe test: `Rscript skills/basicbulk/tests/test_patterns.R`.

## Layout

- `R/normalize.R`: `tidySE()`. `R/utils-DGE.R`: `makeDGEdata()`, `makeSigTibble()`, `makeSigCountTibble()`,
  `joinFeatures()`. `R/GSEA.R`: `getReactome()`, `tidyGSEAlist()`, GSEA plots and table, `makePathwayTibble()`.
  `R/plots-*.R`: PCA, gene (`genePlot()` also plots gene set scores), volcano, MA, DEG bar, heatmap, UpSet and
  logFC-correlation plots.
- `R/utils-knit.R`: report helpers (the `print*Tibble()` tables, `printPlots()`, `printHeatmaps()`).
  `R/utils-slides.R`: cache helpers for the skill's reports and slides (`save_cache()`, `load_if_fresh()`,
  `assert_fresh()`).
- `data/`: the example objects `airwayDGE` (`DGEdata`) and `airwayGSEA` (tidied GSEA), documented in `R/data.R` and
  rebuilt by `data-raw/airwayExample.R`.
- `vignettes/basicBulk.Rmd`: the only vignette (BiocStyle; airway data from `augere.de::loadExampleDataset()`), with
  its bibliography in `vignettes/references.bib`.
- `skills/basicbulk/` and `.claude-plugin/marketplace.json`: the Claude Code plugin (marketplace and plugin both named
  `basicbulk`; skill `/basicbulk:basicbulk`). These and this file are `.Rbuildignore`d. DESCRIPTION Suggests `ggh4x`,
  `furrr`, `future`, `parallelly`, `patchwork`, `png`, `org.Hs.eg.db` and `org.Mm.eg.db` only for the skill: do not drop
  them as unused.

## Data contract

Functions take a `SummarizedExperiment` with a `counts` assay and `rowData` columns `ID` (a real column, not just
rownames: many joins use `by = "ID"`) and `symbol`. `tidySE()` adds colData `Group`, `Sample`, `lib.size` and
`norm.factors`, `rowData$retained`, and the assays `logCPM` and (with `adjustBatch = TRUE`) `ComBat`.

## tidySE mirrors augere.de

Keep these matched to augere.de when editing:

- `filterByExpr(group = Group, large.n = 0)` on the same recoded `Group` used for the design, filtering inside the
  DGEList, then TMM factors (`edgeR::normLibSizes()`); that DGEList feeds voom.
- `logCPM` is `voomWithQualityWeights()`'s `E`, and `ComBat` the voom `E` of the ComBat-seq counts; both are
  `NA`-padded to all genes by row identity. The plots use these voom values, not augere.de's `normalized` export.
- `.designMatrix()` builds safe term names with `make.names(unique = TRUE)` but extracts values by the original colData
  name.

## augere API

- `runVoom(x, groups=, comparisons=, block=, assay=, output.dir=, save.results=)` returns `list(results, normalized)`;
  `results` holds one `DataFrame` per contrast, named ``"Increase in `X` over `Y`"``, which `makeDGEdata()` turns into
  `DGEdata` (a named list of tibbles, the input gg4way expects).
- `runContrast(x, sets, groups=, comparison=, block=, methods=, output.dir=, save.results=)`: `sets` is the second
  positional argument, `comparison` is singular, and there is no seed argument (call `set.seed()` first).
- Pass a fresh `tempfile()` as `output.dir` per call. A named comparison vector is read as two named group sets ordered
  by name, which can silently reverse the direction.
- `getReactome()` is MSigDB C2:CP:REACTOME (`db_species = "HS"`); the native mouse collection is `db_species = "MM"`,
  `collection = "M2"`. Do not substitute one for the other silently.

## Examples

Examples use real data: an SE built as the vignette does (`se0 <- augere.de::loadExampleDataset()`,
`SummarizedExperiment::rowData(se0)$ID <- rownames(se0)`, `se <- tidySE(se0, group = "cell", block = "dex")`) and/or
`data(airwayDGE)`/`data(airwayGSEA)`, plus the call, kept as short as the function allows. The msigdbr-based examples
stay runnable, not in `\donttest{}`, which BiocCheck counts as not runnable.

## Code style

Match the existing code exactly: magrittr `%>%`/`%T>%`; fully qualified calls (`dplyr::filter()`, and
`basicBulk::fromList2()` for the package's own exported functions); purrr (row-bind with `purrr::map()`/`imap()` +
`purrr::list_rbind()`); pipe conditionals as `{if(cond){ ... }else{ . }}`; bare column names and `.`, not `.data`;
braced, cuddled `if(){}else{}`; `dplyr::recode_values()`/`dplyr::replace_values()` for value recoding; `== TRUE`;
`%||%`; `forcats::as_factor()`; no unnecessary temporary variables; continuation arguments aligned under the first
argument. Colours are `MetBrewer::met.brewer()` calls, not pasted hex values. A Suggests package used in package code
is guarded with `stopifnot(rlang::is_installed("pkg"))` where it is used. Use the current tidyverse/ggplot2 calls:
`ggplot2::labs()` for titles, `tidyr::separate_wider_delim()`, `scales::label_*()`, and guide styling through the
guide's `theme =` argument.

## Gotchas

- `library(basicBulk)` attaches SummarizedExperiment (Depends) but not the tidyverse, so vignette, example and skill
  code must qualify tidyverse calls (a bare `filter()` runs `stats::filter()`).
- `basicHeatmap()` drops all-`NA`/`NaN` rows after `scaleRows()`. Overlapping pathway sets can desync the annotation
  length: dedupe with `dplyr::distinct(ID, .keep_all = TRUE)`. Its row clustering depends on the input row order for
  large gene sets, so a symbol `geneSet` is put in its own order only with `cluster_rows = FALSE`.
- The `print*Tibble()` and `printGSEAtibble()` helpers return a datatable: print it on its own pipe line, never through
  `%T>%`, which silently drops it. Under BiocStyle pass `BiocStyle = TRUE`.
- `prepareCor()` requires `cd` (colData whose `Group` levels are the contrast groups); `corTriplot()` draws a triangle
  only from three or more contrasts (two give a single cell).
- `genePlot(DGEdata =)` runs no tests: its brackets are the `DGEdata` rows with FDR < 0.05 and
  `abs(LogFC) > logFCcutoff`, stacked in `DGEdata` order (contrasts against the first `Group` level, in level order,
  nest by span).
- `makeSigTibble()`'s `contrast` is a factor of every `DGEdata` name in C order (`sort(method = "radix")`, which matches
  `printSigTibble()`'s `dplyr::arrange()` in any locale), so contrasts without DEGs stay downstream as zero-count rows.
- Plots drop filtered genes by their `NA` values, never by `rowData$retained`.
- Under tidySummarizedExperiment, a `dplyr::filter()` on a SummarizedExperiment that mixes rowData, colData and assay
  conditions errors, and an assay-only filter can return a tibble: filter after `summarise()` or on
  `tibble::as_tibble()` output.
- Gene set scores: `purrr::map()` the sets through `joinFeatures(score = TRUE, assay = ...)`, stack them with
  `purrr::list_rbind(names_to = "symbol")` and plot with `genePlot(genes = unique(.$symbol))`, passing the same `assay`
  to both (the score column is named after it). Wrap long set names at 20 characters.
- `dotPlot()` maps the group to `fill`, so a legend `guides()` call must set `fill` as well as `color`; a later
  `guides(fill = ...)` replaces the wrappers' guide, so give its `guide_legend()`
  `theme = ggplot2::theme(legend.byrow = TRUE)` to keep the row order. It jitters with
  `position_jitter(width = jitter, height = 0)`: `ggstripchart(jitter =)` alone also moves points vertically.
- Key caches on their inputs only (`cache_key()` of the objects that determine the result), never on package versions or
  function code; force a recompute when an upgrade could change the results.
- A ComplexUpset installed before ggplot2 4 makes `basicUpset()` fail with "Can't merge the `axis.text` theme
  element"; reinstall ComplexUpset.

## Vignette

- BiocStyle draws every figure at a fixed width (wider with `fig.wide = TRUE`, about half with `fig.small = TRUE`), so
  `fig.width` sets only the shape and the relative size of text and points: most figures are 6 in wide, and a wider one
  shrinks its text.
- Do not truncate `Term` before `GSEAplotList()`: `GSEAplot()` truncates only the axis labels (45 characters), so
  distinct terms stay apart.
- `comparisonList` pairs each `Group` level with the first level as the reference (`c(X, reference)`), in level order,
  as the skill builds it; chunks that take a contrast by position (`extract2(1)`) depend on that order.
- One `corTibble` (both orders of each pair) feeds `corTriplot()`; the correlation table and the 4-way grid use one
  order (`dplyr::filter(x < y)`). The grid is 2 x 2 at 10 x 10 in: each gg4way legend needs about 4.2 in, so three
  across clips it.
- The genes of interest are hardcoded because that chunk precedes the Differential expression section.
- References are pandoc citations (`[@key]`) to `vignettes/references.bib`, listed under the closing `# References`
  heading. Take each record from the package's `citation()` or from the paper's DOI record, and brace-protect names in
  titles (`{voom}`, `{RNA-seq}`) so citeproc's title case keeps them. R package entries carry no version
  (`note = {R package}`), and the text cites packages instead of linking them. `vignettes/basicBulk.csl` is pandoc's
  default style without `subsequent-author-substitute`, so a repeated author prints in full instead of "———"; keep its
  CC BY-SA notice.
- Check rendered output headlessly with `chromote` (on a server, add `--no-sandbox` with `chromote::set_chrome_args()`).

## Claude Code skill

- The skill's defaults are deliberate (for example `basicVolcano(style = "enhanced")` and `set.seed(100)` before the
  GSEA loop); do not change them to match package defaults.
- `fgsea` is not in DESCRIPTION: the skill reaches it only through `augere.gsea::runPrecomputed(methods = "fgsea")` and
  has the user install it.
- furrr multicore workers need `suppress.plots = TRUE`; a shared graphics device crashes with "PDF_endpage".
- knitr evaluates inline `r` in the YAML header before any chunk, so report titles are literal strings.
