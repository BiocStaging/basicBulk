# Package Reference

Use this reference for detailed function signatures, parameter names, default values, and return values when generating basicBulk code.

## `augere.de::runVoom()`

limma-voom differential expression on the raw counts of a SummarizedExperiment. Rules that apply to all three augere calls in this file:

- **No `...` in the signature.** An unknown or misspelled argument name is an error, never silently ignored.
- **Pass `assay = "counts"` to `runVoom()` and `runContrast()`.** The default `assay = 1` is the first assay, which is not guaranteed to hold the counts. `runPrecomputed()` has no `assay` argument.
- **`output.dir` always receives a generated `report.Rmd`** (the report is knit there, with the working directory switched into it for the run), even when nothing is saved. Pass a fresh `output.dir = tempfile()` per call so the analysis folder stays clean and parallel workers never share a directory, plus `save.results = FALSE` (the skill saves its own `RData/`).

| Parameter | Description | Default |
|-----------|-------------|---------|
| `x` | SummarizedExperiment with raw counts, genes in rows (the skill passes `se`, the `tidySE()` output) | required |
| `groups` | colData column with the factor of interest (the skill: `"Group"`) | required unless `design` + `contrasts` |
| `comparisons` | A comparison vector or a list of them. Length 2 = first group over second (sets the LogFC sign); two distinct element names = average the groups on each side; length 1 = test a covariate; length 3+ = ANOVA-like (`F` + `LogFC.*` columns, which the basicBulk DGE helpers do not handle) | required unless `design` + `contrasts` |
| `covariates` | Numeric colData columns, fitted as continuous terms | `NULL` |
| `block` | Nuisance factor columns, fitted as additive fixed effects | `NULL` |
| `subset.factor` / `subset.levels` | Keep only the samples whose `subset.factor` value is in `subset.levels` | `NULL` |
| `subset.groups` | Also drop samples of groups that no comparison names, so the fit only sees the compared groups | `TRUE`; skipped for covariate comparisons and for `design` + `contrasts` |
| `design` | Custom design: a formula (evaluated with `model.matrix()` on `colData(x)`), a function of `x`, or a numeric matrix with one row per (retained) sample | `NULL`; used **only when `contrasts` is also given**, and then `groups`, `comparisons`, `block`, `covariates` and `subset.groups` are ignored (write a fixed batch term into the design). `dc.block` and `subset.factor` still apply |
| `contrasts` | Custom contrast, or a list of them: a character expression for `limma::makeContrasts()` (names design columns), a named numeric vector/matrix, or a function of the design matrix | `NULL`; used only with `design` |
| `dc.block` | colData column passed to `limma::duplicateCorrelation()`; the consensus correlation goes into `lmFit()`. For correlated samples (e.g. repeated measures) that cannot be a `block` term | `NULL` |
| `robust` | Robust empirical Bayes in `eBayes()` | `TRUE` |
| `trend` | Variance trend in `eBayes()` (rarely needed with voom weights) | `FALSE` |
| `quality` | `voomWithQualityWeights()` (per-sample weights) instead of plain `voom()` | `TRUE` |
| `lfc.threshold` | LogFC threshold; above 0 the test switches to `limma::treat()` | `0` |
| `assay` | Assay holding the counts | `1`; pass `"counts"` |
| `row.data` | rowData columns to prepend to each result | `NULL`; leave it, `makeDGEdata()` joins the rowData itself |
| `metadata` | Extra metadata stored with saved results | `NULL` |
| `output.dir` | Folder for the generated report (and results, if saved) | `"voom"`; pass `tempfile()` |
| `save.results` | Write the results into `output.dir` | `TRUE`; pass `FALSE` |
| `author` / `dry.run` / `suppress.plots` | Report author (system user when `NULL`); write the report without running it (returns `NULL`); skip the report's plots | `NULL` / `FALSE` / `FALSE` |

Returns: `list(results, normalized)`.

- `results` is a list with one `DataFrame` per comparison or contrast. Rows are **all** genes of `x` in their original order (rownames = `rownames(x)`); columns `LogFC`, `AveExpr`, `t`, `PValue`, `FDR`, `B` (`B` is absent when `lfc.threshold > 0`). Genes removed by `filterByExpr()` stay in the table with `NA` statistics; `makeDGEdata()` drops them.
- Result names: ``"Increase in `X` over `Y`"`` for a length-2 comparison; the list name when `comparisons` is a named list; `"Contrast i"` for unnamed custom contrasts (the list name when `contrasts` is named). `makeDGEdata()` rewrites the default form to `"X vs Y"` by applying the same edits to every name: it removes `"Increase in "`, replaces ``"` over `"`` (between the backticked group names) with ``"` vs `"`` and removes the backticks. Helpers that split contrast names on `" vs "` (`prepareCor()`, the `genePlot()` statistics, `basicHeatmap(contrast = )`) need that form.
- `normalized` is a copy of `x` (sample-subset like the fit) with a `logCPM` assay of `edgeR::cpm(log = TRUE)` values over all genes (not the `tidySE()` voom values of the same name), rowData `retained`, and colData `lib.size`/`norm.factors`. Do not use it downstream; plots use the `tidySE()` assays on the full `se`.

## `augere.gsea::runContrast()`

Gene set tests that refit the voom model on the expression values, for **one** comparison or contrast per call. It shares `runVoom()`'s arguments `x`, `groups`, `covariates`, `block`, `subset.factor`, `subset.levels`, `subset.groups`, `design`, `dc.block`, `robust`, `quality`, `trend`, `assay`, `metadata`, `output.dir`, `author`, `dry.run`, `save.results` and `suppress.plots`, with the same meanings and the same rules (no `...`, `assay = "counts"`, `output.dir = tempfile()`, `save.results = FALSE`). It has no `lfc.threshold` or `row.data`, and no seed argument.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sets` | Gene sets: a named list (or `CharacterList`) of gene IDs matching `rownames(x)`. **Second positional argument** | required |
| `comparison` | **One** comparison vector, in the forms of a single `runVoom()` comparison minus the ANOVA-like one; a multi-comparison list errors | required unless `design` + `contrast` |
| `contrast` | **One** custom contrast, used with `design` (a multi-contrast list errors) | `NULL` |
| `methods` | Any of `"camera"`, `"fry"`, `"mroast"`, `"romer"` | `c("fry", "camera")` |
| `annotation` | `mcols(sets)` columns to copy into each result | `NULL`; keep it, `getReactome()` returns a plain list with no `mcols` |
| `output.dir` | Report folder | `"contrast"`; pass `tempfile()` |

Loop over the comparison list with `purrr::map()`/`purrr::imap()` and name the outer list yourself (see `references/design-patterns.md`). Call `set.seed(100)` once before the loop: only `mroast` and `romer` use random numbers, and without a seed each call draws its own.

Returns: a named list of `DataFrame`s, one per method (`$camera`, `$fry`, ...). Each row is a gene set (rownames = set names, in `sets` order) with `NumGenes` (set members that passed filtering), `Direction`, `PValue`, `FDR`, and method-specific extras. **`Direction` is `"Up"`/`"Down"` for camera and fry** (limma's own values; only romer is lowercase), which is what the basicBulk GSEA plots expect, so the default camera flow needs no Direction fix.

**The skill passes `methods = gsea_methods`**, defined once in the shared setup as `gsea_methods <- "camera"`. camera is a competitive test: it asks whether a set's genes are more differentially expressed than the genes outside it, adjusting for correlation between genes in the set. The package default `c("fry", "camera")` also runs fry, a fast self-contained test (are the set's genes differentially expressed at all?), but `tidyGSEAlist(GSEAmethod = gsea_methods)` keeps only the element it names, so fry would be computed and thrown away. To use fry instead, set `gsea_methods <- "fry"` so the call and the tidy step name the same single method.

## `augere.gsea::runPrecomputed()`

Pre-ranked gene set tests on precomputed per-gene statistics (the `DGEdata` tibbles) rather than a refit on expression. An opt-in replacement for `runContrast()`; the recipe is in `references/optional-sections.md` ("Stats-based GSEA"). It takes **no design arguments and no `assay` argument** (no `...` either), so it is design-agnostic: one `purrr::map()` over `DGEdata`.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `x` | data.frame/`DataFrame` of per-gene statistics with **rownames = gene IDs** matching `sets` (`tibble::column_to_rownames("ID")` on a `DGEdata` tibble). Pass **all** tested genes, never just the DEGs; rows with `NA` in the used columns are dropped | required |
| `sets` | As in `runContrast()` (second positional) | required |
| `rank.field` | Column the rank-based methods rank on | **no default, required**; the skill passes `"t"` (the signed moderated t) |
| `sign.field` | Signed effect column; when given, the rank statistic becomes `abs(rank.field) * sign(sign.field)` for every rank-based method | `NULL`; the skill passes `"LogFC"` (with `rank.field = "t"` this changes nothing, since t and LogFC share a sign) |
| `rank.sqrt` | Square-root the rank statistic's magnitude before re-signing (for F-like statistics) | `FALSE` |
| `methods` | Any of `"cameraPR"`, `"fgsea"`, `"geneSetTest"`, `"hypergeometric"`, `"goseq"` | `c("hypergeometric", "goseq", "fgsea", "cameraPR")`; always pass it. The skill uses `"cameraPR"` (default) or `"fgsea"` (needs the fgsea package) |
| `alternative` | `"mixed"`, `"up"`, `"down"` or `"either"` | `"mixed"`; pass **`"either"`**, the only setting that returns a `Direction` column |
| `signif.field` / `signif.threshold` | Significance column and cutoff that define "significant" genes | no defaults; read only by `"hypergeometric"`/`"goseq"`, omit for cameraPR/fgsea |
| `goseq.bias`, `fgsea.leading.edge`, `*.args` | Method-specific options | `"AveExpr"`, `FALSE`, `list()` |
| `annotation`, `metadata`, `author`, `output.dir`, `dry.run`, `save.results`, `suppress.plots` | As in `runContrast()` | `output.dir = "precomputed"` (pass `tempfile()`); `save.results = TRUE` (pass `FALSE`) |

fgsea and geneSetTest use random numbers: call `set.seed(100)` first.

Returns: a named list of `DataFrame`s, one per method (`$cameraPR`, `$fgsea`), one row per gene set (rownames = set names) with `NumGenes`, `PValue`, `FDR`, `Direction` (only with `alternative = "either"`), and extras (`ES`/`NES` for fgsea, `NumSig` for hypergeometric/goseq). **`Direction` is lowercase `"up"`/`"down"` for cameraPR and fgsea.** After `tidyGSEAlist(GSEAmethod = <method>)`, title-case it with `dplyr::mutate(Direction = stringr::str_to_title(Direction))` before plotting or caching: `topGSEAplot()` and `GSEAplotList()` factor `Direction` on `c("Up", "Down")`, so lowercase values become `NA`.

## Key `basicBulk` Functions

Shared conventions:

- DGE helpers take the `DGEdata` (from `makeDGEdata()`): a named list of tibbles, one per contrast, named `"X vs Y"`. Where it is an argument, the argument is also named `DGEdata`.
- Expression plots (`basicHeatmap()`, `genePlot()`, `geneHistogram()`, `PCAgrid()`, `getHVGs()`) and `joinFeatures()` default to `assay = "logCPM"`. The skill passes `assay = plotAssay` (`"logCPM"`, or `"ComBat"` with `adjustBatch = TRUE`); those assays are `NA` for filtered genes, which the plots drop and `joinFeatures(score = TRUE)` skips.
- The `print*` table helpers **return** a DT datatable without printing it: give each its own line or display chunk (`sigTibble %>% basicBulk::printSigTibble()`), never a `%T>%` tee, which discards the table. Their `BiocStyle = FALSE` default keeps the per-column filter boxes; leave it `FALSE` in html_document reports (`TRUE` is only for BiocStyle vignettes, where those boxes misalign).

| Function | Purpose |
|----------|---------|
| `tidySE(se, group = "Group", covariates, block, groups, adjustBatch = FALSE)` | Normalize a SummarizedExperiment (needs a `counts` assay and unique rownames) and return it with added assays over all genes and samples: `logCPM` (the genes are filtered with `edgeR::filterByExpr()`, TMM factors come from `edgeR::normLibSizes()`, and `limma::voomWithQualityWeights()` gives `E`; `NA` for filtered genes) and, when `adjustBatch = TRUE`, `ComBat` (voom of `sva::ComBat_seq()`-adjusted counts, needs a colData `Batch` column, `NA`-padded the same way). Adds rowData `retained` and colData `Group` (factor built from `group`; `groups` = named `c(old = "new")` vector that recodes and orders the levels, control first), `Sample` (= colnames), `lib.size`, `norm.factors`. Downstream functions also need a rowData `ID` column equal to the rownames and a `symbol` column |
| `makeDGEdata(se, voomList)` | Turn `runVoom()$results` into the `DGEdata`: per contrast, `ID` (from rownames) + the rowData columns (joined by `ID`) + the statistics, with `NA` (filtered) genes dropped; names rewritten from ``Increase in `X` over `Y` `` to `X vs Y` |
| `makeSigTibble(DGEdata, FDRcutoff = 0.05, logFCcutoff = 1)` | Stack the DEGs (`FDR < FDRcutoff` and `abs(LogFC) > logFCcutoff`) of all contrasts into one tibble with `contrast` (a factor whose levels are all the `DGEdata` names, so a contrast without DEGs keeps its level) and `direction` (`Up`/`Down`) columns |
| `makeSigCountTibble(sigTibble)` | Up/Down/total DEG counts per contrast, the input for `DEGbar()`. Every `contrast` level gets a row, so a contrast without DEGs shows as 0; to count only some contrasts, subset the `DGEdata` with `magrittr::extract()` before `makeSigTibble()`, which keeps their 0 rows |
| `printSigTibble(sigTibble, rearrange = TRUE, BiocStyle = FALSE)` | Searchable DT table of DEGs (contrast, symbol, description, LogFC, PValue, FDR; CSV/Excel buttons). `rearrange = FALSE` keeps the input order (e.g. an overlap table; a `dupe_count` column shows as `Hits`) |
| `DEGbar(sigTabyl, x = "contrast", addCounts = TRUE, flip = TRUE)` | Diverging Up/Down barplot of DEG counts from `makeSigCountTibble()`, with totals |
| `basicVolcano(DGEdata, FDRcutoff, logFCcutoff, lab = NULL, nrow, ncol, style = "basic")` | Volcano plots faceted by contrast, with each contrast's Up/Down DEG counts. `style = "enhanced"` (the skill's choice) is basicBulk's **own** ggplot2 styling (four significance colours plus cutoff lines), not the EnhancedVolcano package. `lab = TRUE` labels every DEG; a symbol vector labels those genes. The counts are a separate label layer keyed only by `contrast`, so re-faceting needs the count-layer step in `references/design-patterns.md` |
| `basicMA(DGEdata, FDRcutoff, logFCcutoff, pointShape, pointSize)` | MA plots (LogFC against average expression) faceted by contrast, DEGs coloured by direction |
| `basicHeatmap(se, sigTibble = NULL, geneSet = NULL, assay = "logCPM", contrast = NULL, pathways = FALSE, show_row_names = FALSE, draw = TRUE, ...)` | ComplexHeatmap of row Z-scores with columns split by `Group`. Genes come from `sigTibble$ID` (optionally one `contrast`, which also keeps only that contrast's two groups) or from `geneSet` (Ensembl IDs or symbols). `pathways = TRUE` takes a `makePathwayTibble()` tibble and splits the rows by pathway. Rows that are all `NA` after scaling (filtered or constant genes) are dropped; returns `NULL` with fewer than 2 genes. `...` goes to `ComplexHeatmap::Heatmap()` |
| `getHVGs(se, assay = "logCPM", n = 100, symbols = TRUE)` | The `n` most variable genes of `assay`, most variable first, as symbols (or Ensembl IDs with `symbols = FALSE`); needs MatrixGenerics |
| `basicUpset(DGEdata, directionSplit = TRUE, FDRcutoff, logFCcutoff, n_intersections = 40, sort_sets = "descending")` | ComplexUpset plot of DEG overlaps across contrasts; sets are split into `<contrast> Up`/`Down` unless `directionSplit = FALSE`. Underscores in contrast names become hyphens |
| `genePlot(se, genes, DGEdata = NULL, assay = "logCPM", ref.group = NULL, logFCcutoff = 1, plotType = "dot", nrow, ncol)` | Expression of each gene (`genes` = symbols) by `Group`, one facet per gene; `plotType = "box"` for boxplots. Passing `DGEdata` adds significance brackets for contrasts with `FDR < 0.05` and `abs(LogFC) > logFCcutoff`, which needs contrast names `X vs Y` whose `X`/`Y` are `se$Group` levels and an SE `se`. `se` can also be a tibble of gene set scores (stacked `joinFeatures(score = TRUE)` output with the set names in `symbol`). Rows that are `NA` in the assay are dropped, with a warning naming any gene or set left with no values, and an error when none are left |
| `geneHistogram(se, geneSymbol = NULL, assay = "logCPM")` | Histogram of per-gene mean expression, faceted by `Group`; `geneSymbol` marks that gene's mean with a dashed line |
| `PCAgrid(se, assay = "logCPM", removeVar = 0.1, pal = NULL, ...)` | PCAtools PC1/PC2 biplot (coloured by `Group`) beside the top PC1/PC2 gene loadings. Genes with `NA` values are dropped before the PCA; `removeVar` removes that fraction of lowest-variance genes; `...` reaches `PCAtools::biplot()` (e.g. `colby`, `shape = "Batch"`) |
| `prepareCor(corTibble, DGEdata, cd, format = "wide")` | LogFC correlation (all genes, via `gg4way::getCor()`) for each contrast pair in `corTibble` (columns `x`, `y` = `DGEdata` names; `name` = `"<x label> vs <y label>"`), returned as a `corrr` `cor_df`. Pipe the `corTibble` in as the first argument. **`cd` is required**: pass `cd = SummarizedExperiment::colData(se)`. The labels are factored against `levels(cd$Group)`, so they must match the stripped contrast tokens |
| `corTriplot(corRes)` | Lower-triangle (`corrr::shave()`) blue-white-red `r` plot of a `prepareCor()` result; a triangle needs at least 3 contrasts |
| `getReactome(organism)` | Reactome pathways from MSigDB C2:CP:REACTOME via msigdbr, as a plain named list (`REACTOME_*` names, unversioned Ensembl gene IDs). `organism` is `"Homo sapiens"` or `"Mus musculus"` only (anything else errors); mouse gets the human sets mapped by orthology |
| `tidyGSEAlist(GSEAlist, GSEAmethod = "camera")` | Take the `GSEAmethod` element of each contrast's result and return a named list of tibbles sorted by FDR: `Database` (`REACTOME`), `Term` (set name without the prefix, `_` as spaces), `NumGenes`, `Direction`, `PValue`, `FDR` (plus method extras). `GSEAmethod` must be one method that was run |
| `topGSEAplot(GSEAlist, n_terms = 2)` | Dot plot of the top `n_terms` Up and Down terms of each contrast (FDR < 0.05) across all contrasts; size = -log10 FDR, fill = `Direction` |
| `GSEAplotList(GSEAlist, n_terms = 5)` | Named list of per-contrast barplots of the top Up and Down terms |
| `printGSEAtibble(GSEAlist, BiocStyle = FALSE)` | Searchable DT table of every term with FDR < 0.05 across contrasts |
| `makePathwayTibble(pathways, gene.sets)` | Member genes of the chosen sets as a tibble (`Pathway` factor, `ID`) for `basicHeatmap(pathways = TRUE)`. `pathways` are full set names (`REACTOME_` + underscores), e.g. rebuilt from `Term` with `paste0("REACTOME_", .) %>% stringr::str_replace_all(" ", "_")`. Genes in several chosen sets repeat, so apply `dplyr::distinct(ID, .keep_all = TRUE)` before the heatmap |
| `joinFeatures(se, features, symbol = FALSE, score = FALSE, long = TRUE, assay = "logCPM")` | The assay values of `features` (Ensembl IDs, or symbols with `symbol = TRUE`) joined to the colData, one row per gene and sample. `score = TRUE` averages the features into one score per sample, skipping genes that are `NA` in the assay; the value column is named after `assay`. Stack several sets with `purrr::map()` + `purrr::list_rbind(names_to = "symbol")` and plot them with `genePlot()` |
| `printPlots(.x, .y, depth = 2)` | Write a markdown heading `.y` and knit-print plot `.x`; with `purrr::iwalk()` in a `results = "asis"` chunk under a `{.tabset}` header it makes one tab per plot |
| `printHeatmaps(.x, .y, depth = 2)` | The same for heatmaps |
| `cache_key(...)` | Content hash of the objects a cache depends on |
| `save_cache(file, key, ...)` | Save named objects together with their key |
| `load_if_fresh(file, key)` | Load the cache and return `TRUE` only when its stored key equals `key`; otherwise load nothing and return `FALSE`, so the caller recomputes |
| `assert_fresh(artifact, source, regen_cmd = NULL)` | Stop if `artifact` is missing or older than `source` (file timestamps) |
| `fmt_e(x)` | Format numbers as `"1.2e-04"` for inline R; a raw number in scientific notation renders as math |
