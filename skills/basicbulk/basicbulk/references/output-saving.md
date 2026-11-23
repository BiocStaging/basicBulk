# Output Saving

Saving outputs is on by default (template mode saves at every milestone; guide mode at each checkpoint; opt-out only). The setup chunk creates the directory structure (SKILL.md Section 1); use the patterns below for RData objects and plots.

## Directory Structure

The setup chunk creates:

```
RData/
tables/
plots/
plots/QC/
plots/DGE/
plots/DGE/ContrastVolcanos/
plots/DGE/4way/
plots/heatmaps/
plots/GSEA/
plots/scores/
```

For subdirectories not in this tree (e.g., `plots/target_histograms/`), use `create.dir = TRUE` in `ggplot2::ggsave()`.

## RData Save Points

Save objects with `save()` at each workflow milestone. These enable loading intermediate state without re-running upstream code.

### After Setup (Section 1)

```r
save(group, groups, organism, outliers, se0, PROJECT_TITLE, STUDY_DESCRIPTION,
     file = "RData/settings.RData")
```

`group`/`groups` (the grouping column and its levels) are defined in the Setup chunk (SKILL.md Section 1). Include any extra setup variables relevant to the study (e.g., `palette`).

### After Normalization (Section 2)

```r
save(se, file = "RData/se.RData")
```

### After DGE (Section 4)

```r
save(voomList, DGEdata, sigTibble,
     file = "RData/voom.RData")
```

Also include `sigCountTibble`, `contrastList` or `designTibble` when they exist (`sigCountTibble` is assigned only in the design-pattern downstream visualization in `references/design-patterns.md`; SKILL.md Section 5 pipes `basicBulk::makeSigCountTibble()` inline, leaving no object to save).

### After Heatmaps (Section 6)

```r
save(col, top_annotation, file = "RData/heatmapAnno.RData")
```

Only when heatmap annotation objects are built manually with ComplexHeatmap.

### After GSEA (Section 10)

The GSEA-over-Reactome step (`augere.gsea::runContrast()` looped per contrast over `basicBulk::getReactome(organism)`) is the slowest part of a render, and `rmarkdown::render()` re-runs every chunk from scratch — there is no knitr chunk caching. So `RData/GSEA.RData` is worth caching. But a GSEA cache that can go **stale** is dangerous: if an upstream input changes (an outlier sample dropped, a covariate / `block` / `dc.block` changed, a contrast edited, the input counts revised, or the gene sets changed) and the cache is reused anyway, the Pathways section silently reports last run's results. This can happen with no one "watching" — during refinement the analyst iterates with the report for many turns without re-invoking the skill, and a plain mtime check misses in-place edits and out-of-order re-saves. So guard the cache by **content hash**, not timestamp, and make it safe-by-default.

Use the package guard (`basicBulk::cache_key()` / `load_if_fresh()` / `save_cache()`). The Pathways compute chunk becomes:

```r
run_gsea <- FALSE   # committed default: cache-load unless an input changed; the guard below decides.
## gsea_methods (= "camera") is defined ONCE in setup (references/design-patterns.md Shared Template Setup) —
## the SAME object passed to the loop's methods =, to tidyGSEAlist(GSEAmethod = gsea_methods), and to this
## cache key, so a fry opt-in can't drift between the computation and the cache. Do not redefine it here.

## key = the identity-determining inputs that change the GSEA result: cheap, cross-session-stable
## METADATA plus the REALIZED counts values. User inputs (an .rds, a counts CSV) can change under the
## same sample and gene IDs, so colData + rownames alone would miss a revised counts table; as.matrix()
## realizes the values (dimnames included) so a DelayedArray can't hash on its backing pointer. You MUST
## also include three things a colData hash alone misses:
##  1. The contrast spec object YOUR design pattern defines (references/design-patterns.md):
##       Pattern 1      -> comparisonList
##       Pattern 2      -> designTibble
##       Patterns 3/4   -> comparisonList  AND the subset/split column NAMES (subset.factor / split vars)
##       custom design  -> design + contrastList
##  2. The model-argument SELECTIONS: block, covariates, and dc.block are column-NAME choices you pass to
##     runContrast()/runVoom() (as are the Pattern 3/4 subset/split columns in 1). dc.block is a
##     duplicateCorrelation argument, NOT a design column, so hashing colData values alone would miss a
##     switch like dc.block = "Batch" -> "Animal". Pass the actual argument values (a string or NULL each).
##  3. The GSEA methods SELECTION passed to runContrast() (canonical idiom = "camera"). A switch to
##     "fry"/c("fry","camera") changes the tested result, so key on the SAME gsea_methods object the loop
##     uses (from setup) — a hardcoded literal here could silently drift from a fry switch in the loop.
##     A multi-method run tidies with a separate scalar gsea_tidy_method; add that object to the key too.
## Key on those INPUTS only — NO automatic package/engine-version token. cache_key() hashes just its arguments
## (rlang::hash), so reuse across an upgrade assumes unchanged GSEA behaviour: review the relevant engine /
## dependency / result-tidying changes and force a recompute (run_gsea <- TRUE) when results may change or the
## impact is uncertain. Manual invalidation is the maintenance responsibility that replaces a version token.
gsea_key <- basicBulk::cache_key(
  SummarizedExperiment::colData(se),                     # the sample set — dropping an outlier changes this
  rownames(se),                                          # the gene set
  as.matrix(SummarizedExperiment::assay(se, "counts")),  # the realized counts values (dimnames included)
  comparisonList,                                        # <- contrast spec for THIS pattern (see 1 above; P2: designTibble; custom: design + contrastList)
  BLOCK, COVARIATES, DC_BLOCK,                           # model-arg selections (column names or NULL) — catches a block/covariate/dc.block switch; P3/4: add the subset/split column names
  gsea_methods,                                          # GSEA methods — same object passed to the loop's methods = (see 3 above)
  gene.sets                                              # the gene-set collection (hashed as an object; no version token — see the key note above)
)

if (run_gsea || !basicBulk::load_if_fresh("RData/GSEA.RData", gsea_key)) {
  GSEAlist <- ... # set.seed(100) + the runContrast() loop (per the detected design pattern)
  basicBulk::save_cache("RData/GSEA.RData", gsea_key, gene.sets = gene.sets, GSEAlist = GSEAlist)
}
```

This key is the contract for the expression engine (`augere.gsea::runContrast()`). The opt-in stats-based engine (`augere.gsea::runPrecomputed()`) keys on its own inputs — see `references/optional-sections.md` → Stats-based GSEA.

`load_if_fresh()` loads `gene.sets`/`GSEAlist` **only** when the stored key matches the current inputs; on any mismatch (or a cache with no stored key, or a missing file) it returns `FALSE` and loads nothing, so the chunk recomputes — stale can never be silently reused. `run_gsea <- TRUE` is a manual force-recompute escape hatch on top of the guard, never the only defense. The guard lives in the report code, so it fires for anyone who renders (analyst at the console, an unattended re-render) regardless of whether the skill is loaded. (The slide deck only *reads* the cache — it can't recompute — so it uses the lighter `basicBulk::assert_fresh()` mtime check + `load()`; see `references/slides.md`.) Chain it: any downstream cache keyed on GSEA output should key on the *freshly-loaded* objects so staleness can't propagate. This is the same guard the optional-section gates use (`references/optional-sections.md` → Gating expensive optional sections). When iterating on sections *below* Pathways, this guarded load already yields `gene.sets`/`GSEAlist` from the cache without re-running the GSEA loop.

## Plot Saving Patterns

All ggplot-based plots save into the `plots/` subdirectories using `%T>%` for inline saves or `purrr::iwalk()` for batch saves. Heatmaps are the exception — `basicBulk::basicHeatmap()` returns a `ComplexHeatmap`, so it saves via a `png`/`pdf` graphics device (see "Heatmaps: graphics device" below), not `ggsave`.

**Dimensions are not authored here.** `references/plot-dimensions.md` is the single source. The examples below use `WIDTH`/`HEIGHT` placeholders for the save `width`/`height`; when generating, replace them with the plot's dimensions from `plot-dimensions.md` (a plot's save dims equal its chunk `fig.width`/`fig.height`).

### Single Plot: `%T>%` ggsave

The plot expression evaluates, saves to disk, and continues the pipe:

```r
se %>%
  basicBulk::PCAgrid(assay = plotAssay,
                     colby = "Group",
                     selectLab = NA,
                     alpha = 0.8) %T>%
  ggplot2::ggsave("plots/QC/PCA.png", plot = ., width = WIDTH, height = HEIGHT)
```

Wrap complex ggplot chains in `(...)` before `%T>%`:

```r
(sigCountTibble %>%
   basicBulk::DEGbar()) %T>%
  ggplot2::ggsave("plots/DGE/DEG_barplot.png", plot = ., width = WIDTH, height = HEIGHT)
```

### Faceted or Grid Plot: `%T>%` ggsave

```r
DGEdata %>%
  purrr::set_names(names(.) %>% stringr::str_remove(" vs .*")) %>%
  basicBulk::basicVolcano(style = "enhanced") %T>%
  ggplot2::ggsave("plots/DGE/volcano_grid.png", plot = ., width = WIDTH, height = HEIGHT)
```

### Batch Save: `purrr::iwalk()` with glue paths

For plot lists where each element saves to a named file:

```r
GSEAlist %>%
  basicBulk::GSEAplotList() %T>%
  purrr::iwalk(~ ggplot2::ggsave(glue::glue("plots/GSEA/{.y}_GSEA.png"),
                                 plot = .x, width = WIDTH, height = HEIGHT))
```

### Per-Contrast Save: single-element sublist with `purrr::imap()`

Build a list of individual plots from `DGEdata[.x]` (a one-element sublist), then batch-save:

```r
DGEdata %>%
  seq_along() %>%
  purrr::map(~ DGEdata[.x] %>%
               basicBulk::basicVolcano(style = "enhanced")) %>%
  magrittr::set_names(names(DGEdata) %>% stringr::str_remove(" vs .*")) %T>%
  purrr::imap(~ ggplot2::ggsave(glue::glue("plots/DGE/ContrastVolcanos/{.y}_volcano.png"),
                                plot = .x, width = WIDTH, height = HEIGHT))
```

`basicVolcano()` takes the full `DGEdata`; use `DGEdata[.x]` (a one-element list) for per-contrast plots, not `DGEdata[[.x]]` (a bare tibble).

### Heatmaps: graphics device, not `ggsave`

`basicBulk::basicHeatmap()` returns a `ComplexHeatmap` object, not a ggplot, so it cannot be saved with `ggplot2::ggsave()`. Save heatmaps by opening a graphics device, building the heatmap (which draws to the open device), then closing it with `dev.off()`. Match `width`/`height` to the heatmap chunk's `fig.width`/`fig.height` (see `references/plot-dimensions.md`).

The all-DEGs heatmap → one `png`:

```r
png("plots/heatmaps/DEG heatmap.png", res = 300, units = "in", width = WIDTH, height = HEIGHT)
allHeatmap <- se %>%
  basicBulk::basicHeatmap(sigTibble = sigTibble,
                          assay = plotAssay)
dev.off()
```

The per-contrast `heatList` → one multi-page `pdf`:

```r
pdf("plots/heatmaps/extended complex heatmaps.pdf", width = WIDTH, height = HEIGHT)
heatList <- names(DGEdata) %>%
  purrr::set_names(stringr::str_remove(., " vs .*")) %>%
  purrr::map(~ se %>%
               basicBulk::basicHeatmap(se = .,
                                       sigTibble = sigTibble,
                                       contrast = .x,
                                       assay = plotAssay))
dev.off()
```

The pathway heatmap (Pathways → Heatmap tab) → one `png` the same way: `png("plots/heatmaps/pathway heatmap.png", res = 300, units = "in", width = WIDTH, height = HEIGHT)` around its build, then `dev.off()`.

Wrap the **build** chunk (the `allHeatmap`/`heatList` assignments) in the device — not a separate save step — so the objects are still available for the print chunk that displays them as tabs. Give this build chunk `include=FALSE` (not `fig.show="hide"`): the device + `basicHeatmap` draw emits "Saving … image"-style device messages that `fig.show="hide"` leaves visible, whereas `include=FALSE` suppresses them while still running the chunk so `allHeatmap`/`heatList` remain in scope for the print chunk.

### New Subdirectory: `create.dir = TRUE`

When saving into a path not created by the setup chunk:

```r
c("GENE_SYMBOL_1", "GENE_SYMBOL_2") %>%
  purrr::set_names() %>%
  purrr::map(~ se %>%
               basicBulk::geneHistogram(geneSymbol = .x,
                                        assay = plotAssay)) %T>%
  purrr::iwalk(~ ggplot2::ggsave(glue::glue("plots/target_histograms/{.y}_histogram.png"),
                                 plot = .x, width = WIDTH, height = HEIGHT,
                                 create.dir = TRUE))
```

## Plot Save Path Reference

Default file paths by section. **Dimensions are not repeated here** — `references/plot-dimensions.md` is the single source for both chunk `fig.*` and save `width`/`height`. `N` = number of contrasts.

| Section | Plot | Save Path |
|---------|------|-----------|
| PCA | PCAgrid | `plots/QC/PCA.png` |
| PCA | Batch PCAs (assays `c("logCPM", "ComBat")`) | `plots/QC/{tidyName}.png` |
| DGE | DEGbar | `plots/DGE/DEG_barplot.png` |
| DGE | basicVolcano grid | `plots/DGE/volcano_grid.png` |
| DGE | basicVolcano per-contrast | `plots/DGE/ContrastVolcanos/{.y}_volcano.png` |
| DGE | basicMA grid | `plots/DGE/MA_grid.png` |
| DGE | basicHeatmap (all DEGs) | `plots/heatmaps/DEG heatmap.png` (via `png`) |
| DGE | basicHeatmap (per-contrast list) | `plots/heatmaps/extended complex heatmaps.pdf` (via `pdf`) |
| Overlaps | basicUpset (UpSet tab) | `plots/DGE/UpSet.png` |
| Overlaps | Overlap Expression (genePlot) | `plots/DGE/overlap expression.png` |
| DGE | gg4way (per-pair tabs) | `plots/DGE/4way/{.y}.png` |
| DGE | corTriplot (Overview) | `plots/DGE/LogFC Correlations.png` |
| DGE | gg4way combined grid (optional, small N) | `plots/DGE/4way/combined_grid.png` |
| Genes | genePlot | `plots/DGE/{gene} expression.png` |
| Genes | geneHistogram | `plots/target_histograms/{.y}_histogram.png` |
| GSEA | topGSEAplot | `plots/GSEA/top GSEA.png` |
| GSEA | GSEAplotList | `plots/GSEA/{.y}_GSEA.png` |
| GSEA | basicHeatmap (pathway expression, Pathways → Heatmap tab) | `plots/heatmaps/pathway heatmap.png` (via `png`) |
| Gene Set Scores | genePlot (gene set scores) | `plots/scores/gene set scores.png` |

## When to Save

- **Template mode:** Always include RData saves and `%T>%` ggsave calls by default. Every generated `.Rmd` should save objects at each milestone and plots at each visualization step.
- **Guide mode:** Save after each checkpoint (normalization/PCA, DGE, GSEA) — both the RData save and plot saves for that section.
- **Opt-out only:** Only omit saves if the user explicitly asks not to save.
