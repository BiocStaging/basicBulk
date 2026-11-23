# Optional Analysis Sections

Read this reference before generating optional sections such as LogFC Correlations (gg4way + corTriplot), gene plots, or gene set scores.

Use `plotAssay` for optional plotting sections. When batch adjustment is enabled, `plotAssay` should already be `"ComBat"`; otherwise it should be `"logCPM"`. Pass it as `assay = plotAssay` to the basicBulk plotting function on `se`. Statistics always use `se` with `assay = "counts"`.

## Gating expensive optional sections (toggle + cache + freshness guard)

Some optional sections are slow or need remote services, so recomputing them on every render is wasteful — but a cache that can go stale is worse than no cache. Gate such a section so it renders safely anywhere and can never silently reuse an out-of-date result. The Pathways GSEA cache (`references/output-saving.md` → After GSEA) and Stats-based GSEA below both use this pattern; reuse it for any new expensive section.

- **Toggle.** A committed `run_<section> <- FALSE` that the analyst flips to `TRUE` to force a (re)generation. When `FALSE`, the section loads its cache instead of recomputing whenever that cache is fresh.
- **Cache keyed on its inputs.** Save with `basicBulk::save_cache("RData/<section>.RData", key, ...)` where `key <- basicBulk::cache_key(<the actual invalidating input objects>)`, and load with `fresh <- basicBulk::load_if_fresh("RData/<section>.RData", key)`. `load_if_fresh()` loads **only** when the stored key matches the current inputs; otherwise it returns `FALSE` and loads nothing (fail-safe: recompute). This survives stale caches even when file mtimes lie (in-place edits, out-of-order saves) — it is the same mechanism as the GSEA guard in `references/output-saving.md`.
- **Render-safe gating.** These sections need only in-memory data or a network service, so they can recompute anywhere: `if (run_<section> || !basicBulk::load_if_fresh(...)) { compute; basicBulk::save_cache(...) }` — recompute when the cache is stale, since you can. Downstream chunks `eval` on the loaded object's existence (`eval = exists("<obj>")`) so the section renders empty and error-free when its inputs/cache are absent.
- **Chain the guards.** Build each cache's `key` from its *freshly-computed upstream* objects (never from objects loaded out of another cache), so staleness cannot propagate through a cache-of-a-cache.

## LogFC Correlations

This section examines how the per-gene logFC profiles agree across contrasts. It pairs two views: a `gg4way` plot for each pair of contrasts (the per-pair detail), and — when there are enough contrasts to make a matrix meaningful — an overview correlation triangle (`basicBulk::corTriplot()`) that condenses every pairwise *r* into one figure. Auto-include whenever there are two or more contrasts.

R Markdown header: `# LogFC Correlations` — add `{.tabset .tabset-pills}` only when the section has `##` children (i.e. at N ≥ 3, per the routing below). Any chunk that calls `purrr::iwalk(basicBulk::printPlots)` needs `results="asis"`.

Routing by number of contrasts (`N = length(DGEdata)`):

| `N` | What to generate |
|-----|------------------|
| 2 | A single `gg4way` plot. Drop `{.tabset .tabset-pills}` from the header (there are no `##` children — no `## Overview plot`, no per-pair tabs) and no correlation triangle, since a 2×2 correlation matrix is just one number. Put the *r* in the title instead (`gg4way::getCor()`). |
| ≥ 3 | `## Overview plot` tab (the `corTriplot()` triangle) **plus** one `## <x> vs <y>` tab per contrast pair (the loop below). |

### Per-pair gg4way (the detail view)

This is the default detail view — give each pair of contrasts its own `gg4way` plot. Do not collapse the pairs into a single combined figure by default; one plot per pair is what lets the reader see each comparison at full size.

At **N ≥ 3**, loop the whole pairwise set into tabs. Capture the named list as `fourWayList` first (the combined-grid extra below reuses it), then walk it:

```r
fourWayTibble <- tidyr::crossing(
  x = names(DGEdata),
  y = names(DGEdata)
) %>%
  dplyr::filter(x < y) %>%
  dplyr::mutate(name = glue::glue("{stringr::str_remove(x, ' vs .*')} vs {stringr::str_remove(y, ' vs .*')}"))

fourWayList <- fourWayTibble %>%
  purrr::pmap(~ with(list(...),
                     gg4way::gg4way(DGEdata,
                                    x = x,
                                    y = y,
                                    logFC = "LogFC",
                                    FDR = "FDR"))) %>%
  purrr::set_names(fourWayTibble$name)

fourWayList %>%
  purrr::iwalk(basicBulk::printPlots)
```

`dplyr::filter(x < y)` keeps one entry per unordered pair (drops both the self-pairs and the mirror duplicates), so the tab count is `choose(N, 2)`. When the user specifies custom cutoffs, pass `logFCcutoff = logFC_cutoff` and `FDRcutoff = FDR_cutoff` to `gg4way::gg4way()`.

At **N = 2** there is a single pair, so skip the tabset and the `printPlots` loop — render the one plot directly and put the *r* in its title:

```r
gg4way::gg4way(DGEdata,
               x = names(DGEdata)[1],
               y = names(DGEdata)[2],
               logFC = "LogFC",
               FDR = "FDR") %>%
  {. + ggplot2::labs(title = glue::glue("{stringr::str_remove(names(DGEdata)[1], ' vs .*')} vs {stringr::str_remove(names(DGEdata)[2], ' vs .*')} (r = {gg4way::getCor(.)})"))}
```

### Overview correlation triangle (`## Overview plot`, N ≥ 3)

`basicBulk::prepareCor()` runs `gg4way::getCor()` over every contrast pair and returns a `corrr` `cor_df`; `basicBulk::corTriplot()` shaves it to a lower triangle and renders a blue–white–red `r` heatmap. The correlation uses the logFC of **all** genes, not just the DEGs.

Keep this chunk **compact** — it is a small overview by default (dimensions in `references/plot-dimensions.md`), and at low/moderate contrast counts you should not enlarge it. Only when a genuinely high contrast count makes the cells unreadable should you resize it up — or better, split the overview by control (one smaller triangle per control, below) or rely on the per-pair gg4way tabs.

```r
tidyr::crossing(x = names(DGEdata),
                y = names(DGEdata)) %>%
  dplyr::mutate(name = glue::glue("{stringr::str_remove(x, ' vs .*')} vs {stringr::str_remove(y, ' vs .*')}")) %>%
  dplyr::filter(x != y) %>%
  basicBulk::prepareCor(DGEdata = DGEdata,
                        cd = SummarizedExperiment::colData(se)) %>%
  basicBulk::corTriplot()
```

`filter(x != y)` drops only the self-comparisons (a contrast against itself, trivially `r = 1`); `prepareCor()` builds the full symmetric matrix and `corTriplot()` keeps the lower triangle.

**`cd` must carry a `Group` column whose levels match the stripped contrast tokens.** `prepareCor()` factors the contrast names (after stripping `" vs .*"`) against `levels(cd$Group)`; `cd` has no default, so always pass it (normally `SummarizedExperiment::colData(se)`). In a compound-group design — where the contrasts are built on a combined factor (e.g. `sampleGroup` like `"DrugA 50nM"`) but `colData(se)$Group` holds only the bare molecule (`"DrugA"`) — the stripped tokens fall outside `levels(Group)`, every cell goes `NA`, and `corTriplot()` errors on the empty matrix. Pass the matching factor explicitly: `... %>% basicBulk::prepareCor(DGEdata = DGEdata, cd = SummarizedExperiment::colData(se) %>% as.data.frame() %>% dplyr::mutate(Group = COMPOUND_GROUP_FACTOR))`.

### Split by control (ask only when it applies)

A contrast's logFC is measured relative to its control, so the control structure shapes how the overview triangle should be read:

- When several contrasts **share one control**, their logFCs sit on a common baseline — correlations among them compare treatment effects directly.
- When contrasts have **different controls**, a single triangle correlates across differing baselines as well as treatment effects.

Neither is universally the "right" or "wrong" view — which one the analyst wants is project-dependent. So detect the structure and **ask** rather than guess. Derive each contrast's control from the text after `" vs "`:

```r
controls <- names(DGEdata) %>% stringr::str_extract("(?<= vs ).*")
splittable <- names(table(controls))[table(controls) >= 2]
```

If `length(splittable) >= 2` (two or more distinct controls each anchoring at least two contrasts — e.g. a dose-matched design with one control per dose), ask the user whether to split the overview into **one triangle per control** (each built from only that control's contrasts, via the same `prepareCor()` recipe on the subset) or keep a single combined triangle. Otherwise — a single shared control (one triangle), each contrast with its own control (no group has ≥ 2 members), or a lone control with one straggler — don't ask; emit the single overall triangle. A per-control split needs ≥ 2 contrasts in the group to form a triangle, which is why the `>= 2` guard matters.

**When the split (or the single-baseline limitation) is in play, say so in the narrative** — a one-line note under the section header that the correlations are read within/across a shared control — and carry the same caveat into the slide take-home if a deck is built (see `references/slides.md`). A correlation figure with a hidden baseline confound is easy to over-read; the note is what keeps the reader honest.

### Optional combined grid (small N only)

The per-pair tabs are the default. As an extra, when there are a few contrasts (`3 <= N <= 4`, so 3–6 panels) the pairwise `gg4way` plots can also be arranged into one combined figure for a single at-a-glance comparison. This reuses `fourWayList` from the per-pair loop above, which only exists at N ≥ 3 — so the combined grid applies only there (at N = 2 there is a single plot and nothing to combine). Above N = 4 the panels get too cramped to read, so stay with per-tab only. Save under `plots/DGE/4way/combined_grid.png`. The grid uses `patchwork`, a basicBulk Suggests package.

**Can the panels share one legend? First work out what varies between them.** `gg4way` colours each panel by *that panel's own x/y contrasts*, on a deterministic four-level scale you can rely on:

| gg4way level | colour |
|---|---|
| Not Significant | `grey80` |
| Significant in Both | `mediumblue` |
| Significant in the **x** contrast | `firebrick` |
| Significant in the **y** contrast | `forestgreen` |

So each panel's legend reads "Significant in `<x-contrast>`" / "Significant in `<y-contrast>`" — a meaning tied to *that* panel. Whether a shared legend is honest depends on whether that meaning is identical in every panel:

- **Standard pairwise grid (`fourWayList`): each panel is a different pair of contrasts, so `firebrick`/`forestgreen` mean something different per panel.** One shared legend would mislabel every panel but the first, which is misleading — keep **per-panel legends**:

```r
patchwork::wrap_plots(fourWayList, nrow = NROW)   # NROW = ceiling(choose(N, 2) / 2) for two columns
```

  Leaving `patchwork::plot_layout(guides = "collect")` off (or on) is safe here: `collect` only merges guides that render *identically*, so mismatched per-panel legends stay separate rather than silently collapsing into a wrong one.

- **Split / crossed grid: every panel shares a common comparison and only some *other* factor varies between them** (built from a split/crossed `DGEdata`, see Pattern 4 in `references/design-patterns.md`). Here a shared legend is correct *once you relabel each panel's colour scale to that common meaning*, so all panels render the same guide and patchwork collects them into one:

```r
splitFourWayList %>%
  purrr::imap(~ .x +
                ggplot2::labs(title = glue::glue("{.y} (r = {gg4way::getCor(.x)})")) +
                ggplot2::scale_color_manual(
                  values = c("Not Significant"          = "grey80",
                             "Significant in Both"      = "mediumblue",
                             "Significant in X_LEVEL"   = "firebrick",     # THIS panel's real x-level string
                             "Significant in Y_LEVEL"   = "forestgreen"),  # THIS panel's real y-level string
                  labels = c("Not Significant", "Significant in Both",
                             "Significant in X_COMMON", "Significant in Y_COMMON"),  # shared display meaning
                  drop = FALSE, name = NULL,
                  guide = ggplot2::guide_legend(nrow = 2, override.aes = list(size = 3)))) %>%
  patchwork::wrap_plots(nrow = NROW) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
```

Prefer **`patchwork` + `plot_layout(guides = "collect")`** over `ggpubr::ggarrange(common.legend = TRUE)`: `ggarrange` renders solid black filler boxes flanking the merged legend whenever the legend is narrower than the combined figure (a short or multi-row legend on side-by-side panels). `basicBulk::printPlots` (via `knitr::knit_print`) handles patchwork objects, so it drops straight into a `results="asis"` tabset loop. If you must keep `ggarrange`, make the shared legend `nrow = 1` / full-width so it fills the strip and the black boxes don't appear.

**Why the relabel works, and the two shapes it takes.** patchwork collects by the *rendered* guide (title + labels + glyphs), so panels merge only when their `labels=` match. Name the `values` vector to *this* panel's real gg4way level strings (`scale_color_manual` maps colour→meaning by name, not position) and set `drop = FALSE` so every panel carries all four levels even when one is empty. Then move whatever *differs* between panels out of the legend — two shapes:

- **Map to a higher-level / common axis** — collapse the per-panel contrast labels to a constant shared meaning (e.g. "Significant in `X_COMMON`" / "Significant in `Y_COMMON`" in every panel).
- **Drop the varying token** — when panels are the same comparison across a varying factor, strip that factor from `labels=` so the legend shows only the shared comparison, and carry it in each panel's title and axis labels (`ggplot2::labs()`) instead.

The varying factor is data-dependent — dose is only one example; it may be timepoint, region, cell type, assay, genotype, and so on. **Read it off the contrast/panel names**: the token shared across panels is the legend's common meaning, the token that differs is what moves to the title/axis. Work this out from the contrast structure yourself; only when it is genuinely ambiguous which factor is shared vs varying — or whether a shared legend belongs at all — ask the user rather than guess.

## Gene Plots

Use these only when the user provides genes of interest or requests expression histograms.

```r
# Genes of interest
se %>%
  basicBulk::genePlot(genes = GENES_OF_INTEREST,
                      plotType = "dot",
                      assay = plotAssay)

# Expression histogram
se %>%
  basicBulk::geneHistogram(geneSymbol = "GENE_SYMBOL",
                           assay = plotAssay)
```

When the user specifies custom cutoffs, pass `logFCcutoff = logFC_cutoff` to `basicBulk::genePlot()`.

## Stats-based GSEA (`runPrecomputed`)

The default GSEA engine is `augere.gsea::runContrast()` (SKILL.md §10 Pathways, looped per design pattern in `references/design-patterns.md`). **Only when the user asks for stats-based / pre-ranked GSEA**, swap that compute step for `augere.gsea::runPrecomputed()`. This is an opt-in *replacement*, not an additional section: it produces the same tidied `GSEAlist` (once `Direction` is title-cased, below), so everything below it in §10 (`topGSEAplot()` / `GSEAplotList()` / `printGSEAtibble()`) runs unchanged. Don't run both engines unless the user explicitly wants the comparison.

The two engines differ only in the compute step — pick on the distinguishing axis, *what the test ranks on*:

| | `runContrast()` (default) | `runPrecomputed()` (opt-in) |
|---|---|---|
| Input | `se` with `assay = "counts"` — **refits the voom model** per contrast | the already-built `DGEdata` (per-contrast stats tibbles; no `assay` argument) |
| What it tests on | the original expression values | the **precomputed per-gene statistics** — pre-ranked by `rank.field` (`"t"`; the `DGEdata` carries a signed moderated `t`) |
| Methods | `c("fry", "camera")` | `cameraPR` (**default**) or `fgsea` |
| Design-dependence | **per design pattern** (needs `groups`/`comparison` or `design`/`contrast`, plus `block`) | **design-agnostic** — one `purrr::map` over `DGEdata`, no design args |
| Tidy | `tidyGSEAlist()` (`GSEAmethod = "camera"`) | `tidyGSEAlist(GSEAmethod = GSEAMETHOD)` |

`cameraPR` is the **default** stats method because it is the pre-ranked analogue of the default `camera` and tracks it closely; `fgsea` is consistently more liberal, so offer it only when the user asks.

Because the stats engine consumes `DGEdata`, **build `DGEdata` before this step** — it is the input. In the default flow `DGEdata` and `GSEAlist` are built in the pattern-specific GSEA/DGE code (`references/design-patterns.md`), not the shared setup that defines `se`/`plotAssay`/`gene.sets`; for the stats engine just order `DGEdata` first, then run `runPrecomputed` over it.

```r
GSEAMETHOD <- "cameraPR"   # pre-ranked default; "fgsea" only on request

## Cache under the same freshness guard as the default engine (SKILL.md §10): key on INPUTS only — the DGEdata
## it consumes, the gene sets, and the method — with NO package/engine-version token (same rule as the default
## engine; force a recompute on an engine change that could alter results, see output-saving.md -> After GSEA).
run_gsea <- FALSE
gsea_key <- basicBulk::cache_key(DGEdata, gene.sets, GSEAMETHOD)

if (run_gsea || !basicBulk::load_if_fresh("RData/GSEA.RData", gsea_key)) {
  set.seed(100)
  GSEAlist <- DGEdata %>%
    purrr::map(~ .x %>%
                 dplyr::filter(!is.na(ID)) %>%
                 tibble::column_to_rownames("ID") %>%
                 augere.gsea::runPrecomputed(sets = gene.sets,   # getReactome(organism), same collection as the default engine
                                             methods = GSEAMETHOD,
                                             rank.field = "t",
                                             sign.field = "LogFC",
                                             alternative = "either",
                                             output.dir = tempfile(),
                                             save.results = FALSE),
               .progress = TRUE) %>%
    basicBulk::tidyGSEAlist(GSEAmethod = GSEAMETHOD) %>%
    purrr::map(~ dplyr::mutate(.x, Direction = stringr::str_to_title(Direction)))
  basicBulk::save_cache("RData/GSEA.RData", gsea_key, gene.sets = gene.sets, GSEAlist = GSEAlist)
}
```

A few things keep this correct:

- **`alternative = "either"`** is what makes `runPrecomputed` return the `Direction` column that `tidyGSEAlist()` and the dot/bar plots split on — don't drop it. `runPrecomputed` returns it lowercase (`up`/`down`), while `topGSEAplot()` / `GSEAplotList()` key on `Up`/`Down`, so title-case it right after `tidyGSEAlist()`, before the cache and the plots.
- **Ranking is by `rank.field = "t"`, which must be passed** (`runPrecomputed` has no default for it). The `DGEdata` carries a signed moderated `t`. With a `sign.field` supplied, the pre-rank statistic becomes `abs(t) * sign(LogFC)` — but since `t = LogFC / SE` with `SE > 0`, `sign(t) == sign(LogFC)` for every gene, so this equals the raw signed `t`. `sign.field = "LogFC"` is kept by convention; it is a no-op here for that reason, not because of `rank.sqrt` (which only swaps the magnitude to `sqrt(abs(t))` when `TRUE`). Pass `column_to_rownames("ID")` so the gene IDs are the rownames `runPrecomputed` matches the gene sets against.
- **`tidyGSEAlist(GSEAmethod = ...)` must name the same method you ran** (`"cameraPR"`/`"fgsea"`) — it extracts that element from each result. The bare `tidyGSEAlist()` defaults to `"camera"`, which is the contrast engine's method and would error here.
- The tidied output has the **identical shape** to the contrast engine (`Database, Term, NumGenes, Direction, PValue, FDR`) once `Direction` is title-cased, which is why the §10 downstream needs no changes. Save the `GSEAlist` to `RData/GSEA.RData` exactly as the default engine does (see `references/output-saving.md`) so the slide deck and downstream re-renders load it the same way.

### Report prose

Mirror the default Pathways intro (Report Prose Style — descriptive, methodological, tool linked, no interpretation) but capture the one thing that differs: the test ranks on the **precomputed differential-expression statistics**, not the expression values, and uses **cameraPR** rather than camera. The default contrast intro reads:

> Gene Set Enrichment Analysis (GSEA) was performed using the REACTOME pathways gene set (MSigDB C2:CP:REACTOME, obtained with `basicBulk::getReactome()`) through [augere.gsea](https://bioconductor.org/packages/augere.gsea). The testing utilized the original gene expression values, which improves statistical performance. A competitive test using camera was utilized, as this is the [recommended approach for testing a large number of gene sets](https://support.bioconductor.org/p/103504/#103572).

The stats-engine intro keeps the first and last sentences and swaps the middle clause for the pre-ranked statistic:

> Gene Set Enrichment Analysis (GSEA) was performed using the REACTOME pathways gene set (MSigDB C2:CP:REACTOME, obtained with `basicBulk::getReactome()`) through [augere.gsea](https://bioconductor.org/packages/augere.gsea). The testing was performed on the precomputed per-gene differential expression statistics (the moderated t-statistic), ranking all tested genes for a pre-ranked enrichment test. A competitive pre-ranked test using cameraPR was utilized, as this is the [recommended approach for testing a large number of gene sets](https://support.bioconductor.org/p/103504/#103572).

If the user runs `fgsea` instead, name fgsea in the last sentence in place of cameraPR.

## Gene Set Scores

A gene set score is one value per sample for each chosen pathway: `basicBulk::joinFeatures(score = TRUE)` averages
the pathway's genes in the assay, skipping the genes that are `NA` there (those removed by `filterByExpr`), and
`basicBulk::genePlot()` plots the stacked scores, faceted by pathway with `Group` on the x axis. Pass
`assay = plotAssay` to both calls: `joinFeatures()` names the score column after the assay, and `genePlot()` reads that
column, so the batch-adjusted values are used when `adjustBatch = TRUE`. A pathway with no expressed genes is dropped
with a warning.

Ask the user if they want scores for specific pathways. Either reuse the vector chosen for the Pathways heatmap, name them
explicitly, or — when the top enrichments are dominated by a **redundant family** (e.g. many near-duplicate
ribosome/translation pathways crowding out distinct themes) — select a diverse, non-redundant set programmatically:

```r
# Greedy, redundancy-aware pathway selection. Ranks significant pathways
# (FDR < 0.05) by best FDR across contrasts, then keeps a pathway only if its
# gene-membership Jaccard overlap with EVERY already-kept pathway is below
# jaccard_cutoff. This keeps one representative per redundant family and fills the
# remaining slots, up to n, with the most distinct themes. Returns FEWER than n if
# too few distinct pathways exist (diversity is preferred over hitting the count;
# raise jaccard_cutoff or name pathways explicitly if you need more). n is the
# analysis's choice. REACTOME names are assumed (basicBulk::getReactome()).
selectDiversePathways <- function(GSEAlist, gene.sets, n = 6, jaccard_cutoff = 0.2) {
  ranked <- GSEAlist %>%
    dplyr::bind_rows(.id = "contrast") %>%
    dplyr::filter(FDR < 0.05) %>%
    # group_by Direction keeps a pathway significant in both directions as two
    # candidates; the duplicate is harmlessly rejected later by self-Jaccard = 1.
    dplyr::group_by(Term, Direction) %>%
    dplyr::summarise(minFDR = min(FDR), .groups = "drop") %>%
    dplyr::arrange(minFDR) %>%
    dplyr::mutate(setName = Term %>%
                    paste0("REACTOME_", .) %>%
                    stringr::str_replace_all(" ", "_"))

  # The display Term -> collection setName round-trip should be lossless. Warn
  # (rather than silently drop) if a significant pathway fails to map back to
  # gene.sets, so a future collection change can't quietly skew the selection.
  unmatched <- setdiff(ranked$setName, names(gene.sets))
  if (length(unmatched) > 0) {
    cli::cli_warn(c("!" = "selectDiversePathways dropped {length(unmatched)} pathway name(s) not found in gene.sets:",
                    unmatched %>% purrr::set_names(rep("*", length(.)))))
  }

  ranked <- ranked %>%
    dplyr::filter(setName %in% names(gene.sets))

  memberList <- gene.sets[ranked$setName] %>%
    as.list() %>%
    purrr::map(as.character)

  jaccard <- function(a, b) {
    if (length(union(a, b)) == 0) return(0)
    length(intersect(a, b)) / length(union(a, b))
  }

  # Greedy walk down the FDR-ranked list; sequential dependency on `kept` plus an
  # early break make an explicit loop clearer here than purrr::reduce.
  kept <- character(0)
  for (sn in ranked$setName) {
    if (length(kept) == 0 ||
        max(purrr::map_dbl(kept, ~ jaccard(memberList[[sn]], memberList[[.x]]))) < jaccard_cutoff) {
      kept <- c(kept, sn)
    }
    if (length(kept) >= n) break
  }
  kept
}

scorePathways <- selectDiversePathways(GSEAlist, gene.sets, n = 3)
```

Or reuse the heatmap's pathways (`scorePathways <- pathways`), or name them explicitly:

```r
scorePathways <- c("REACTOME_PATHWAY_NAME_1",
                   "REACTOME_PATHWAY_NAME_2")
```

Then plot and save the scores (dimensions in `references/plot-dimensions.md`):

```r
if (length(scorePathways) > 0) {                             # a selection with no significant pathways would otherwise fail
  (scorePathways %>%
     magrittr::extract(gene.sets, .) %>%
     purrr::set_names(~ .x %>%
                        stringr::str_remove("^REACTOME_") %>%
                        stringr::str_replace_all("_", " ") %>%
                        stringr::str_wrap(20)) %>%
     purrr::map(~ basicBulk::joinFeatures(se, features = .x, score = TRUE, assay = plotAssay)) %>%
     purrr::list_rbind(names_to = "symbol") %>%
     basicBulk::genePlot(genes = unique(.$symbol),
                         assay = plotAssay,
                         plotType = "dot") +
     ggplot2::labs(y = "Gene set score")) %T>%
    ggplot2::ggsave("plots/scores/gene set scores.png",
                    plot = .,
                    width = WIDTH,
                    height = HEIGHT)
}
```

The panel labels drop the `REACTOME_` prefix and wrap at 20 characters, and the panels keep the order of
`scorePathways`. Keep the parentheses around the plot: `%T>%` binds tighter than `+`, so without them `ggsave()`
receives only the `labs()` labels object and fails (see `references/code-style.md`).

## Secondary-Factor Aesthetics (crossed designs)

When `Group` is a crossed `PRIMARY_FACTOR × SECONDARY_FACTOR` (e.g. molecule × dose), per-sample
plots read best with the **primary factor on x + colour** and the **secondary factor encoded by
point shape** — rather than putting the full crossed `Group` (all combinations) on the x-axis.
Add `PRIMARY_FACTOR` and `SECONDARY_FACTOR` as their own colData columns in setup, alongside the
combined `Group`. Applies to PCA and to `genePlot` (gene plots and gene set scores; it wraps `dotPlot`) below.

**PCA** — `PCAgrid()` forwards `shape` to `PCAtools::biplot()`; coerce the shape column with
`forcats::as_factor()` first (shape can't be continuous; it keeps numeric levels in numeric order and an existing
factor's order):

```r
se %>%
  dplyr::mutate(SECONDARY_FACTOR = forcats::as_factor(SECONDARY_FACTOR)) %>%
  basicBulk::PCAgrid(assay = plotAssay,
                     colby = "PRIMARY_FACTOR",
                     shape = "SECONDARY_FACTOR",
                     selectLab = NA,
                     alpha = 0.8)
```

**`genePlot()` (it wraps `dotPlot`)** fixes the point shape internally, so a bare
`+ ggplot2::aes(shape = ...)` won't take effect. Drop the fixed param on the point layer, then
map shape and supply fillable shapes (21 = circle, 24 = triangle). For `genePlot`, re-point
`Group` to the primary factor so its x / colour follow it. Do not pass `DGEdata =` when `Group`
is re-pointed: its stats branch parses the contrast names against `Group`.

```r
p <- se %>%
  dplyr::mutate(Group = PRIMARY_FACTOR,
                SECONDARY_FACTOR = forcats::as_factor(SECONDARY_FACTOR)) %>%   # shape can't be continuous
  basicBulk::genePlot(genes = GENES_OF_INTEREST,
                      assay = plotAssay)

p$layers[[1]]$aes_params$shape <- NULL             # free the fixed shape so the mapping applies

(p +
   ggplot2::aes(shape = SECONDARY_FACTOR) +
   ggplot2::scale_shape_manual(values = c(21, 24)) +   # one fillable shape (21-25) per SECONDARY_FACTOR level
   ggplot2::guides(fill = ggplot2::guide_legend(theme = ggplot2::theme(legend.byrow = TRUE),
                                                override.aes = list(shape = 21))))
```

For gene set scores, stack the scores as in Gene Set Scores, then re-point `Group` on the score tibble (it carries the
colData columns) before `genePlot()`:

```r
p <- scorePathways %>%
  magrittr::extract(gene.sets, .) %>%
  purrr::set_names(~ .x %>%
                     stringr::str_remove("^REACTOME_") %>%
                     stringr::str_replace_all("_", " ") %>%
                     stringr::str_wrap(20)) %>%
  purrr::map(~ basicBulk::joinFeatures(se, features = .x, score = TRUE, assay = plotAssay)) %>%
  purrr::list_rbind(names_to = "symbol") %>%
  dplyr::mutate(Group = PRIMARY_FACTOR,
                SECONDARY_FACTOR = forcats::as_factor(SECONDARY_FACTOR)) %>%   # shape can't be continuous
  basicBulk::genePlot(genes = unique(.$symbol),
                      assay = plotAssay)

p$layers[[1]]$aes_params$shape <- NULL

(p +
   ggplot2::aes(shape = SECONDARY_FACTOR) +
   ggplot2::scale_shape_manual(values = c(21, 24)) +
   ggplot2::labs(y = "Gene set score") +
   ggplot2::guides(fill = ggplot2::guide_legend(theme = ggplot2::theme(legend.byrow = TRUE),
                                                override.aes = list(shape = 21))))
```

`scale_shape_manual(values = c(21, 24))` covers a two-level secondary factor; extend the vector
for more levels. The point geom is `layers[[1]]` in current basicBulk — if a future version prepends
a layer, target the `GeomPoint` layer instead. A future package option (a `shape =` argument on
`genePlot()`) would remove the `aes_params` surgery entirely; not available yet.
