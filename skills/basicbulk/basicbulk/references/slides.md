# basicBulk Analysis Slide Deck

Read this reference before generating a Quarto RevealJS slide deck from a completed basicBulk
analysis. The deck summarises findings from an already-rendered analysis report whose caches live
in `RData/`.

## YAML Front Matter

```yaml
---
title: "PROJECT — STUDY TITLE"
subtitle: "Study type · Cell/tissue · Organism"
author: "AUTHOR"
date: today
format:
  revealjs:
    theme: simple
    slide-number: true
    width: 1920
    height: 1080
    margin: 0.04
    fig-align: center
    embed-resources: true          # single self-contained HTML — no companion _files/ dir
    include-in-header:
      text: |
        <style>
        /* CSS goes here — injected into <head>, not a slide (see CSS section) */
        </style>
    title-slide-attributes:
      data-background-color: "#f5f5f5"
execute:
  echo: false
  warning: false
  message: false
---
```

`embed-resources: true` makes the output a single portable HTML with all images base64-encoded.
This is the **only** output format.

## Setup Chunk: Staleness Guards + Live Values

All headline statistics must trace to a live read of the `RData/` caches. Never hardcode a count
or p-value in slide prose. The setup chunk loads the caches, guards their freshness, and derives
every number that appears on slides.

````r
```{r setup, include=FALSE}
# Working directory: set explicitly so the chunk can run from anywhere.
root_dir <- "/path/to/PROJECT"
setwd(root_dir)

suppressPackageStartupMessages({ library(basicBulk); library(magrittr) })

# ── Staleness guards ──────────────────────────────────────────────────────────
# settings.RData is written in Setup when se0 is loaded — it is the
# canonical "raw data is this fresh" anchor. Compare all artifacts against it.
# assert_fresh() and fmt_e() are basicBulk functions — no inline definitions needed.
basicBulk::assert_fresh("RData/se.RData", "RData/settings.RData",
                        'rmarkdown::render("report.Rmd")')
basicBulk::assert_fresh("RData/voom.RData", "RData/settings.RData",
                        'rmarkdown::render("report.Rmd")')
basicBulk::assert_fresh("RData/GSEA.RData", "RData/settings.RData",
                        'rmarkdown::render("report.Rmd")')

# ── Load caches ───────────────────────────────────────────────────────────────
invisible(lapply(c(
  "RData/settings.RData",   # group, groups, organism, outliers, se0, PROJECT_TITLE, STUDY_DESCRIPTION
  "RData/se.RData",         # se
  "RData/voom.RData",       # voomList, DGEdata, sigTibble
  "RData/GSEA.RData"        # gene.sets, GSEAlist
), load, envir = globalenv()))

# ── Pre-format p-values ───────────────────────────────────────────────────────
# Raw numerics in scientific notation (2e-07) get knitr's math formatter and
# render as broken 2^{-7}. Always pre-format with basicBulk::fmt_e() before embedding in prose.

# ── Derive live values ────────────────────────────────────────────────────────
n_samples   <- nrow(SummarizedExperiment::colData(se0))
n_contrasts <- length(DGEdata)

deg_counts <- sigTibble %>%
  dplyr::count(contrast, direction)
# Example: n_down_contrast1 <- deg_counts %>% dplyr::filter(...) %>% dplyr::pull(n)

# GSEA top term per contrast
top_gsea <- GSEAlist %>%
  dplyr::bind_rows(.id = "contrast") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::group_by(contrast) %>%
  dplyr::slice_min(FDR, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(contrast, Term, Direction, FDR) %>%
  dplyr::mutate(FDR = purrr::map_chr(FDR, basicBulk::fmt_e))
```
````

## CSS — via `include-in-header`

Embed CSS using the `text` sub-key of `include-in-header` in the YAML front matter. This
injects the `<style>` block directly into `<head>` and is picked up by `embed-resources: true`,
so the output is a **single portable HTML** with no companion files needed.

**Do not** use `css: |` in YAML (URL-encodes the block as a broken `<link>` href — nothing
loads). **Do not** use a `{=html}` block in the document body (RevealJS treats pre-`##` content
as a blank slide). **Do not** use an external `.css` file (not embedded by `embed-resources:
true` — sharing just the HTML loses all styles).

```yaml
format:
  revealjs:
    embed-resources: true
    include-in-header:
      text: |
        <style>
        .reveal h2 { font-size: 1.35em; margin-top: 0; margin-bottom: 0.25em; }
        .reveal h3 { font-size: 1.05em; color: #444; margin-top: 0.3em; margin-bottom: 0.15em; }
        .reveal p  { font-size: 0.88em; line-height: 1.4; margin: 0.25em 0; }
        .reveal li { font-size: 0.88em; line-height: 1.4; margin: 0.15em 0; }
        .reveal ul { margin: 0.2em 0 0.2em 1.2em; padding: 0; }
        .reveal .small { font-size: 0.72em; }
        .reveal .takehome {
          font-size: 0.88em; font-weight: bold; color: #1a3a6a;
          border-left: 4px solid #1a3a6a; padding: 0.2em 0.6em;
          margin-bottom: 0.4em; line-height: 1.35;
        }
        .two-col {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 1.5em;
          align-items: start;
        }
        .two-col-wide {
          display: grid;
          grid-template-columns: 1.1fr 0.9fr;  /* wider left col */
          gap: 1.5em;
          align-items: start;
        }
        table { font-size: 0.76em; border-collapse: collapse; }
        td, th { padding: 0.15em 0.45em; }
        .reveal section > * { margin-top: 0; }
        </style>
```

Adjust font sizes / gap / grid proportions as needed. Add more `grid-template-columns` variants
(e.g. `2fr 1fr`) for other asymmetric layouts.

## Slide Structure

### Deck opening

Open with two context slides before any findings, so the audience knows the system and trusts the
data before they see a result. Recommended default — omit them only if the user explicitly wants
a findings-only deck.

1. **Overview** — what the study is. Pull the study title and description recorded in Setup
   (`PROJECT_TITLE`, `STUDY_DESCRIPTION`, loaded from `settings.RData`), the design variables (the
   `Group` factor and any dose / genotype / timepoint columns), per-group sample sizes (live from
   `colData(se0)`, never hardcoded), and the model-system / treatment context (cell line or tissue,
   `organism`, treatment duration). This is a **context slide — no take-home**; it orients the
   audience, it doesn't conclude.
2. **Sample QC** — did the data come out clean. Show the report's PCA figure (`plots/QC/PCA.png`;
   when `adjustBatch = TRUE`, the default vs batch-adjusted pair the report saved as
   `plots/QC/{tidyName}.png`, side by side in a two-column layout) and give the per-sample library
   sizes as numbers, read live from `se$lib.size` — no new plot. Its take-home is **read from the
   actual PCA / library-size result and the QC Outlier Scan outcome** (SKILL.md), not assumed: if no
   sample stands out, say so; if any sample was flagged, name it and state whether it was dropped
   (the `outliers` vector) or kept. Never default to "all samples pass QC" — that sentence is wrong
   the moment one doesn't, and it's exactly the unverified claim the take-home rule below warns
   against.

### Section headers
Use `##` for all slides (flat structure). Section-break slides use background colour via attribute:

```markdown
## Main Analysis {background-color="#003366" color="white"}
## Summary {background-color="#1a3a1a" color="white"}
```

**Do not use `#` top-level headers** — they create nested RevealJS sections, and `#/N` fragment
addressing in Chrome headless (used for screenshots) only reaches section roots, not content slides
nested inside them.

### Two-column layout

```markdown
::: {.two-col}

::: {}
Left column content — text, table, bullet list.
:::

::: {}
![Caption.](plots/DGE/figure.png){height=600}
:::

:::
```

### Take-home on every finding slide

Every finding slide (not section headers, not overview slides) leads with a bold take-home:

```markdown
::: {.takehome}
Take-home: One-sentence plain-language message the audience keeps after 30 seconds.
:::
```

Two guardrails:

1. The take-home summarises what the slide's own data shows — never a stronger claim. Verify
   direction, scope, and precision against the actual contrast/model result before writing it.
2. **No safety or risk verdicts.** Do not write "safe", "risky", "low risk", "of concern", or
   "manageable". Summarise the biology and evidence, name what warrants further investigation,
   and let the scientist judge. This applies to DEG findings, GSEA enrichment, or any
   other result — the model can be wrong and a safety claim in a slide deck misleads.

### Summary slides

Focus on findings only. Do not include a files / outputs table — the audience is interested in
the biology, not the file names.

## Figure Sizing

**Check the source PNG dimensions before setting `{height=}` or `{width=}`** (all three print width
then height; use whichever is available — `sips` is macOS-only, `identify`/`png::readPNG` work on Linux):

```sh
sips -g pixelWidth -g pixelHeight plots/DGE/volcano_grid.png   # macOS
identify -format '%w %h\n' plots/DGE/volcano_grid.png          # Linux (ImageMagick)
```

```r
d <- png::readPNG("plots/DGE/volcano_grid.png");  dim(d)[2:1]  # width, height — anywhere R runs
```

At 1920×1080 with `margin: 0.04`, the usable body is ~1843×980 px, minus title (~65px) and
take-home (~55px) ≈ **860 px body height**.

For a two-column layout, each column is ~900 px wide.

| Aspect ratio | Constraint | Set |
|---|---|---|
| Landscape (W > H, e.g. 1800×1200 = 1.5:1) | width-clamped at ~900px | set `height=600` — comes out 900px wide |
| Portrait (H > W, e.g. 2100×2400 = 0.875:1) | height-limited | set `height=700`, width = 0.875×700 = 613px — fits column |
| Very wide (2.25:1, e.g. 2700×1200) | width-clamps to ~900px, height collapses to 400px | set `width=900` and accept the constrained height |

For **full-width** single-image slides, max usable is ~1843px wide; set `height=` to respect the
860px body budget.

If bumping `{height=}` does nothing, the image is width-clamped — widen the column or regenerate
the PNG with a taller aspect ratio.

**Two layout traps that cost the most iterations:**

- **Content placed *below* a tall grid clips off the bottom.** The slide is a fixed 1080 px, so anything below a tall `.two-col`/grid can silently run off the bottom edge. Fix by moving that content **into a column of the grid** (or trimming the text) — **not** by shrinking the figure. Have the screenshot-reading subagent (see Visual QC Sweep) explicitly inspect the bottom ~10 % of every slide for cut-off content.
- **Do not add a deck-wide `.reveal img { height: auto }`** to "fix" sizing — it overrides per-image `{height=}` attributes and makes tall images overflow. Size tall/composite figures by `{height=}` (never a blanket CSS rule).

**Shared figures — the report produces the PNGs the deck embeds.** Deck slides `![](plots/…png)` figures that a **report chunk `ggsave`s** to that same path, so the two are coupled: (a) restyling a shared plot in the report (a new legend, changed dims) silently changes the deck too — re-render **both**, and note the deck's HTML keeps a stale baked-in copy until it is itself re-rendered; (b) if the deck needs a *different* variant (e.g. no legend, to save slide space), give it its **own** `ggsave` to a distinct path rather than editing the shared one; (c) embed height-only (`{height=}`, per the aspect table above) so the shared PNG's aspect is preserved wherever it lands.

## Rendering

The Quarto CLI is a prerequisite (see SKILL.md "Environment Setup"). With `quarto` on `PATH`:

```sh
quarto render PROJECT_slides.qmd
```

quarto auto-detects R from `PATH`, so when the R on `PATH` is the one with basicBulk installed no
`QUARTO_R` override is needed (if a chunk can't find a package, run `quarto check knitr` — it should
report the libPaths holding basicBulk).

On a **laptop (macOS)** without quarto on `PATH`, use the quarto bundled with RStudio Desktop:

```sh
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render PROJECT_slides.qmd
```

Either way, outputs a single self-contained `PROJECT_slides.html`. No companion `_files/` directory.

## Visual QC Sweep (required after any layout change)

After any change to figure placement, `{height=}`, `{width=}`, or slide splitting: render →
screenshot every slide → dispatch a subagent that can read images to read the PNGs and return a
punch-list.

**Do NOT Read PNGs in the main thread** — images accumulate in context and cannot be evicted
until `/compact`.

Use **decktape** (not raw Chrome headless `#/N`). Chrome headless fires the screenshot before
RevealJS lazy-loads images, producing blank slides regardless of `--virtual-time-budget`. Decktape
properly waits for RevealJS to initialize.

Note: by default `npx decktape` downloads its own Google Chrome for Testing (~630 MB) on first run via
Puppeteer. To reuse an existing Chrome and skip that download, see Prerequisites below; otherwise the
one-time download is harmless.

```sh
# 1. Screenshot all slides as a PDF
NPM_CONFIG_CACHE=/tmp/npm-cache npx --yes decktape reveal \
  "file:///abs/path/to/PROJECT_slides.html" \
  /tmp/PROJECT_deck.pdf \
  --size 1920x1080

# 2. Convert PDF pages to PNGs (Ghostscript)
gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r150 \
  -dFIXEDMEDIA -dPDFFitPage \
  -dDEVICEWIDTHPOINTS=1920 -dDEVICEHEIGHTPOINTS=1080 \
  -sOutputFile=/tmp/slide_%02d.png \
  /tmp/PROJECT_deck.pdf

# 3. Dispatch a subagent that can read images to read the PNGs and return a punch-list
```

For the subagent prompt: give it the absolute paths to each PNG, the expected layout per slide,
and ask it to flag **overflow** (content cut at bottom), **too-small** (labels illegible at 3 m),
or **OK**. Fix issues and re-sweep until clean — a clean deck is the expected outcome.

**Prerequisites.** The sweep needs **Node/npm** (`npx`), a **Chrome** binary, and **Ghostscript**
(`gs`). To point decktape at an existing Chrome and skip the ~630 MB Chrome-for-Testing download, use
`PUPPETEER_SKIP_DOWNLOAD=true` + `--chrome-path`:

```sh
NPM_CONFIG_CACHE=/tmp/npm-cache PUPPETEER_SKIP_DOWNLOAD=true \
npx --yes decktape reveal \
  "file:///abs/path/to/PROJECT_slides.html" \
  /tmp/PROJECT_deck.pdf \
  --size 1920x1080 \
  --chrome-path /path/to/google-chrome
```

On a sandboxed/container node where Chrome's sandbox can't initialize, also pass
`--chrome-arg=--no-sandbox --chrome-arg=--disable-gpu --chrome-arg=--disable-dev-shm-usage` (harmless if
not needed). Without an existing Chrome, the plain command above works as-is (decktape downloads its own
Chrome once). Either way the `embed-resources: true` HTML is portable, so the deck can also be rendered on
one machine and swept on another.

## Gotchas

- **Three wrong ways to add CSS; one right way.** `css: |` YAML URL-encodes the block as a
  broken `<link>` href — nothing loads. A `{=html}` body block before the first `##` creates a
  blank slide in RevealJS. An external `.css` file is not embedded by `embed-resources: true` —
  sharing just the HTML loses all styles. **The only correct approach:** `include-in-header:
  text: |` in the YAML (see CSS section above) — injected into `<head>`, embedded, no blank slide.
- **`#` top-level sections nest slides.** Chrome headless `#/N` only reaches the section root,
  not content slides nested within. Use `##` for everything; add `{background-color=...}` for
  section-break colour.
- **Inline p-values break without `fmt_e()`.** A raw numeric like `4.02e-07` passed through
  knitr renders as 4.02^{-7}. Pre-format with `basicBulk::fmt_e()`
  (`formatC(signif(x, 2), format = "e", digits = 1)`) — no inline definition needed.
