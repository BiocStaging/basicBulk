---
name: basicbulk
description: Use for bulk and pseudobulk RNA-seq analysis and reporting with the basicBulk R package (DGE via augere.de, GSEA via augere.gsea) on a SummarizedExperiment or a counts table plus sample sheet — normalization, PCA, differential gene expression (DGE, voom-limma), volcano/MA plots, DEG heatmaps, overlaps/UpSet, LogFC correlations, GSEA/Reactome pathway enrichment and gene-set scores, and an optional Quarto slide deck. Always use for anything naming basicBulk, or for explicit /basicbulk:basicbulk invocations (template mode, or /basicbulk:basicbulk guide). For a generic bulk RNA-seq or DGE request that does not mention basicBulk, still activate, then ask the user whether to use the basicBulk workflow before proceeding. Do NOT use for single-cell/scRNA-seq, ChIP-seq, proteomics, generic plotting or file conversion, or dashboards unrelated to a DGE workflow.
---

# basicBulk Bulk RNA-seq Analysis

Guide bioinformaticians through bulk RNA-seq analysis using basicBulk with the augere.de / augere.gsea pipelines. Two modes: template generation and interactive step-by-step guidance.

## Reuse First

When writing ANY analysis code, do NOT reinvent the wheel — take it for a drive. Follow this reuse priority:

1. **FIRST, follow this basicBulk skill and use the package's own functions** (`basicBulk::`, `augere.de::`, `augere.gsea::`) — the documented idioms almost always cover the task.
2. **THEN fall back to the package vignette** for the real-world idiom: `vignette("basicBulk")` (if the installed package has no built vignette, read `vignettes/basicBulk.Rmd` in the package source).

Coding STYLE must always follow this skill's own style rules — the inline conventions in this `SKILL.md` plus the extended `references/code-style.md` ("Code Style Reference").

**Propagate this rule (do it early, proactively).** Early in any basicBulk analysis, make sure the project's own `CLAUDE.md` (the analysis memory file — see **Analysis CLAUDE.md** below) carries this reuse-first priority, creating the `CLAUDE.md` if none exists. Build on and augment whatever logic is already there — never replace it. The point: the project `CLAUDE.md` auto-loads every session, so landing the rule there makes it persist for later follow-ups even when this skill isn't re-invoked, and for any user — not just the author. Permission behavior: if you are ALREADY creating the project `CLAUDE.md`, or have ALREADY been granted permission to edit it this session, just ADD the rule directly — do NOT ask again. Only ask first when editing the project `CLAUDE.md` would otherwise require permission you don't already have.

## Activation

This skill assumes the basicBulk workflow (augere.de / augere.gsea). If the request is an explicit `/basicbulk:basicbulk` invocation, or names basicBulk, proceed directly. If instead the request is a **generic** bulk RNA-seq or DGE request that does **not** mention basicBulk, first ask the user whether they want to use the basicBulk workflow before proceeding — do not assume it.

## Modes

- **Template:** `/basicbulk:basicbulk` or `/basicbulk:basicbulk template` - Collect inputs, inspect data via R, and generate a complete `.Rmd`.
- **Guide:** `/basicbulk:basicbulk guide` - Pair-program through the analysis interactively. Rules:
  - **Execution:** Claude runs each chunk via R, shows the output, and adapts if errors occur.
  - **Pacing:** Group related chunks into logical units based on dependencies (e.g., `augere.de::runVoom` + `makeDGEdata` + `makeSigTibble` together since they form a pipeline).
  - **Checkpoints:** Pause and ask the user after three key points: (1) normalization/PCA, (2) DGE results, (3) GSEA results. Flow freely within those blocks.
  - **Navigation:** The user can say "skip to pathways", "go back to PCA", or "redo that with different parameters" at any checkpoint.
  - **Adaptation:** If the user changes a parameter mid-guide (e.g., different control group, add a covariate), regenerate all downstream chunks.

## Audience

Bioinformaticians who know R and Bioconductor but may not know every basicBulk function. Explain which functions to use and why. Assume R and Bioconductor familiarity.

## Reference Files

Reference files live in `references/` relative to this `SKILL.md`. They do not auto-load. Claude must Read them when the task reaches the matching generation step.

- `references/plot-dimensions.md` - Default plot chunk dimensions and side-effect chunk options. Read references/plot-dimensions.md before generating R Markdown plot chunks.
- `references/package-reference.md` - Detailed function signatures, parameter lists, defaults, and return values. Read references/package-reference.md before generating calls to less familiar basicBulk, `augere.de`, or `augere.gsea` functions.
- `references/design-patterns.md` - Four contrast patterns, custom design matrix guidance, and complete DGE/GSEA code examples. Read references/design-patterns.md before generating differential expression or GSEA code.
- `references/optional-sections.md` - Optional LogFC Correlations (gg4way + corTriplot), gene plot, gene set score, stats-based GSEA (`augere.gsea::runPrecomputed`), and secondary-factor aesthetics (crossed designs) sections, plus the cache guard for expensive optional sections. Read references/optional-sections.md before generating any optional analysis section.
- `references/code-style.md` - Extended style examples and edge cases. Read references/code-style.md before generating complex tidyverse pipelines, saving code, custom plotting code, tidy-eval code, or when unsure about style.
- `references/output-saving.md` - RData save points, `%T>%` ggsave patterns, batch save idioms, and default file paths per section. Saving is on by default (opt-out only). Read references/output-saving.md before generating output in **template mode**, at **each guide-mode checkpoint**, or whenever the user asks to save outputs, objects, or plots.
- `references/slides.md` - The Quarto RevealJS slide deck: YAML, CSS, the setup chunk (staleness guards + live values), slide structure, figure sizing, rendering, and the visual QC sweep. Read references/slides.md before generating any `.qmd` deck (Section 13).

## Required Workflow and Question Order

Do not generate a template before completing data inspection and design classification.

1. Detect the execution environment and configure R/Pandoc behavior: check the R packages the analysis needs (see Environment Setup).
2. Ask for the input: a SummarizedExperiment (an `.rds` path, or an object in an R session, which is saved to `.rds` first), or a counts table plus a sample sheet to build into one (see Data Inspection). Also ask once for a study title and a one-line study description (they become `PROJECT_TITLE` and `STUDY_DESCRIPTION` in Section 1).
3. Inspect the dataset with R and summarize sample metadata, group values, likely design variables, and possible controls, and check the input contract (see Data Inspection).
4. Confirm the group column, factor level order, and control group. Control must be the first factor level.
5. Determine organism from the SE metadata (if recorded) or ask the user: `"Homo sapiens"` or `"Mus musculus"`, the organisms `basicBulk::getReactome()` supports (mouse Reactome sets are the human MSigDB collection mapped by orthology). Then complete the `symbol` column if it is missing (Data Inspection → Gene annotation).
6. Ask about outliers, metadata subsetting, blocking variables, numeric covariates, repeated-measures `dc.block`, and batch handling. (Outliers are also caught proactively at the QC checkpoint, where you scan the QC metrics yourself rather than relying on the user to name them — see **QC Outlier Scan** in Section 2.) For a Batch column, the batch-handling question offers two options that differ only in the stats parameter — Option A (`dc.block = "Batch"`) or Option B (fixed `block = "Batch"`); both use `adjustBatch = TRUE` in `tidySE`. See Batch Handling.
7. Ask whether the user wants custom logFC or FDR cutoffs. If not specified, omit cutoff arguments.
8. Ask whether to strip control group names from plot labels. Default is to strip: volcano, MA, heatmap, and Overlaps (UpSet) plots use `purrr::set_names(names(.) %>% stringr::str_remove(" vs .*"))` inline for cleaner facet/tab titles. If the user asks to keep full names, omit the `set_names()` call on those plots. `DGEdata` itself always retains full contrast names (needed by the LogFC Correlations section and other functions).
9. Ask which optional sections to include: gene plots, gene set scores, stats-based GSEA (replaces the default GSEA engine), and a slide deck. (The LogFC Correlations section auto-includes when there are 2+ contrasts — do not ask. If the contrasts span ≥ 2 controls each anchoring ≥ 2 contrasts, ask only the split-by-control question; see the LogFC Correlations section below.)
10. Classify the experimental design as Pattern 1, Pattern 2, Pattern 3, Pattern 4, or custom design; then generate or guide through sections in the required order below.

## Analysis CLAUDE.md

Keep a `CLAUDE.md` in the analysis directory that records the decisions a later session can't re-derive
from the code — the group / control choice, outlier calls **and the rationale for each**, batch handling,
custom contrasts, the file layout, and the render order. This is what lets the next session pick the
analysis up without reverse-engineering it.

- **If the directory has no `CLAUDE.md`, create one** as you set the analysis up.
- **If one already exists, ask before updating it** — it may carry decisions you didn't make.
- Refresh it at the natural checkpoints (after the QC / outlier decisions, after the design is fixed,
  after the report and slides render) so it tracks the analysis rather than going stale.
- When the running session has the `claude-md-management:claude-md-improver` skill available, use it
  after editing the file.

This is the analysis's own design record — separate from this skill.

## Environment Setup

On invocation, detect and configure the environment.

### Step 1: Check the Environment

Probe the packages with `requireNamespace()` rather than assuming them:

```bash
Rscript -e 'pkgs <- c("basicBulk", "augere.de", "augere.gsea", "rmarkdown")
for(p in pkgs) cat(p, ":", requireNamespace(p, quietly = TRUE), "\n")
cat("pandoc:", rmarkdown::pandoc_available(), "\n")'
# "Rscript: command not found"  -> no R on PATH at all
```

Once the design and optional sections are settled (Required Workflow step 10), re-run the probe with the extras the chosen analysis uses added to `pkgs`:

| Used by | Extra packages |
|---------|----------------|
| Counts + sample sheet input (symbol mapping) | `org.Hs.eg.db` / `org.Mm.eg.db` |
| Batch correction (`adjustBatch = TRUE`) | `sva` |
| Pattern 3/4 nested facets | `ggh4x` |
| Pathways heatmap tab with two or more pathways (the recipe default), or a palette for more than 20 categories (`pals::glasbey()`) | `pals` |
| Parallel GSEA (many contrasts) | `furrr`, `future`, `parallelly` |
| Stats-based GSEA with `fgsea` | `fgsea` |
| Combined gg4way grid (optional) | `patchwork` |
| Slide deck | `png`, plus the Quarto CLI (not an R package; check with `quarto --version`) |

### Installation

If a package is missing, show:

```r
if(!requireNamespace("BiocManager", quietly = TRUE)){utils::install.packages("BiocManager")}
BiocManager::install("basicBulk", dependencies = TRUE)
```

`dependencies = TRUE` also installs the extras above, which basicBulk suggests, except `fgsea` (`BiocManager::install("fgsea")`); augere.de and augere.gsea (the DGE and GSEA engines) are installed with basicBulk itself.

If `basicBulk::basicUpset()` fails with "Can't merge the `axis.text` theme element" while the other plots render, reinstall ComplexUpset (`utils::install.packages("ComplexUpset")`): a copy built under ggplot2 3 stores stale themes that ggplot2 4 cannot merge.

### Pandoc for Rendering

If `rmarkdown::pandoc_available()` is `FALSE`, silently search for the RStudio-bundled pandoc (macOS):

```bash
find /Applications/RStudio.app -name "pandoc" -type f 2>/dev/null | head -1
```

If found, set `RSTUDIO_PANDOC` before rendering:

```r
Sys.setenv(RSTUDIO_PANDOC = "DETECTED_DIR")
rmarkdown::render("report.Rmd")
```

Do not ask the user about this; detect and use it.

### Rendering

To render the report, run from the analysis directory:

```bash
Rscript -e 'rmarkdown::render("report.Rmd")'
```

Output is `report.html`, next to the `.Rmd` in the working directory. Because the GSEA-over-Reactome render is multi-minute, run it as a background task and poll for completion rather than blocking — the same applies to the slide-deck render in `references/slides.md`.

### Output Saving

Always create the output directory tree in setup (the `purrr::walk(dir.create, ...)` pipe at the top of the Section 1 chunk) unless the user explicitly asks not to: `RData/`, `tables/` and `plots/` with subdirectories. Creating it saves nothing. These paths are all relative to the R working directory, so set the working directory once to a project directory the user chooses before creating the tree, and keep every `save()`/`ggsave()`/`png()` path relative to it. This is the same `setwd(root_dir)` contract the slides setup chunk uses (see `references/slides.md`), so the report and slides write into the same tree. Saving outputs is on by default (opt-out only): in template mode every generated `.Rmd` includes RData saves at each milestone and plot saves at each visualization step (`%T>%` ggsave for ggplot-based plots; a `png`/`pdf` graphics device for heatmaps — see `references/output-saving.md`); in guide mode save after each checkpoint. Read `references/output-saving.md` for RData save points, `%T>%` ggsave patterns, and default file paths before generating output (template mode and each guide-mode checkpoint), or whenever the user asks to save. Only omit saves if the user explicitly opts out. This does not affect `augere.de::runVoom()` or `augere.gsea::runContrast()`, which always use `output.dir = tempfile(), save.results = FALSE` (a fresh directory per call).

## Data Inspection

Run R from the shell to inspect the dataset and infer experimental design. Ask the user for the input: a SummarizedExperiment (an `.rds` path, or an object in an R session), or a counts table plus a sample sheet.

### Load and Inspect

Both input routes end in the same contract, which the inspection snippet below checks:

- unique, non-empty gene IDs (rownames) and sample IDs (colnames);
- a `counts` assay of finite, non-negative numbers;
- a `rowData` column `ID` equal to `rownames(se0)` (a real column — downstream joins use it) and a `symbol` column;
- unversioned Ensembl gene rownames (`ENSG…`, `ENSMUSG…`), the IDs the `basicBulk::getReactome()` gene sets use;
- sample metadata matched and ordered to the count columns by an explicit sample key.

**Branch A — a SummarizedExperiment.** If the SE exists only as an object in an R session, have it saved first (`saveRDS(OBJECT, "SE_PATH_HERE.rds")`): the report renders in a fresh R process, which can only load it from disk.

```bash
Rscript -e '
library(tidySummarizedExperiment)
se0 <- readRDS("SE_PATH_HERE.rds")
cat("=== colData columns ===\n")
print(names(SummarizedExperiment::colData(se0)))
cat("\n=== Sample count ===\n")
print(nrow(SummarizedExperiment::colData(se0)))
cat("\n=== Group values ===\n")
print(if(is.factor(se0$GROUP_COL)) levels(se0$GROUP_COL) else unique(se0$GROUP_COL))
cat("\n=== Numeric columns ===\n")
cd <- as.data.frame(SummarizedExperiment::colData(se0))
print(names(cd)[sapply(cd, is.numeric)])
cat("\n=== Categorical columns (few unique values) ===\n")
cat_cols <- names(cd)[sapply(cd, function(x) is.factor(x) || (is.character(x) && length(unique(x)) <= 20))]
for(col in cat_cols) cat(col, ":", paste(unique(cd[[col]]), collapse=", "), "\n")
cat("\n=== Input contract ===\n")
ids_ok <- function(x) !is.null(x) && !anyDuplicated(x) && all(!is.na(x) & nzchar(x))
has_counts <- "counts" %in% SummarizedExperiment::assayNames(se0)
counts <- if(has_counts) as.matrix(SummarizedExperiment::assay(se0, "counts"))
cat("assays:", SummarizedExperiment::assayNames(se0), "\n")
cat("counts assay, finite and non-negative:", has_counts && is.numeric(counts) && all(is.finite(counts)) && all(counts >= 0), "\n")
cat("unique non-empty IDs | genes:", ids_ok(rownames(se0)), "| samples:", ids_ok(colnames(se0)), "\n")
cat("rowData columns:", names(SummarizedExperiment::rowData(se0)), "\n")
cat("rowData ID == rownames:", identical(SummarizedExperiment::rowData(se0)$ID, rownames(se0)), "\n")
cat("unversioned Ensembl gene IDs:", all(grepl("^ENS[A-Z]*G[0-9]+$", rownames(se0))), "\n")
print(head(rownames(se0)))
'
```

Replace `SE_PATH_HERE` with the user's path and `GROUP_COL` with the group column. Default group column is `"Group"`, but ask if unclear.

Fix only what needs no guessing — add `SummarizedExperiment::rowData(se0)$ID <- rownames(se0)` when `ID` is missing, and add `symbol` once the organism is confirmed (Gene annotation, below) — then re-save the `.rds`. For anything else (no `counts` assay, non-count values, duplicated or empty IDs, versioned or non-Ensembl rownames), stop and ask the user how to proceed.

**Branch B — a counts table plus a sample sheet.** Read the counts with the gene IDs as rownames (`row.names = 1`) and `check.names = FALSE` so the IDs and sample names survive unchanged, match the sheet to the count columns by an explicit key column, and reject duplicated or unmatched keys (adapt the readers to the file type):

```bash
Rscript -e '
library(SummarizedExperiment)
counts <- as.matrix(utils::read.csv("COUNTS_PATH.csv", row.names = 1, check.names = FALSE))
sheet <- utils::read.csv("SHEET_PATH.csv", check.names = FALSE)
key <- as.character(sheet[["KEY_COL"]])
stopifnot("duplicated sample keys" = !anyDuplicated(key),
          "duplicated count columns" = !anyDuplicated(colnames(counts)),
          "sample keys do not match the count columns" = setequal(key, colnames(counts)))
rownames(sheet) <- key
se0 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts),
                                                  colData = sheet[colnames(counts), , drop = FALSE])
SummarizedExperiment::rowData(se0)$ID <- rownames(se0)
saveRDS(se0, "SE_PATH_HERE.rds")
'
```

Then run the Branch A inspection on the saved `.rds`, so both branches are checked against the same contract.

**Gene annotation (after the organism is confirmed, Required Workflow step 5).** When `symbol` is missing, map it from the Ensembl IDs and fill unmapped genes with their ID, then re-save:

```bash
Rscript -e '
library(basicBulk)
organism <- "ORGANISM"
se0 <- readRDS("SE_PATH_HERE.rds")
symbols <- {if(organism == "Homo sapiens"){org.Hs.eg.db::org.Hs.eg.db}else{org.Mm.eg.db::org.Mm.eg.db}} %>%
  AnnotationDbi::mapIds(keys = rownames(se0), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first") %>%
  unname()
SummarizedExperiment::rowData(se0)$symbol <- dplyr::case_when(is.na(symbols) ~ rownames(se0),   # NA symbols fall back to the ID
                                                              .default = symbols)
saveRDS(se0, "SE_PATH_HERE.rds")
'
```

### Control Group Detection

Auto-suggest the control group by matching factor levels case-insensitively against: PBS, Vehicle, Saline, NTC, siNTC, Untreated, Control, WT, Wild-type, Lipofectamine, Mock, Scramble, Sham.

- **Single match:** Suggest to user for confirmation.
- **Multiple partial matches:** For example, `siNTC_1`, `siNTC_2`; flag as potential averaged-control design requiring Pattern 2 or custom contrasts.
- **No match:** Ask the user which level is the control.

### Design Parameter Detection

| Signal | Parameter | Example |
|--------|-----------|---------|
| Categorical column like Sex, Gender | `block` | `block = "GENDER"` |
| Column named Batch | Ask the user: Option A (`dc.block = "Batch"`) or Option B (fixed `block = "Batch"`); both use `adjustBatch = TRUE` in `tidySE` | See Batch Handling |
| Numeric column like Dose, Takedown | `covariates` | `covariates = "Dose"` |
| Same subject in multiple conditions, such as Animal | `dc.block` | `dc.block = "Animal"` |
| Multiple controls or need for averaged contrasts | Custom `design`/`contrasts` | Pattern 2 or custom |

`dc.block` note: when used, `block` should only specify one field.

### Batch Handling

When a Batch column is detected, recode the source variable to a colData column named exactly `Batch` in Section 1 setup (`se0$Batch <- se0$SOURCE_COL %>% forcats::as_factor()`) — both options below depend on that name, and `adjustBatch = TRUE` errors if the column is missing. Then **present the user with two coherent options and ask which to use** — do not assume one. Neither is universally right; it is a per-study call (see the trade-offs below).

**Normalization and visualization are the same for both options.** `adjustBatch` only affects *visualization* — it builds the `"ComBat"` assay, which the statistics never touch (`augere.de::runVoom()`/`augere.gsea::runContrast()` always refit on the raw `"counts"` assay, passed as `assay = "counts"`). So it is orthogonal to the `dc.block`-vs-`block` choice, and the batch-removed PCA is an equally useful diagnostic whichever option you pick. In **both** options:

- **Normalization:** `basicBulk::tidySE(adjustBatch = TRUE)` — and do **not** pass `block = "Batch"` to `tidySE`. `tidySE` strips `"Batch"` from the ComBat design itself, so a `block = "Batch"` here would not change any stored assay, the ComBat-adjusted one included, nor the stats result. (A *non-batch* `block`/`covariate` such as Sex or Dose still goes in `tidySE` normally.)
- **Visualization:** `plotAssay <- "ComBat"` (passed as `assay = plotAssay`); show both the default and batch-adjusted PCAs side by side — the batch-adjusted panel shows the data with batch *removed*, so the audience can see whether batch was driving the separation. (For Option B, ComBat-removal of batch *means* is a close visual analog to the fixed block; for Option A it is a looser one, since `dc.block` models a correlation rather than a mean shift — fine as a diagnostic either way.) The plain `"logCPM"` panel is the right *uncorrected* baseline precisely because it is **not** design-residualized — `tidySE`'s `"logCPM"` assay is voom's `E` (normalized log2-CPM), which ignores the design — so the comparison reads as "raw batch structure (logCPM) vs batch removed (ComBat)," not two differently-modeled views. That is also why batch is never baked into the `"logCPM"` assay (and why `block = "Batch"` in `tidySE` would not produce a batch-corrected visual — see Normalization above).
- **ComBat fallback:** `sva::ComBat_seq()` needs batch reasonably crossed with Group. If it errors or the batch is near-confounded with Group, set `adjustBatch = FALSE` and visualize plain `"logCPM"` (`plotAssay <- "logCPM"`, single PCA) instead — the statistics are unaffected either way.

**The options differ in one place only: the statistics parameter.**

**Option A — `dc.block`.** `augere.de::runVoom()` and `augere.gsea::runContrast()` get `dc.block = "Batch"` (no `block`). Batch enters the model as a consensus intra-block correlation passed to `lmFit` — not a design column.

**Option B — fixed block.** `runVoom()` and `runContrast()` get `block = "Batch"`. Batch is a fixed effect in the design. (For a **custom design**, `block =` is bypassed — add `+ Batch` to the model matrix instead; see `references/design-patterns.md`.)

**Trade-off to state when asking:** Option A's `duplicateCorrelation` treats batch as a random-effect-style correlation; Option B's fixed-effect block is the more conservative choice when batch is orthogonal to the design (`duplicateCorrelation` can be anticonservative there). Pass the chosen parameter consistently to both `runVoom()` and `runContrast()`; both always use the raw `"counts"` assay (`assay = "counts"`), never ComBat-adjusted values.

## Experimental Design Classification

**Pattern 1 is the default.** If the data suggests a different pattern might apply, confirm with the user before proceeding — do not auto-select Pattern 2, 3, or 4.

- **Pattern 1: Single model, simple comparisons.** Default. Single control group, straightforward pairwise comparisons, standard blocking.
- **Pattern 2: Nested split via `designTibble`.** Groups encode both a biological context and treatment in their names (e.g., `NeuronDrug`, `AstroPBS`). Each comparison is an explicit pair in the `designTibble`.
- **Pattern 3: Single-variable subset via `subset.factor`.** Same comparisons run separately for each level of one variable (e.g., Drug vs PBS subset by CellType). Group column stays clean.
- **Pattern 4: Crossed multi-factor split.** Contrasts are crossed across two variables. For example, Drug vs Vehicle compared within each combination of brain Region and animal Cohort. Uses `tidyr::crossing()` to iterate over all combinations.
- **Custom design matrix.** Use when you must **average ≥2 controls** into one baseline (a single control instead routes to Pattern 3 for one factor, or Pattern 4 for two), for interaction terms, or for complex contrasts. The `design-patterns.md` "Custom Design Matrix" section has a routing table and single-factor vs crossed (e.g. multiple controls × dose) recipes.

When Pattern 2, 3, or 4 is identified as a possibility, ask the user to confirm the pattern and the variable(s) to split/subset by. When Pattern 2 is selected, warn: "Pattern 2 support is experimental — please review the generated code carefully before running." When Pattern 3 is selected, warn: "Pattern 3 support is experimental — please review the generated code carefully before running." When Pattern 4 is selected, warn: "Pattern 4 support is highly experimental — the generated visualization code may need manual adjustment."

Read `references/design-patterns.md` before generating `comparisonList`, `voomList`, `DGEdata`, `sigTibble`, `gene.sets`, `GSEAlist`, or custom `design`/`contrasts` code.

## Report Prose Style

The rendered HTML report is a reference document, not an essay — keep to this house style. Concretely:

- **Section intros are 1–2 sentences, descriptive and methodological** — what the analysis does and what
  a plot shows ("The loadings represent the most positively (red) and negatively (blue) correlated genes
  with each PC."), with the tool/package cited and linked. Not a paragraph, not a subtitle, not a blockquote.
- **No interpretation in the report prose.** Don't rank findings, call a result strong / clean /
  concerning, draw biological conclusions, or give a safety verdict in the report text. The tables and
  figures carry the results; the reader draws the conclusions. Interpretation is the slide deck's job,
  worked out with the user — see Section 13 and the report→slides ordering there.
- **Inline `r` only for quantities** — counts, sample numbers, gene-set sizes — never for prose
  conclusions, and kept parenthetical.
- Plain text — no runs of **bold** / *italic* emphasis, no pep.
- **Cite the external engine that performs each step, not just the basicBulk wrapper.** A basicBulk function usually delegates to a well-known engine; name (and link) that engine in the methods sentence — the wrapper is the *how*, the engine is what belongs in the methods description:

  | Step | Cite (engine) |
  |------|---------------|
  | DGE | limma-voom (via augere.de) |
  | Batch adjustment | ComBat-seq / sva |
  | PCA | PCAtools |
  | Heatmaps | ComplexHeatmap |
  | Pathways | camera + REACTOME (MSigDB C2:CP:REACTOME via msigdbr; via augere.gsea) |

  Caution: `basicBulk::basicVolcano(style = "enhanced")` is basicBulk's **own** ggplot2 volcano styling — attribute the volcano step to basicBulk/ggplot2, **not** the EnhancedVolcano package.

## Workflow Sections

Generate code for each section in this order. In guide mode, show each chunk with explanation. In template mode, assemble into a complete `.Rmd`.

### Cutoff Parameters

When the user specifies custom logFC or FDR cutoffs, define `logFC_cutoff` and `FDR_cutoff` once after normalization and pass them consistently to all of these functions, which default to `logFCcutoff = 1` and `FDRcutoff = 0.05`:

- `basicBulk::makeSigTibble(logFCcutoff = logFC_cutoff, FDRcutoff = FDR_cutoff)`
- `basicBulk::basicVolcano(logFCcutoff = logFC_cutoff, FDRcutoff = FDR_cutoff)`
- `basicBulk::basicMA(logFCcutoff = logFC_cutoff, FDRcutoff = FDR_cutoff)`
- `basicBulk::basicUpset(logFCcutoff = logFC_cutoff, FDRcutoff = FDR_cutoff)`
- `basicBulk::genePlot(logFCcutoff = logFC_cutoff)`
- `gg4way::gg4way(logFCcutoff = logFC_cutoff, FDRcutoff = FDR_cutoff)`

Do not change `augere.de::runVoom(lfc.threshold = ...)`; that controls the limma treat test, which is a different statistical test. If the user does not specify custom cutoffs, do not define cutoff variables and omit these parameters entirely.

### 1. Setup

```r
library(basicBulk)
library(tidySummarizedExperiment)

c("RData", "tables", "plots/QC", "plots/DGE/ContrastVolcanos", "plots/DGE/4way",
  "plots/heatmaps", "plots/GSEA", "plots/scores") %>%
  purrr::walk(dir.create, recursive = TRUE, showWarnings = FALSE)

PROJECT_TITLE <- "TITLE"               # asked once (Required Workflow step 2)
STUDY_DESCRIPTION <- "DESCRIPTION"

se0 <- "SE_PATH_HERE.rds" %>%
  readRDS() %>%
  # Optional: subset by metadata
  dplyr::filter(FILTER_CONDITION)

se0$Group <- se0$GROUP_SOURCE_COL %>%
  factor(levels = c("CONTROL", "TREATMENT1", "TREATMENT2"))

group <- "Group"
groups <- se0$Group %>%
  levels()

organism <- "ORGANISM"   # "Homo sapiens" or "Mus musculus"

outliers <- c("SAMPLE_ID")  # NULL if none
if(!is.null(outliers)){
  se0 <- se0[, !colnames(se0) %in% outliers]
}
```

`se0` always loads from the `.rds` checked in Data Inspection (both input branches end there). `PROJECT_TITLE` and `STUDY_DESCRIPTION` feed the report YAML (see Template Mode), `RData/settings.RData`, and the deck's Overview slide.

### 2. Normalization

R Markdown header: `# Normalization` — a **flat** header, no `{.tabset}`: one sentence of methodological prose (what `tidySE()` adds, citing voom and, with batch correction, ComBat-seq), then the normalization chunk. `tidySE()` draws no plot, and a one-tab tabset would render as a lone orphan pill.

```r
se <- se0 %>%
  basicBulk::tidySE(group = "Group",
                    block = BLOCK,
                    covariates = COVARIATES,
                    adjustBatch = ADJUST_BATCH)

plotAssay <- "logCPM"

gene.sets <- basicBulk::getReactome(organism)

gsea_methods <- "camera"
```

Only include `block`, `covariates`, and `adjustBatch` when they apply. If there is no blocking, covariates, or batch correction, generate just `basicBulk::tidySE(group = "Group")`. For a Batch column, set `adjustBatch = TRUE` and do **not** put `"Batch"` in `block` — this is the same for Batch Handling Option A and Option B (the options differ only in the stats parameter, not here); it builds the `"ComBat"` assay and strips `"Batch"` from the ComBat design itself. (Fallback: if `sva::ComBat_seq()` errors or Batch is near-confounded with Group, use `adjustBatch = FALSE`.) A non-batch `block`/`covariates` (e.g. Sex, Dose) applies regardless.
Generate `plotAssay <- "ComBat"` whenever `adjustBatch = TRUE` (both batch options), and `plotAssay <- "logCPM"` only with no batch correction (or the ComBat fallback). Define `plotAssay`, `gene.sets`, and `gsea_methods` once after normalization; these shared variables are used by all downstream DGE, GSEA, and visualization code, all on the one `se`. `gsea_methods <- "camera"` (see §10 Pathways) is the single source the GSEA loop's `methods =`, its `tidyGSEAlist(GSEAmethod = ...)`, and the GSEA cache key all reference, so they can't drift. When custom cutoffs are user-specified, define `logFC_cutoff` and `FDR_cutoff` once after this setup and pass them to cutoff-aware plotting/table functions; otherwise do not define those objects or pass cutoff arguments. The setup chunk creates the output folders by default — see `references/output-saving.md` for save patterns. This does not change `augere.de::runVoom()`/`augere.gsea::runContrast()` arguments, which always use `assay = "counts", output.dir = tempfile(), save.results = FALSE`.

### QC Outlier Scan

Don't just render the PCA and move on — **read it**. Once the PCA (Section 3) and the library sizes are
available, proactively scan them for outlier samples and surface them to the user. A sample that is an
outlier on a *quality* axis (rather than a biological one) can distort normalization and every
downstream contrast, and the author expects these caught and raised — not silently passed, and not
silently dropped. Cover the metrics the report actually has:

- **Library size** (`se$lib.size`, added by `tidySE()`): read depth / library size.
- **PCA**: a sample sitting far from its group or with another group (a possible sample swap), or driving a PC
  largely on its own.

For each flagged sample, give the user a **quantitative rationale** — which metric, and how far from the
cohort it sits (a robust distance such as a median/MAD or z, e.g. "~6 MAD below the median library size,
~3× smaller than the next sample"). The read is "how unusual is this sample on this axis," not a hard
cutoff. Then **stop and ask** — removing a sample or correcting its label is the user's call, never
yours: present the evidence with a recommendation (remove, or keep and note; for a likely sample swap,
also correct the label) and act only on their confirmation, the same "ask, don't assume" stance as the
design-pattern and batch decisions. A confirmed removal goes in the `outliers <-` vector in Section 1.
If the user keeps a flagged sample, note it factually in one line in the report (e.g. that a low-depth
sample was retained) — descriptively, not interpreted. A sample can be an outlier on one axis and fine
on every other, so name the axis, so the user can weigh it.

### 3. PCA

R Markdown header: `# PCA`.

```r
se %>%
  basicBulk::PCAgrid(assay = plotAssay,
                     colby = "Group",
                     selectLab = NA,
                     alpha = 0.8)
```

For a crossed design (e.g. molecule × dose), colour by the primary factor and shape by the secondary one instead of putting the full crossed `Group` on `colby` — `PCAgrid(assay = plotAssay, colby = "PRIMARY_FACTOR", shape = "SECONDARY_FACTOR")`. See `references/optional-sections.md` (Secondary-Factor Aesthetics).

When `adjustBatch = TRUE` (the default for **both** batch options, Option A and Option B), show both default and batch-adjusted PCAs as tabset pills, colored by both Group and Batch (this needs the `"ComBat"` assay). Only if the ComBat fallback applied — `adjustBatch = FALSE`, no ComBat assay — emit just the single `plotAssay` PCA above instead:

```r
PCAtibble <- tidyr::crossing(
  "colby" = c("Batch", "Group") %>%
    forcats::as_factor() %>%
    forcats::fct_rev(),
  "assay" = c("logCPM", "ComBat") %>%
    forcats::as_factor() %>%
    forcats::fct_rev()
) %>%
  dplyr::mutate(tidyName = glue::glue("{colby} PCA ({tidyAssay})",
                                      tidyAssay = dplyr::recode_values(assay,
                                                                       "logCPM" ~ "Default",
                                                                       "ComBat" ~ "Batch-adjusted")))

PCAtibble %>%
  purrr::pmap(~ with(list(...),
                     se %>%
                       basicBulk::PCAgrid(assay = assay,
                                          colby = colby,
                                          selectLab = NA,
                                          alpha = 0.8))) %>%
  purrr::set_names(PCAtibble$tidyName) %>%
  purrr::iwalk(basicBulk::printPlots)
```

The batch-adjusted PCA section uses `purrr::iwalk(basicBulk::printPlots)`, so the parent header must be `{.tabset .tabset-pills}` and the chunk must use `results="asis"`.

### 4. Differential Gene Expression

R Markdown header: `# Differential Gene Expression` — a **flat** header, no `{.tabset}` (a one-tab tabset renders as a lone orphan pill; the DEG-count bar and contrast counts live under Volcano Plots, not here). A plain `## Significant Genes` **subheading** is fine here — it is not a `{.tabset}`, so it does not orphan-pill.

Structure the section as **intro prose -> a compute chunk -> a `## Significant Genes` display chunk**, keeping compute SEPARATE from display:

1. A 1-2 sentence methodological intro under the header (the model, and what each contrast compares).
2. A compute chunk (use the detected pattern from `references/design-patterns.md` for `augere.de::runVoom()`) that builds and **saves** `DGEdata` and `sigTibble` — do NOT print the table from here. A buried `%T>% basicBulk::printSigTibble()` inside the compute chunk gives no subheading and no caption, and renders nothing (the print helpers return the table widget rather than printing it); that is the shape to avoid.
3. A `## Significant Genes` subsection: a one-line caption naming the thresholds, then a dedicated display chunk that calls `basicBulk::printSigTibble()` (its searchable DEG table is the section's single output).

```r
# compute chunk (build + save; no printing) — runVoom per the detected pattern (Pattern 1 shown)
voomList <- se %>%
  augere.de::runVoom(groups = "Group",
                     comparisons = comparisonList,
                     assay = "counts",
                     output.dir = tempfile(),
                     save.results = FALSE) %>%
  magrittr::extract2("results")
DGEdata   <- se %>% basicBulk::makeDGEdata(voomList)
sigTibble <- DGEdata %>% basicBulk::makeSigTibble()

# display chunk, under a `## Significant Genes` header preceded by a caption line
# e.g. "Genes passing FDR < 0.05 and |LogFC| > 1 (makeSigTibble defaults; state your actual cutoffs)."
sigTibble %>% basicBulk::printSigTibble()
```

### 5. Volcano Plots

Pass `assay = plotAssay` to every expression plot below and in the later sections, including heatmaps, gene plots, and gene set scores (the basicBulk plotters default to `assay = "logCPM"`, which leaves out the ComBat adjustment when `adjustBatch = TRUE`). Statistics always use the raw `"counts"` assay (`assay = "counts"`).

**Pattern 3/4 note:** The visualization code below is for Pattern 1/2. For Pattern 3, read `references/design-patterns.md` "Pattern 3 Downstream Visualization" — `tidyContrasts()` + `ggh4x::facet_nested()` for the volcano and MA grids (the volcano is `basicVolcano()` plus the count-layer step, which adds the `tidyContrasts()` columns to its DEG-count label layer before re-faceting), `ggplot2::facet_wrap()` for its other plots, and pre-filtered heatmaps. For Pattern 4, read "Pattern 4 Downstream Visualization" — `basicVolcano()` plus the same count-layer step, `tidyContrasts()` + `ggh4x::facet_nested()`, and pre-filtered heatmaps.

R Markdown header: `# Volcano Plots {.tabset .tabset-pills}`. DEG Magnitude, Volcano Grid, and MA Grid each get a `##` tab.

```r
## DEG Magnitude
sigTibble %>%
  basicBulk::makeSigCountTibble() %>%
  basicBulk::DEGbar()

## Volcano Grid
DGEdata %>%
  purrr::set_names(names(.) %>% stringr::str_remove(" vs .*")) %>%
  basicBulk::basicVolcano(style = "enhanced")   # basicBulk's OWN ggplot2 volcano styling — NOT the EnhancedVolcano package

## MA Grid
DGEdata %>%
  purrr::set_names(names(.) %>% stringr::str_remove(" vs .*")) %>%
  basicBulk::basicMA()
```

**After the report renders, look over the MA grid** (`plots/DGE/MA_grid.png`) for any contrast whose DEGs are
mostly low abundance (bunched at the low end of the Average logCPM axis), and flag it to the user with the evidence.

### 6. Heatmaps

R Markdown header: `# Heatmaps {.tabset .tabset-pills}`. Split into two chunks: first builds the list with `include=FALSE`, second prints as tabset pills with `results="asis"`. The build chunk uses `include=FALSE`, not `fig.show="hide"`: it wraps the `png()`/`pdf()` device saves (see `references/output-saving.md`), which print "Saving … image" device messages — `fig.show="hide"` suppresses the figures but not those messages, so they leak into the report. `include=FALSE` still runs the chunk, so `allHeatmap`/`heatList` remain available to the print chunk.

```r
# Chunk 1: include=FALSE
allHeatmap <- se %>%
  basicBulk::basicHeatmap(sigTibble = sigTibble,
                          assay = plotAssay)

heatList <- names(DGEdata) %>%
  purrr::set_names(stringr::str_remove(., " vs .*")) %>%
  purrr::map(~ se %>%
               basicBulk::basicHeatmap(se = .,
                                       sigTibble = sigTibble,
                                       assay = plotAssay,
                                       contrast = .x))
```

```r
# Chunk 2: results="asis"
c(list("All DEGs" = allHeatmap),
  heatList) %>%
  purrr::iwalk(basicBulk::printHeatmaps)
```

**Customizing the column annotations.** `basicHeatmap()` annotates the columns by `Group` by default. To
add another annotation track (e.g. Sex, Dose, Batch) with your own colours, pass a `top_annotation` — a
`ComplexHeatmap::HeatmapAnnotation` with a `col` list that keys each track to its factor levels so the
colours map deterministically (`ggsci::pal_d3()` over `levels(Group)`; a named hex vector for the
secondary track). Build it in the `include=FALSE` chunk:

```r
col <- list("Group" = ggsci::pal_d3()(n = nlevels(se$Group)) %>%
              purrr::set_names(levels(se$Group)),
            "Sex"   = c("F" = "#8624F5", "M" = "#1FC3AA"))   # named to the track's levels

top_annotation <- ComplexHeatmap::HeatmapAnnotation("Group" = se$Group,
                                                    "Sex"   = se$Sex,
                                                    col = col,
                                                    annotation_legend_param = list(
                                                      labels_gp = grid::gpar(fontsize = 9)),
                                                    border = TRUE)

allHeatmap <- se %>%
  basicBulk::basicHeatmap(sigTibble = sigTibble,
                          assay = plotAssay,
                          top_annotation = top_annotation)
```

`vignette("basicBulk")` ("Column annotations") is the authoritative example; the Pattern 3/4 heatmaps in
`references/design-patterns.md` use the same construction to annotate a subset/split column. `save(col,
top_annotation, ...)` per `references/output-saving.md` since the objects are built by hand, and pass the
same `top_annotation` to the per-contrast `heatList` map if you want the tracks on every tab.

**After the report renders, look over the heatmaps** (the All DEGs heatmap across all samples, and each contrast's
heatmap across its samples) for a sample that doesn't match the rest of its group, or looks like or clusters with
another group: a possible outlier or sample swap. To see which sample it is, redraw that heatmap to a temporary file
with `show_column_names = TRUE`. Flag it to the user with the evidence, as in the QC Outlier Scan. An outlier on the
heatmaps alone is worth flagging; one that also stands out on the PCA is a strong outlier.

### 7. Overlaps

R Markdown header: `# Overlaps {.tabset .tabset-pills}`. Open with a one-line intro ("This section examines the overlaps of the DGE profiles across contrasts."). Three `##` tabs — UpSet, Table, Overlap Expression — so the reader sees the intersections, the genes in them, and what those shared genes look like across groups. (The tabs are plain `##` with a single plot/table each, so they do not call `printPlots`/`printHeatmaps` and do not need `results="asis"`.)

```r
## UpSet
DGEdata %>%
  purrr::set_names(names(.) %>% stringr::str_remove(" vs .*")) %>%
  basicBulk::basicUpset(sort_sets = "descending")
```

```r
## Table
# get_dupes already orders by descending dupe_count (the "Hits" column), so the top overlaps
# sit first — no extra arrange needed. The contrast mutate is cosmetic (cleaner labels); it does
# not drive the dedup, which keys on (ID, direction). printSigTibble(rearrange = FALSE) preserves
# the get_dupes order. If the user opted to keep full contrast names (Step 8), drop the mutate.
sigTibble %>%
  dplyr::mutate(contrast = stringr::str_remove(contrast, " vs .*")) %>%
  janitor::get_dupes(ID, direction) %>%
  basicBulk::printSigTibble(rearrange = FALSE)
```

```r
## Overlap Expression
# A few top shared genes — the genes that recur across contrasts — shown as expression dot plots.
topOverlap <- sigTibble %>%
  dplyr::mutate(contrast = stringr::str_remove(contrast, " vs .*")) %>%
  janitor::get_dupes(ID, direction) %>%
  dplyr::arrange(FDR) %>%
  dplyr::distinct(symbol) %>%
  dplyr::slice(seq_len(min(N_TOP_OVERLAP, dplyr::n()))) %>%   # data-driven cap (e.g. 20); never over-slice
  dplyr::pull(symbol)

if (length(topOverlap) > 0) {                                # a study with no overlaps would otherwise fail genePlot's facet_wrap(~ symbol)
  se %>%
    basicBulk::genePlot(genes = topOverlap, assay = plotAssay, plotType = "dot", ncol = 3)
}
```

For a crossed design (e.g. molecule × dose), build the Overlap Expression plot per-sample with the primary factor on x and the secondary factor by shape rather than the full crossed `Group` — see `references/optional-sections.md` (Secondary-Factor Aesthetics). Save each tab's output with `%T>% ggplot2::ggsave()` per `references/output-saving.md`.

### 8. LogFC Correlations

R Markdown header: `# LogFC Correlations` — add `{.tabset .tabset-pills}` only at ≥ 3 contrasts (when the section has `##` children); at exactly 2 contrasts the section is a single plot with no tabs, so drop the pills (see `references/optional-sections.md` routing). Auto-include when there are two or more contrasts; only omit if the user explicitly asks not to. This section pairs a per-contrast-pair `gg4way` view with an overview correlation triangle (`basicBulk::corTriplot()`, fed by `basicBulk::prepareCor(DGEdata = DGEdata, cd = SummarizedExperiment::colData(se))` with the pair tibble piped in as its first argument; `cd` is required) that appears once there are ≥ 3 contrasts. At exactly 2 contrasts it collapses to a single `gg4way` (no triangle). Keep the overview triangle **compact** — it is small by default and should not be enlarged at low/moderate contrast counts; resize it up only when a high contrast count makes the cells unreadable (see `references/plot-dimensions.md` for dimensions, and `references/optional-sections.md` for the split-by-control alternative at high N). When the contrasts span two or more controls that each anchor ≥ 2 contrasts (e.g. a dose-matched design), ask the user whether to split the overview into one triangle per control — and surface that read-with-care caveat in the narrative (and the slide take-home, if a deck is built). Read `references/optional-sections.md` before generating this section.

### 9. Gene Plots

Only include when the user provides genes of interest or asks for expression histograms. Read `references/optional-sections.md` before generating this section.

### 10. Pathways

R Markdown header: `# Pathways {.tabset .tabset-pills}`. The chunk for `augere.gsea::runContrast()` uses `fig.show="hide"`. Top GSEA Terms, Significant Pathways, and Heatmap get `##` headers. `basicBulk::GSEAplotList()` plus `purrr::iwalk(basicBulk::printPlots)` generates `##` tabs directly (the GSEA barplots), so do not add an extra explicit header before that output. `topGSEAplot` dimensions are N-dependent — read `references/plot-dimensions.md` before assigning chunk dimensions.

`gene.sets` is defined after normalization (Section 2, the Shared Template Setup in `references/design-patterns.md`) alongside `plotAssay` and `gsea_methods`, so it is already available here.

Use the detected pattern from `references/design-patterns.md` to loop `augere.gsea::runContrast()`. `runContrast()` has no seed argument: call `set.seed(100)` once before the loop.

**Pass `methods = gsea_methods` (`"camera"`) in every `runContrast()` call** (defined once in setup; the design-pattern loops show this). The package default `methods = c("fry", "camera")` also runs the additional deterministic `fry` test (limma's analytic approximation to the `mroast` rotation gene-set test — not a random-rotation test), but `basicBulk::tidyGSEAlist()` defaults to `GSEAmethod = "camera"` and extracts only that, so `fry` is discarded in the standard flow and compounds badly over many contrasts. Opt into `fry` only when you genuinely need its self-contained test — then set `gsea_methods` and pass a scalar `GSEAmethod` to the tidy step.

**Many contrasts (dozens+): parallelize the GSEA loop with furrr.** `camera` uses no RNG, so parallelizing is safe: serial `purrr::imap()` and multicore `furrr::future_imap()` should give identical results — confirm with `identical(serial, par)`, and re-validate exact identity when you change platform, package versions, or method. Use the serial idiom by default; switch to the furrr form only when you have dozens of contrasts (fork plan + `furrr_options(seed = TRUE)` + restore-the-plan-on-exit). It **must degrade to serial** when parallelism isn't available or worthwhile — `furrr`, `future` and `parallelly` are not hard basicBulk dependencies (Suggests only; an install may lack them), fork-based `future::multicore` is unsupported on Windows, and with only 1–2 cores `workers` collapses to 1 (furrr overhead, no gain). So gate on an executable `parallel_ok` and keep the serial `purrr::imap()` loop in the `else`. Swap `purrr::imap()` for `furrr::future_imap()` (or `purrr::pmap()` → `furrr::future_pmap()` for the Pattern 2 loop), keeping every `runContrast()` argument (`methods = gsea_methods`, `assay = "counts"`, plus any `block`/`covariates`/`dc.block` your pattern uses) unchanged. Each call keeps its own `output.dir = tempfile()`: `runContrast()` writes its report into `output.dir`, so concurrent tasks sharing one directory would overwrite each other.

```r
parallel_ok <- requireNamespace("furrr", quietly = TRUE) &&
               requireNamespace("future", quietly = TRUE) &&
               requireNamespace("parallelly", quietly = TRUE) &&
               parallelly::supportsMulticore() &&             # fork multicore (FALSE on Windows/RStudio)
               parallelly::availableCores(omit = 1) > 1L      # leave a core for the parent; 1-2 cores -> serial

set.seed(100)   # runContrast() has no seed argument
if (parallel_ok) {
  oldplan <- future::plan(future::multicore, workers = parallelly::availableCores(omit = 1))
  GSEAlist <- tryCatch(
    comparisonList %>%
      purrr::set_names(purrr::map_chr(., magrittr::extract, 1)) %>%
      furrr::future_imap(~ augere.gsea::runContrast(se,
                                                    sets = gene.sets,
                                                    groups = "Group",
                                                    comparison = .x,
                                                    methods = gsea_methods,
                                                    assay = "counts",
                                                    output.dir = tempfile(),   # a unique directory per task
                                                    save.results = FALSE,
                                                    suppress.plots = TRUE),    # no shared graphics device across forked workers
                         .options = furrr::furrr_options(seed = TRUE),   # defensive; camera/fry use no RNG
                         .progress = TRUE) %>%
      basicBulk::tidyGSEAlist(GSEAmethod = gsea_methods),
    finally = future::plan(oldplan)   # restore prior plan on success OR error — no leaked multicore
  )
} else {
  GSEAlist <- comparisonList %>%      # serial fallback (identical result, just slower)
    purrr::set_names(purrr::map_chr(., magrittr::extract, 1)) %>%
    purrr::imap(~ augere.gsea::runContrast(se,
                                           sets = gene.sets,
                                           groups = "Group",
                                           comparison = .x,
                                           methods = gsea_methods,
                                           assay = "counts",
                                           output.dir = tempfile(),
                                           save.results = FALSE),
                .progress = TRUE) %>%
    basicBulk::tidyGSEAlist(GSEAmethod = gsea_methods)
}
```

`fry` is *also* deterministic — limma implements it as a fast analytic approximation to the `mroast` rotation test, with no rotations/RNG of its own — so `methods = c("fry", "camera")` would parallelize just as safely; the identity check above is simply stated for camera-only runs. `mroast` and `romer`, by contrast, use random rotations: `furrr_options(seed = TRUE)` makes future executions reproducible across `future` plans, but does **not** by itself guarantee equality with a serial `purrr` run, and `augere.gsea` draws a fresh seed for each `mroast`/`romer` run from the session's random-number stream — so validate that combination separately before parallelizing `mroast`/`romer`.

**Cache the GSEA compute behind the content-hash freshness guard** so a full render reuses it safely and never serves a stale Pathways section: wrap the `runContrast()` loop in the `basicBulk::cache_key()` / `load_if_fresh()` / `save_cache()` pattern (see `references/output-saving.md` → After GSEA), keyed on the sample set (`colData(se)`), the genes (`rownames(se)`), the realized count values (`as.matrix(SummarizedExperiment::assay(se, "counts"))` — an input file can change without changing sample or gene names), your pattern's contrast spec (Pattern 1 `comparisonList`; Pattern 2 `designTibble`; Patterns 3/4 `comparisonList` + the subset/split column names; custom `design` + `contrastList`), the model-argument selections (`block`/`covariates`/`dc.block`/`subset.*` values — a colData hash alone misses a `dc.block` column switch), the `methods` selection (key on the *same* object the loop passes to `methods =` — the shared `gsea_methods <- "camera"` — so a `camera` ↔ `fry` switch can't drift between the call and the cache key), the tidy step's `GSEAmethod`, and `gene.sets`. Key on the INPUTS only — **no automatic package/engine-version hash**: reuse across an upgrade assumes unchanged GSEA behaviour, so review the relevant engine / dependency / result-tidying changes and force a recompute (`run_gsea <- TRUE`) when results may change or the impact is uncertain. Manual invalidation is the maintenance responsibility that replaces an automatic version token.

**GSEA engine — default vs opt-in.** The default engine is `augere.gsea::runContrast()` (expression/contrast-based, looped per the detected design pattern). If the user requests **stats-based / pre-ranked GSEA**, replace that compute step with `augere.gsea::runPrecomputed()` — a single design-agnostic map over the already-built `DGEdata`, `methods = "cameraPR"` by default — and read `references/optional-sections.md` ("Stats-based GSEA"). Its methods report `Direction` in lowercase, so the recipe there title-cases it after the tidy step; with that, it is a drop-in replacement that produces the same tidied `GSEAlist` (and its own report prose), so everything below (`topGSEAplot`/`GSEAplotList`/`printGSEAtibble`) is unchanged.

```r
## Top GSEA Terms
GSEAlist %>%
  purrr::map(~ dplyr::mutate(.x, Term = stringr::str_trunc(Term, 44, side = "right"))) %>%
  basicBulk::topGSEAplot(n_terms = 5)

# No ## header; printPlots generates the ## tabs
# Chunk options: results="asis"
GSEAlist %>%
  basicBulk::GSEAplotList() %>%
  purrr::imap(~ .x + ggplot2::labs(title = .y)) %>%
  purrr::iwalk(basicBulk::printPlots)

## Significant Pathways
GSEAlist %>%
  basicBulk::printGSEAtibble()
```

**Heatmap tab — pathway expression.** Once the GSEA results are in (the guide-mode GSEA checkpoint; in template mode, after the first render), ask the user which pathways to show, suggesting the top significant terms from the Significant Pathways table, then add the tab. `basicBulk::makePathwayTibble()` returns one row per (pathway, gene), so `dplyr::distinct(ID, .keep_all = TRUE)` keeps each gene once, under its first pathway, while `.keep_all` keeps the `Pathway` column that `basicHeatmap(pathways = TRUE)` splits the rows by. Build it inside the `png()` device in an `include=FALSE` chunk, as in Section 6, and print it under the `## Heatmap` header (dimensions from `references/plot-dimensions.md`, save path from `references/output-saving.md`):

```r
# Build chunk: include=FALSE
pathways <- c("TERM 1", "TERM 2") %>%   # Term values from the Significant Pathways table
  paste0("REACTOME_", .) %>%
  stringr::str_replace_all(" ", "_")

png("plots/heatmaps/pathway heatmap.png", res = 300, units = "in", width = WIDTH, height = HEIGHT)
pathwayHeatmap <- pathways %>%
  basicBulk::makePathwayTibble(gene.sets = gene.sets) %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
  basicBulk::basicHeatmap(se = se,
                          sigTibble = .,
                          pathways = TRUE,
                          assay = plotAssay,
                          show_row_names = FALSE)
dev.off()
```

```r
## Heatmap
pathwayHeatmap
```

### 11. Gene Set Scores

Optional. Ask the user if they want scores for specific pathways (or for the ones already chosen for the Pathways heatmap). If the top enrichments are dominated by a redundant family (e.g. many near-duplicate ribosome/translation pathways), the greedy-Jaccard `selectDiversePathways()` selection is a good default for the first pass. Read `references/optional-sections.md` before generating this section.

**After the report renders, review the selected pathways against the top GSEA terms** (the `topGSEAplot` dot plot from the Top GSEA Terms tab above). A diversity heuristic optimises for non-redundancy, not biological salience, so it can pass over a distinct, important theme — e.g. a stress/senescence signal in a perturbation study — in favour of a less-relevant but more "different" pathway. If the selection doesn't reflect the dominant biology, offer the user an explicit named-pathway list (see `references/optional-sections.md`) and swap only on their go. This is a post-hoc call: it needs the rendered results in hand.

### 12. Session Info

```r
utils::sessionInfo()
```

In template mode, wrap in the collapsible CSS button pattern below.

### 13. Slide Deck (optional)

**When to offer:** After completing a template analysis (or when the user explicitly asks for
"slides"), proactively offer a slide deck. "Slides", "deck", or "presentation" are all sufficient
triggers — do not require the user to say "RevealJS".

**Order of work — report first, then a first-draft deck, then interpretation with the user.** Generate
the rendered HTML report (descriptive prose, per Report Prose Style) and then a **first-draft** deck whose
finding slides carry factual `.takehome` lines read straight from the data. The biological story — the
ranking, what matters, the take-home emphasis — is layered into the **slides** in a second pass *with the
user*: present a short summary of the results and your suggested take-homes in chat, and revise the deck
on their confirmation. This keeps interpretation out of the report (where the house style is descriptive)
and puts it where it belongs, the deck, while keeping the user in the loop on every interpretive claim.
The no-safety-verdict rule holds throughout (see `references/slides.md`).

Read `references/slides.md` before generating any `.qmd`. Key rules in brief:

- Output is a single self-contained HTML (`embed-resources: true`); Quarto RevealJS format (needs the Quarto CLI — see Environment Setup).
- CSS goes in the **`text` sub-key of `include-in-header`** in the YAML — injected into `<head>`, embedded by `embed-resources: true`, no blank slide. Never `css: |` inline (broken href), never `{=html}` body block (blank slide), never an external `.css` file (not embedded).
- Use `##` for all slides (flat, no `#` nesting). Section breaks use `{background-color=...}`.
- **Open with two context slides** before findings (recommended default): an **Overview** (the study title and description from `PROJECT_TITLE`/`STUDY_DESCRIPTION` + design variables + live per-group sample sizes + model-system/treatment context; a context slide, no take-home) and a **Sample QC** slide (the report's PCA + live per-sample library sizes from `se$lib.size`; take-home read from the *actual* QC outcome — name any outliers, never assume "all pass"). See `references/slides.md` (Slide Structure → Deck opening).
- Every finding slide leads with a `.takehome` div: one sentence summarising what the data shows. **Do not state whether a finding is safe, risky, concerning, or low-risk** — summarise the biology and evidence only; let the scientist judge.
- All headline numbers are live-read from `RData/` caches in the setup chunk; `basicBulk::fmt_e()` for p-values.
- `basicBulk::assert_fresh()` guards the caches.
- After any layout change: screenshot with **decktape** (not Chrome `#/N`) → a subagent that can read images reads the PNGs → fix → re-sweep until clean (needs Node/`npx`, a Chrome binary, and Ghostscript; see `references/slides.md`).
- Summary slides: findings and biology only — no files/outputs table, no safety conclusions.

## Common Mistakes

Recurring errors when generating basicBulk code. Each is covered in detail in the linked section/reference; this is the quick checklist.

- **`%T>%` binding to the wrong operand.** When a ggplot chain ends in `+ ...` (e.g. `+ ggh4x::facet_nested(...)`, `+ ggplot2::theme(...)`), wrap the whole chain in `(...)` before `%T>% ggplot2::ggsave()`, or `%T>%` binds only the last `+` operand and `ggsave()` fails. See `references/code-style.md`.
- **Missing `results="asis"` / tabset pairing.** Any section whose code calls `purrr::iwalk(basicBulk::printPlots)` or `printHeatmaps` needs both `{.tabset .tabset-pills}` on its header and `results="asis"` on the chunk — these print sub-headers via `cat()`. Applies to heatmaps, LogFC Correlations (Overview + per-pair tabs), GSEA barplots, batch-adjusted PCA. See Section Hierarchy and Chunk Rules.
- **One-tab tabset / orphan pill.** Only add `{.tabset .tabset-pills}` to a section that will have **≥ 2** `##` children. A tabset with a single tab (or a section with one inline output) renders as a lone orphan pill. The DGE section is flat for exactly this reason (its only output is the DEG table). See Section Hierarchy and Chunk Rules.
- **`DGEdata[[.x]]` vs `DGEdata[.x]` for per-contrast plots.** `basicVolcano()`/`basicMA()` take the full `DGEdata`; for a per-contrast plot pass `DGEdata[.x]` (a one-element list), not `DGEdata[[.x]]` (a bare tibble). See `references/output-saving.md`.
- **Re-faceting `basicVolcano()` repeats DEG counts in every panel.** Its per-panel counts are a separate label layer keyed only by `contrast`, so adding `tidyContrasts()` columns to the DGE tibbles isn't enough. Add them to that layer too before `ggh4x::facet_nested()`. See `references/design-patterns.md` (Pattern 4 Volcano and MA Grids).
- **`corTriplot()` erroring on an empty matrix (compound-group designs).** `prepareCor()` factors the stripped contrast tokens against `levels(cd$Group)` (`cd` is required — pass `SummarizedExperiment::colData(se)`). If the contrasts are built on a combined factor (e.g. `sampleGroup` like `"DrugA 50nM"`) but `se$Group` is the bare treatment, the tokens fall outside `levels(Group)`, all cells go `NA`, and `corTriplot()` errors. Pass the matching factor: `corTibble %>% basicBulk::prepareCor(DGEdata = DGEdata, cd = SummarizedExperiment::colData(se) %>% as.data.frame() %>% dplyr::mutate(Group = COMPOUND_GROUP_FACTOR))`. See `references/optional-sections.md` (LogFC Correlations).
- **Treating `runVoom(lfc.threshold = ...)` as the DEG cutoff.** `lfc.threshold` controls the limma treat test (a different statistical test), not the `logFC`/`FDR` thresholds used by `makeSigTibble()`/volcano/MA. Do not repurpose it for custom cutoffs. See Cutoff Parameters.
- **Hardcoding `"logCPM"` downstream, or omitting `assay =`.** Visualization passes `assay = plotAssay` (which may be `"ComBat"`); hardcode `"logCPM"` only in the `plotAssay` assignment or an intentional default-vs-batch comparison that shows both assays side by side (the batch PCA iterates `c("logCPM", "ComBat")` on purpose; do not "correct" it to `plotAssay`). Leaving `assay =` out is the same mistake in another form: the basicBulk plotters default to `assay = "logCPM"`, so with `adjustBatch = TRUE` the plot silently drops the ComBat adjustment.
- **Mapped `aes(shape=)`/`aes(size=)` ignored on `dotPlot`-based plots.** `genePlot()` (via `basicBulk::dotPlot()`), gene set score plots included, sets a fixed point shape/size param, which overrides a later `+ ggplot2::aes(shape = ...)`. To encode a secondary factor (e.g. dose), first free the param — `p$layers[[1]]$aes_params$shape <- NULL` — then add the mapping + `scale_shape_manual()`. See `references/optional-sections.md` (Secondary-Factor Aesthetics).
- **Trusting greedy gene-set-score pathways without a post-render check.** `selectDiversePathways()` (greedy gene-membership Jaccard) is the right default first pass when a redundant family (e.g. ribosome/translation) dominates the top enrichments. But diversity ≠ salience: after the report renders, sanity-check the selected set against the top GSEA terms — greedy selection can drop a distinct, biologically-important theme, and an explicit named list is the fix when it does. See `references/optional-sections.md` (Gene Set Scores).
- **Averaged-control / crossed-design contrast-name collisions.** When a custom contrast embeds a second factor (e.g. dose), put the second-factor token *before* `" vs "` in the contrast name (e.g. `"<treatment> <dose> vs Control"`, not `"<treatment> vs Control <dose>"`), or the default `stringr::str_remove(" vs .*")` label-stripping collapses every contrast to the same label. See `references/design-patterns.md` (Custom Design Matrix).
- **Slide deck CSS — three wrong approaches, one right one.** `css: |` YAML URL-encodes as a broken `<link>` href (nothing loads). A `{=html}` body block before the first `##` creates a blank RevealJS slide. An external `.css` file is not embedded by `embed-resources: true` (styles lost when sharing just the HTML). Use `include-in-header: text: |` in the YAML with a `<style>` block — it injects into `<head>` and is embedded. See `references/slides.md` (CSS section).
- **Slide deck inline p-values render broken.** A raw numeric in scientific notation (e.g. `4.02e-07`) embedded via an inline `` `r x` `` expression gets knitr's LaTeX math formatter and renders as 4.02^{-7}. Pre-format with `basicBulk::fmt_e()` (`formatC(signif(x, 2), format = "e", digits = 1)`) and use the string. See `references/slides.md` (Setup Chunk).

## Testing Notes

`skills/basicbulk/tests/test_patterns.R` validates that the Pattern 1–4 design recipes produce the correct number of contrasts on mock SummarizedExperiments (it fits them with `augere.de::runVoom()`, so it needs augere.de installed). It is a check on recipe correctness, not on whether this skill's behavioral rules (question order, label stripping, save-by-default, design-pattern confirmation) are followed — those are exercised by editing the skill test-first: author the edit with `skill-creator:skill-creator`, then follow `superpowers:writing-skills` (a fresh subagent runs a realistic scenario on the current skill to record the RED baseline before the edit, and again on the edited skill for GREEN).

## Code Style Rules

Generated code must match the author's tidyverse + magrittr style. These are requirements for output, not preferences.

- Use magrittr pipes only: `%>%`, `%T>%`, `%$%`; never base R `|>`.
- Use `%T>%` for side effects that should not break a chain, especially `ggplot2::ggsave()`.
- Use `magrittr::extract2()` instead of `[[` in pipes and `magrittr::extract()` instead of `[` in pipes.
- Use fully qualified namespaces for package functions: `dplyr::filter()`, `dplyr::mutate()`, `ggplot2::ggplot()`, `SummarizedExperiment::colData()`, `basicBulk::basicVolcano()`, etc.
- The only exceptions are:
  - `factor()` is acceptable in simple setup code when explicitly setting analysis group levels. In tidyverse pipelines and other factor manipulation, prefer `forcats::*`.
  - Base helpers such as `grepl()` are acceptable in shell/Rscript inspection snippets. Generated analysis code should use `stringr::*` for string operations.
- Use `purrr::pmap(~ with(list(...), ...))` for iterating over tibble rows with named access.
- Use `purrr::imap()` for named iteration and `purrr::iwalk()` for side-effect iteration such as printing or saving.
- Use `purrr::set_names()` for list names, not `stats::setNames()` or `rlang::set_names()`.
- Add `.progress = TRUE` on long-running `purrr::map()` and `purrr::pmap()` calls.
- Use `forcats::as_factor()` to preserve order in tidyverse pipelines; use `forcats::fct_rev()`, `forcats::fct_drop()`, and `forcats::fct_relevel()` for factor manipulation.
- Relevel the analysis group so the control is `levels(...)[1]`.
- Use `stringr::` functions for string operations in generated analysis code; do not use base `gsub()`, `sub()`, `grep()`, or `grepl()` outside shell/Rscript inspection snippets.
- Use `glue::glue()` for string interpolation and dynamic names; use `glue::glue_collapse()` for vector collapse.
- Use `tibble::tribble()` for inline data and `tidyr::crossing()` for parameter grids.
- Use `ggplot2::theme_bw()` as the base theme.
- Use `ggsci::scale_fill_d3("category20")` and `ggsci::scale_color_d3()` as default palettes.
- Use `ggpubr::rotate_x_text()` for angled x-axis labels.
- Define `plotAssay` once as `"logCPM"` or `"ComBat"` after normalization, then pass `assay = plotAssay` to every basicBulk expression plot (`PCAgrid()`, `basicHeatmap()`, `genePlot()`, `geneHistogram()`, `getHVGs()`, and `joinFeatures()` for gene set scores). Hardcode `"logCPM"` only in the `plotAssay` assignment or an intentional default-vs-batch comparison that shows both assays (the batch PCA).
- Use `tibble::as_tibble()` to convert S4 objects to tibbles for dplyr operations.
- Use `%ni%` for "not in" when defined with `purrr::negate()` around `%in%`.

Read `references/code-style.md` before generating complex tidyverse pipelines, custom plotting code, saving code, tidy-eval code, or when unsure about edge cases.

## Template Mode

When invoked with `/basicbulk:basicbulk` or `/basicbulk:basicbulk template`, generate a complete `.Rmd` file after inspection and design classification.

### `.Rmd` Structure

```yaml
---
title: "TITLE"
subtitle: "DESCRIPTION"
author:
  - name: "Claude Code basicBulk Skill"
output:
  html_document:
    toc: true
    toc_float: true
---
```

Take the title and subtitle from the one answer to Required Workflow step 2, so the report, `RData/settings.RData`, and the deck Overview carry the same title and description — a stale title copied from another study is a silent provenance error on the rendered report. Write that answer as literal strings in the YAML and as the `PROJECT_TITLE`/`STUDY_DESCRIPTION` assignments in the setup chunk (Section 1). Do not reference the setup variables from the YAML with inline `r`: knitr evaluates the YAML before any chunk runs, so a fresh render stops with "object not found".

Include:

```r
knitr::opts_chunk$set(collapse = TRUE,
                      message = FALSE,
                      warning = FALSE,
                      crop = NULL)
```

### Section Hierarchy and Chunk Rules

- Any section whose code calls `purrr::iwalk(basicBulk::printPlots)` or `purrr::iwalk(basicBulk::printHeatmaps)` must have `{.tabset .tabset-pills}` on its R Markdown header and `results="asis"` on its chunk. These functions generate sub-headers via `cat()`, so both are required for tabs to render. Applies to heatmaps, the LogFC Correlations per-pair gg4way tabs, GSEA barplots, and batch-adjusted PCA.
- Add `{.tabset .tabset-pills}` only to a section that will have **≥ 2** `##` children. A one-tab tabset renders as an orphan pill — keep single-output sections flat (e.g. Normalization, DGE).
- Use `fig.width` and `fig.height` from `references/plot-dimensions.md`; Read it before assigning dimensions.
- Use `fig.show="hide"` for chunks that produce plots as side effects, including `augere.de::runVoom()` and `augere.gsea::runContrast()`.
- Use `include=FALSE` for the heatmap build chunks (Section 6 and the Pathways Heatmap tab): they both run computation and emit `png()`/`pdf()` + ggsave "Saving … image" device messages, which `fig.show="hide"` does not suppress. `include=FALSE` still executes the chunk, so the built objects stay available to the print chunk.
- Set `include=FALSE` for pure computation chunks.

### Collapsible Session Info

End the `.Rmd` with this CSS pattern:

```css
button.btn.collapsed:before { content:'Session Info'; display:block; width:100px; }
button.btn:before { content:'Hide'; display:block; width:40px; }
```

```r
knitr::knit_hooks$set(drop1=function(before, options, envir) {
    if (before) {
        paste('<p>',
              '<button class="btn btn-primary collapsed" data-toggle="collapse" data-target="#ce1">',
              '</button>', '</p>',
              '<div class="collapse" id="ce1">',
              '<div class="card card-body">', sep = "\n")
    } else {
        paste("</div>", "</div>", sep = "\n")
    }
})
```

Then use `{r, drop1=TRUE, results="markup"}` for the `utils::sessionInfo()` chunk.

### Template Generation Process

1. Complete the Data Inspection phase to collect all parameters.
2. Classify the experimental design as Pattern 1, Pattern 2, Pattern 3, Pattern 4, or custom.
3. Read the relevant reference files before generating DGE, GSEA, optional sections, dimensions, or complex style-sensitive code.
4. Generate the `.Rmd` using the Write tool, filling in all project-specific values.
5. Present the file to the user for review.

## Package Details

Read `references/package-reference.md` before generating code that depends on detailed function signatures, parameter defaults, or return values.
