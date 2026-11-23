## utils-slides.R -- report/slide cache freshness + integrity helpers.

#' Assert a cached artifact is present and newer than its source
#'
#' @description
#' A cheap first-line freshness backstop for report/slide caches: stop with an
#' informative error if `artifact` is missing, or if its modification time
#' predates `source`. This catches the common "I re-ran an upstream step and
#' forgot to regenerate the cache" case. It is a *timestamp* check only — for a
#' guarantee that a cache reflects its inputs regardless of file mtimes (e.g. in
#' place edits, out-of-order re-saves), use [cache_key()] + [load_if_fresh()].
#'
#' @param artifact Path to the cached file that must be up to date.
#' @param source Path to the input file the artifact is derived from.
#' @param regen_cmd Optional string shown in the error telling the user how to
#'   regenerate the artifact. Defaults to "the step that creates it".
#' @return Invisibly `TRUE` when the artifact is present and fresh; otherwise
#'   stops with an error.
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' data(airwayDGE)
#' settings <- tempfile(fileext = ".RData")
#' save(se0, file = settings)
#' results <- tempfile(fileext = ".RData")
#' save(airwayDGE, file = results)
#' assert_fresh(results, settings)
#' @export
assert_fresh <- function(artifact, source, regen_cmd = NULL) {
  msg <- if (is.null(regen_cmd)) "the step that creates it" else regen_cmd
  if (!file.exists(artifact))
    stop(sprintf("Missing artifact: %s\nRe-run: %s", artifact, msg), call. = FALSE)
  if (!file.exists(source))
    stop(sprintf("Missing freshness source: %s (needed to check %s)\nRe-run: %s",
                 source, artifact, msg), call. = FALSE)
  if (file.mtime(artifact) < file.mtime(source))
    stop(sprintf("Stale artifact: %s predates %s\nRe-run: %s", artifact, source, msg),
         call. = FALSE)
  invisible(TRUE)
}

#' Format a number in one-significant-digit scientific notation
#'
#' @description
#' Formats a p-value (or any numeric) as e.g. `"1.2e-04"`. Raw numeric p-values
#' render inconsistently inside inline knitr `r` expressions; pre-formatting with
#' `fmt_e()` keeps slide/report values stable.
#'
#' @param x A numeric vector.
#' @return A character vector in scientific notation with one decimal place.
#' @examples
#' fmt_e(0.000123)
#' @export
fmt_e <- function(x) formatC(signif(x, 2), format = "e", digits = 1)

#' Content hash of a cache's invalidating inputs
#'
#' @description
#' Computes a stable digest of the objects a cached computation depends on, for
#' use as the freshness key with [save_cache()] / [load_if_fresh()]. Pass the
#' *actual input objects* (e.g. the assay matrix and colData of the normalized
#' `SummarizedExperiment`, the design matrix, the contrast list, the gene sets),
#' so that any change to a hashed input — an outlier sample removed, a covariate /
#' block / `dc.block` changed, a contrast or gene set edited — changes the key and
#' forces a recompute. Key on the inputs only, not on package versions or function
#' code: when an upgrade could change the results, force a recompute instead.
#' Because the *objects* are hashed (not a hand-maintained list of field names), a
#' newly added input is captured automatically.
#'
#' @details
#' Hash data objects, not environments or closures, which do not hash stably.
#' Key on the inputs that determine the result. The sample `colData`, the gene
#' `rownames`, the model/contrast specification and the gene-set collection catch a
#' dropped sample or a changed covariate/block/`dc.block`, but not revised counts
#' under the same sample and gene names: add the count values when the input counts
#' can change. Avoid hashing a DelayedArray / HDF5-backed assay directly — it may
#' hash on its backing pointer rather than its contents; realize it with
#' `as.matrix()` first.
#'
#' @param ... The invalidating inputs to hash (any R objects and/or strings).
#' @return A single character string (the digest).
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' cache_key(SummarizedExperiment::colData(se0), rownames(se0))
#' @importFrom rlang hash
#' @export
cache_key <- function(...) {
  rlang::hash(list(...))
}

#' Save a cache stamped with its freshness key
#'
#' @description
#' Saves named objects to `file` together with the content-hash `key` produced by
#' [cache_key()], so [load_if_fresh()] can later verify the cache matches its
#' current inputs before loading it.
#'
#' @param file Path to write the `.RData` cache to.
#' @param key The freshness key from [cache_key()].
#' @param ... Named objects to store in the cache (e.g. `GSEAlist = GSEAlist`).
#' @return Invisibly, the `file` path.
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' data(airwayGSEA)
#' key <- cache_key(SummarizedExperiment::colData(se0), rownames(se0))
#' file <- tempfile(fileext = ".RData")
#' save_cache(file, key, GSEAlist = airwayGSEA)
#' @export
save_cache <- function(file, key, ...) {
  objs <- list(...)
  if (length(objs) == 0L)
    stop("save_cache(): supply at least one named object to cache", call. = FALSE)
  if (is.null(names(objs)) || any(!nzchar(names(objs))))
    stop("save_cache(): all objects passed via ... must be named", call. = FALSE)
  if (anyDuplicated(names(objs)) > 0L)
    stop("save_cache(): duplicate object names in ...: ",
         paste(unique(names(objs)[duplicated(names(objs))]), collapse = ", "), call. = FALSE)
  if (".cache_key" %in% names(objs))
    stop("save_cache(): '.cache_key' is a reserved name and cannot be a cached object", call. = FALSE)
  e <- new.env(parent = emptyenv())
  for (nm in names(objs)) assign(nm, objs[[nm]], envir = e)
  assign(".cache_key", key, envir = e)
  save(list = c(names(objs), ".cache_key"), envir = e, file = file)
  invisible(file)
}

#' Load a cache only if its stored key matches the current inputs
#'
#' @description
#' The safe-by-default guarded loader. Loads `file` (written by [save_cache()])
#' into `envir` **only** when its stored freshness key is identical to `key`
#' (the current inputs, from [cache_key()]). In every other case — the file is
#' missing, unreadable, carries no stored key (an old or foreign cache), or its
#' key differs from `key` — it returns `FALSE` **without loading anything**, so
#' the caller recomputes. This makes stale caches fail safe: a change to any
#' hashed input forces a recompute rather than silently reusing an old result.
#'
#' @details
#' The guard is enforced wherever this runs — a report chunk, the R console, an
#' unattended re-render — so freshness does not depend on any external process
#' being "aware" of the change. Chain it across dependent caches (each keyed on
#' its freshly-computed upstream) so staleness cannot propagate through a
#' cache-of-a-cache.
#'
#' @param file Path to the `.RData` cache written by [save_cache()].
#' @param key The current freshness key from [cache_key()].
#' @param envir Environment to load the cached objects into on a hit. Defaults to
#'   the caller's environment.
#' @return `TRUE` (and injects the cached objects into `envir`) on a fresh hit;
#'   otherwise `FALSE` and nothing is loaded.
#' @examples
#' se0 <- augere.de::loadExampleDataset()
#' data(airwayGSEA)
#' key <- cache_key(SummarizedExperiment::colData(se0), rownames(se0))
#' file <- tempfile(fileext = ".RData")
#' save_cache(file, key, GSEAlist = airwayGSEA)
#' load_if_fresh(file, key)
#' @export
load_if_fresh <- function(file, key, envir = parent.frame()) {
  if (!file.exists(file)) return(FALSE)
  tmp <- new.env(parent = emptyenv())
  ok <- tryCatch({ load(file, envir = tmp); TRUE }, error = function(e) FALSE)
  if (!ok) return(FALSE)
  if (!exists(".cache_key", envir = tmp, inherits = FALSE)) return(FALSE)
  if (!identical(get(".cache_key", envir = tmp, inherits = FALSE), key)) return(FALSE)
  for (nm in setdiff(ls(tmp, all.names = TRUE), ".cache_key"))
    assign(nm, get(nm, envir = tmp, inherits = FALSE), envir = envir)
  TRUE
}
