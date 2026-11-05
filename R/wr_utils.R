#' Internal helpers: safe parallel Monte Carlo runs and plot windows
#'
#' These helpers replace the old foreach / doSNOW code and the direct calls to
#' \code{dev.new()} that were spread over the package. They use base R
#' (the \pkg{parallel} package that ships with R) only.
#'
#' Parallel design:
#' \itemize{
#'   \item sequential by default (no cluster is created unless more than one core is requested);
#'   \item the number of cores is capped by the number of tasks, by the physical cores minus one,
#'   by the \code{_R_CHECK_LIMIT_CORES_} environment variable (CRAN checks) and by the free RAM
#'   (a rough estimate of the memory needed per worker);
#'   \item the tasks run in batches, results are combined after every batch, so a large
#'   number of simulations never keeps all intermediate results in memory at once (important for
#'   \code{combine = "sum"}, where only one running total is kept);
#'   \item each task gets its own random seed, derived from the seed of the master session, so
#'   results do not depend on the number of cores and \code{set.seed()} before a call gives
#'   reproducible output;
#'   \item socket (PSOCK) workers on every platform (Linux, macOS and Windows), also in
#'   the RStudio GUI; the cluster is always stopped, also when an error occurs.
#' }
#'
#' Options: \code{waverider.n_cores} (default number of cores when \code{run_multicore = TRUE}),
#' \code{waverider.max_cores} (upper limit, default 4), \code{waverider.cluster_type} (no
#' longer used: socket workers are always used), \code{waverider.dev_new} (default of the \code{dev_new}
#' arguments of the plot functions).
#'
#' @name wr_utils
#' @keywords internal
#' @noRd
NULL

# Environment

wr_in_rstudio <- function() {
  identical(.Platform$GUI, "RStudio") || identical(Sys.getenv("RSTUDIO"), "1")
}

#' Free memory in bytes (NA when it cannot be determined)
#' @noRd
wr_ram_available <- function() {
  out <- NA_real_
  try(silent = TRUE, {
    sysname <- Sys.info()[["sysname"]]
    if (sysname == "Linux" && file.exists("/proc/meminfo")) {
      ml <- readLines("/proc/meminfo", n = 10L, warn = FALSE)
      ln <- ml[grepl("^MemAvailable:", ml)]
      if (length(ln) == 0L) {
        ln <- ml[grepl("^MemFree:", ml)]
      }
      out <- as.numeric(gsub("[^0-9]", "", ln[1])) * 1024
    } else if (sysname == "Darwin") {
      vm <- system2("vm_stat", stdout = TRUE, stderr = FALSE)
      page <- as.numeric(gsub("[^0-9]", "", vm[1]))
      get <- function(key) as.numeric(gsub("[^0-9]", "", vm[grepl(key, vm)][1]))
      out <- (get("Pages free") + get("Pages inactive") + get("Pages speculative")) * page
    } else if (sysname == "Windows") {
      r <- system2("powershell",
          c("-NoProfile", "-Command",
              "\"(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory\""),
          stdout = TRUE, stderr = FALSE)
      out <- as.numeric(gsub("[^0-9]", "", r[1])) * 1024
    }
  })
  if (length(out) != 1L || !is.finite(out) || out <= 0) {
    out <- NA_real_
  }
  out
}

#' Number of cores to use
#'
#' @param n_cores Requested number of cores (integer, >= 1); \code{NULL} takes the option \code{waverider.n_cores},
#'   and if that is unset 1 core, or with \code{run_multicore = TRUE} the detected cores minus one, at most
#'   \code{waverider.max_cores} (default 4). The result is capped by the detected cores, \code{n_tasks},
#'   \code{_R_CHECK_LIMIT_CORES_} (2) and the free memory. \code{Default=NULL}.
#' @param run_multicore Ask for an automatic (safe) number of cores when \code{n_cores} is \code{NULL} (logical, TRUE
#'   or FALSE; kept for backwards compatibility). \code{Default=FALSE}.
#' @param n_tasks Number of tasks (numeric); never more cores than tasks are used. \code{Default=Inf}.
#' @param mem_per_worker Estimated memory needed per worker in bytes (numeric); the number of cores is limited to 60
#'   percent of the free memory divided by this value. \code{NA} (unknown) applies no memory limit.
#'   \code{Default=NA_real_}.
#' @param verbose Print a message when the number of cores is reduced (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @noRd
wr_n_cores <- function(n_cores = NULL, run_multicore = FALSE, n_tasks = Inf,
    mem_per_worker = NA_real_, verbose = FALSE) {
  detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  if (is.na(detected) || detected < 1L) {
    detected <- 1L
  }
  max_cores <- getOption("waverider.max_cores", 4L)

  if (is.null(n_cores)) {
    n_cores <- getOption("waverider.n_cores", NULL)
  }
  if (is.null(n_cores)) {
    n_cores <- if (isTRUE(run_multicore)) { min(max(detected - 1L, 1L), max_cores) } else { 1L }
  }
  n_cores <- as.integer(max(1, floor(n_cores)))
  requested <- n_cores

  n_cores <- min(n_cores, detected, n_tasks)

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk != "false") {
    n_cores <- min(n_cores, 2L)
  }

  if (n_cores > 1L && is.finite(mem_per_worker) && mem_per_worker > 0) {
    ram <- wr_ram_available()
    if (!is.na(ram)) {
      mem_cap <- max(1L, floor(0.6 * ram / mem_per_worker))
      if (mem_cap < n_cores) {
        if (verbose) {
          message("WaverideR: number of cores reduced from ", n_cores, " to ", mem_cap,
              " to stay within the free memory (about ",
              round(ram / 2^30, 1), " GB free, about ",
              round(mem_per_worker / 2^30, 2), " GB needed per worker).")
        }
        n_cores <- mem_cap
      }
    }
  }
  n_cores <- as.integer(max(1L, n_cores))
  if (verbose && n_cores < requested && requested > 1L) {
    message("WaverideR: using ", n_cores, " core(s) instead of the ", requested, " requested.")
  }
  n_cores
}

# size in bytes of the objects a task closure carries along (rough estimate)
wr_env_size <- function(FUN) {
  e <- environment(FUN)
  if (is.null(e) || identical(e, globalenv()) || isNamespace(e)) {
    return(0)
  }
  size <- 0
  try(silent = TRUE, {
    for (nm in ls(e, all.names = TRUE)) {
      size <- size + as.numeric(utils::object.size(get(nm, envir = e)))
    }
  })
  size
}

# the task closure only carries FUN, X and the seed (it is serialised to socket workers)
wr_make_task <- function(FUN, X, seed) {
  force(FUN)
  force(X)
  force(seed)
  function(k) {
    set.seed(seed + k)
    FUN(X[[k]])
  }
}

# Parallel lapply

#' Run a function over many tasks (Monte Carlo), sequentially or in parallel
#'
#' @param X Task inputs (vector or list, e.g. \code{seq_len(n_simulations)}); \code{FUN} is called once per element.
#'   Required, no default.
#' @param FUN Function of one argument (one element of \code{X}) that runs one task (function or name of a function).
#'   Required, no default.
#' @param n_cores Number of cores (integer, >= 1); \code{NULL} decides with \code{wr_n_cores} from the options and
#'   \code{run_multicore}. \code{Default=NULL}.
#' @param run_multicore Backwards compatible switch (logical, TRUE or FALSE): \code{TRUE} uses an automatic safe
#'   number of cores when \code{n_cores} is \code{NULL}. \code{Default=FALSE}.
#' @param combine How the results are combined (character, one of \code{"list"}, \code{"cbind"}, \code{"sum"}):
#'   \code{"list"} a list of all results, \code{"cbind"} the results bound as columns, \code{"sum"} the results added
#'   (only the running total is kept). \code{Default="list"}.
#' @param progress Show a text progress bar (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @param seed Integer seed from which the task seeds are derived (task k uses \code{seed + k}); \code{NULL} draws it
#'   from the current random state of the session (so \code{set.seed()} before the call makes the run reproducible).
#'   The random state of the session is restored afterwards. \code{Default=NULL}.
#' @param mem_per_worker Estimated memory (bytes) a worker needs (numeric); \code{NA} estimates it as twice the size of
#'   the objects in the environment of \code{FUN} (at least 100 MB). \code{Default=NA_real_}.
#' @param batch_size Number of tasks per batch (integer, at least the number of cores); results are combined after
#'   every batch. \code{NULL} uses 4 tasks per core. \code{Default=NULL}.
#' @param type Kept for compatibility (character); has no effect: socket (PSOCK) workers are used on every platform.
#'   \code{Default=getOption("waverider.cluster_type", "auto")}.
#' @param verbose Print notes on the number of cores (logical, TRUE or FALSE). \code{Default=FALSE}.
#' @noRd
wr_lapply <- function(X, FUN, n_cores = NULL, run_multicore = FALSE,
    combine = c("list", "cbind", "sum"), progress = FALSE, seed = NULL,
    mem_per_worker = NA_real_, batch_size = NULL,
    type = getOption("waverider.cluster_type", "auto"), verbose = FALSE) {
  combine <- match.arg(combine)
  FUN <- match.fun(FUN)
  X <- as.list(X)
  n <- length(X)
  if (n == 0L) {
    return(NULL)
  }

  if (is.na(mem_per_worker)) {
    mem_per_worker <- max(2 * wr_env_size(FUN), 100 * 2^20)
  }
  nc <- wr_n_cores(n_cores, run_multicore, n, mem_per_worker, verbose)

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max %/% 2L, 1L)
  }
  # the task seeds must not change the RNG state of the session
  old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)

  task <- wr_make_task(FUN, X, seed)

  # socket (PSOCK) workers on every platform: forked workers do not exist on Windows and are not used, so the
  # results and their reproducibility are the same on Linux, macOS and Windows
  type <- "psock"

  if (nc > 1L) {
    if (is.null(batch_size)) {
      batch_size <- nc * 4L
    }
    batch_size <- max(nc, as.integer(batch_size))
  } else {
    batch_size <- 1L
  }

  cl <- NULL
  if (nc > 1L && type == "psock") {
    cl <- parallel::makeCluster(nc)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  run_batch <- function(idx) {
    if (nc == 1L) {
      res <- vector("list", length(idx))
      for (k in seq_along(idx)) {
        res[k] <- list(task(idx[[k]]))
      }
      res
    } else {
      parallel::parLapply(cl, idx, task)
    }
  }

  merge2 <- function(acc, new) {
    if (combine == "list") {
      c(acc, new)
    } else if (combine == "cbind") {
      add <- do.call(cbind, new)
      if (is.null(acc)) {
        add
      } else {
        cbind(acc, add)
      }
    } else {
      add <- NULL
      for (k in seq_along(new)) {
        if (k == 1L) {
          add <- new[[1L]]
        } else {
          add <- add + new[[k]]
        }
      }
      if (is.null(acc)) {
        add
      } else {
        acc + add
      }
    }
  }

  acc <- if (combine == "list") { list() } else { NULL }
  done <- 0L
  while (done < n) {
    idx <- (done + 1L):min(n, done + batch_size)
    acc <- merge2(acc, run_batch(idx))
    done <- max(idx)
    if (progress) {
      utils::setTxtProgressBar(pb, done)
    }
  }
  acc
}

# Plot windows

#' Default of the dev_new arguments
#'
#' A separate window is opened by default only in the RStudio GUI, where the plot pane is too
#' small for the raster images of the wavelet plots ("figure margins too large", clipped text).
#' In every other front end (plain R, ESS, Positron, knitr, Rscript) the plot is drawn on the
#' current device, like the ggplot2 functions do. Change with
#' \code{options(waverider.dev_new = TRUE/FALSE)}.
#' @noRd
wr_dev_new_default <- function() {
  opt <- getOption("waverider.dev_new", NULL)
  if (!is.null(opt)) {
    return(isTRUE(opt))
  }
  interactive() && wr_in_rstudio()
}

#' Open a new plot window on any platform
#'
#' Replaces the direct calls to \code{dev.new(noRStudioGD = TRUE)}. Nothing is opened when R is
#' not interactive (knitr, Rscript, batch), so no stray Rplots.pdf files or blocked sessions
#' are produced. In an interactive session a separate window is tried first; if that is not
#' possible or the window cannot take mouse input, the platform's own screen device and finally
#' the RStudio plot pane (which also supports \code{identify()} / \code{locator()}) are used,
#' silently, so that the user never sees a device message.
#'
#' @param width,height Size of the new plot window in inches (numeric, > 0). \code{Default=7} for both.
#' @param interactive_input The window is used for mouse input with \code{identify()} / \code{locator()} (logical,
#'   TRUE or FALSE): it must then accept mouse input, and outside an interactive session an error is raised instead
#'   of silently drawing nothing. \code{Default=FALSE}.
#' @noRd
wr_open_device <- function(width = 7, height = 7, interactive_input = FALSE) {
  if (!interactive()) {
    if (interactive_input) {
      stop("This function needs mouse clicks and can only be used in an interactive R session.",
          call. = FALSE)
    }
    return(invisible(FALSE))
  }
  usable <- function() {
    if (grDevices::dev.cur() == 1) {
      return(FALSE)
    }
    !interactive_input || isTRUE(grDevices::dev.interactive(orNone = FALSE))
  }
  try_open <- function(fun) {
    before <- grDevices::dev.cur()
    opened <- tryCatch({
      suppressWarnings(suppressMessages(fun()))
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (opened && !usable()) {
      if (grDevices::dev.cur() != before && grDevices::dev.cur() != 1) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      opened <- FALSE
    }
    opened
  }
  ## 1. a separate window
  if (try_open(function() { grDevices::dev.new(width = width, height = height, noRStudioGD = TRUE) })) {
    return(invisible(TRUE))
  }
  ## 2. the platform's own screen device
  sysname <- Sys.info()[["sysname"]]
  devs <- if (identical(sysname, "Windows")) {
    "windows"
  } else if (identical(sysname, "Darwin")) {
    c("quartz", "X11")
  } else {
    "X11"
  }
  for (d in devs) {
    f <- get0(d, envir = asNamespace("grDevices"), inherits = FALSE)
    if (is.function(f)) {
      g <- if (d == "X11") {
        function() { f(width = width, height = height, type = "Xlib") }
      } else {
        function() { f(width = width, height = height) }
      }
      if (try_open(g)) {
        return(invisible(TRUE))
      }
    }
  }
  ## 3. the RStudio plot pane (or the current device)
  if (try_open(function() { grDevices::dev.new() })) {
    return(invisible(TRUE))
  }
  invisible(FALSE)
}


#' @importFrom stats mvfft
#' @importFrom grDevices pdf
#' @noRd
NULL

# Colours of the power scale of the wavelet plots (low to high power), shared by all scalogram plots.
# color_brewer is the family of the palette ("RColorBrewer", "colorRamps", "grDevices", "viridis" or "scico").
.wr_key_cols <- function(palette_name, color_brewer, n.levels) {
  if (color_brewer == "RColorBrewer") {
    return(rev(grDevices::colorRampPalette(wr_brewer_pal(wr_brewer_info[palette_name, 1], palette_name))(n.levels)))
  }
  if (color_brewer == "colorRamps") {
    return(get(paste0("wr_cr_", palette_name), mode = "function")(n = n.levels))
  }
  if (color_brewer == "grDevices") {
    if (palette_name == "rainbow") {
      return(rev(grDevices::rainbow(n = n.levels, start = 0, end = 0.7)))
    }
    if (palette_name %in% c("heat.colors", "terrain.colors", "topo.colors", "cm.colors")) {
      return(rev(get(palette_name, envir = asNamespace("grDevices"), mode = "function")(n = n.levels)))
    }
    return(.wr_hcl_cols(n.levels, palette_name))
  }
  if (color_brewer == "viridis") {
    return(rev(get(paste0("wr_viridis_", palette_name), mode = "function")(n = n.levels, direction = -1)))
  }
  if (color_brewer == "scico") {
    return(rev(wr_scico(n = n.levels, palette = palette_name, direction = -1)))
  }
  stop("unknown color_brewer: ", color_brewer)
}
