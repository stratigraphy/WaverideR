#' @title Precession constant from the RATIO of two cycles in one section
#'
#' @description Scale-free: the ratio of two cycle periods equals the ratio of
#' the two wavelengths, so this works on an untuned depth series just as well
#' as on a tuned record, and no sedimentation rate is needed. Solving
#' T1/T2 = (p + f2)/(p + f1) for the precession frequency p gives
#'     p = (f2 - R f1) / (R - 1),  R = T1/T2,  T1 > T2.
#'
#' @param long,short the two cycles, same units (kyr or m), long one first
#' @param term_long,term_short which terms they are
#' @param age_ma approximate age (Ma), for the secular frequencies and obliquity
#' @param rel_sd relative 1-sigma uncertainty on each measured cycle (default
#'   0.03); or give `record_length` in the same units as the cycles to use the
#'   frequency resolution instead. `nsim = 0` skips the uncertainty entirely
#'   and returns the plain point estimate.
#' @param record_length,mult see cycle_band()
#' @param nsim number of Monte Carlo draws
#' @param method "angular_momentum" or "astrogeo22"
#' @param seed optional RNG seed
#' @param data_dir folder with the data tables
#' @details Either cycle may be an eccentricity term ("e405", "e125", "e95",
#'   ...). Pairing an eccentricity cycle with a precession or obliquity cycle
#'   is often the most robust option in an untuned section, because the
#'   eccentricity terms carry no k and are the best known frequencies.
#' @return list with the Earth-Moon `summary`, `k_kyr`, `k_sd_kyr`,
#'   `k_ci95_kyr`, `ratio`, `samples`, `terms`, `reference` and `method`
#' @export
k_from_ratio <- function(long, short, term_long = "k+s3", term_short = "k+g5",
                         age_ma, rel_sd = 0.03, record_length = NULL, mult = 1,
                         nsim = 20000, method = c("angular_momentum", "astrogeo22"),
                         seed = NULL, data_dir = NULL){
  # ---- internal helper functions (local to this function) ----

  coco_data_dir <- function(){
    if (!is.null(.coco_env$src)){
      src_data <- file.path(dirname(.coco_env$src), "..", "data")
    } else {
      src_data <- NULL
    }
    cand <- c(getOption("cocoR.data"), system.file("extdata", package = "WaverideR"), "data", file.path("cocoR", "data"), src_data)
    for (d in cand){
      if (!is.null(d) && file.exists(file.path(d, "AstroGeo22.txt"))){
        return(d)
      }
    }
    stop("Cannot find AstroGeo22.txt / La04gisi.csv. Reinstall WaverideR or set options(cocoR.data = '/path/to/data').")
  }

  .la04_gisi <- function(data_dir = NULL){
    if (is.null(data_dir)){
      data_dir <- coco_data_dir()
    }
    gs <- readLines(file.path(data_dir, "La04gisi.csv"), warn = FALSE, encoding = "UTF-8")
    gs <- gs[!grepl("%", gs) & nzchar(trimws(gs))]
    gisi <- as.numeric(strsplit(gs[1], ",")[[1]])[1:18]
    list(g = gisi[1:9], s = gisi[10:18])
  }

  # Split a term expression into signed, integer-weighted base terms
  # Accepts "g2-g5", "k+g5", "2(g2-g5)", "(g3-g2)-(g2-g5)",
  # "2(g2-g5)+(g4-g3)", "2(k+g5)" and any alias from .TERM_ALIASES.
  .parse_term_expr <- function(term){
    tm <- gsub("\\s+", "", as.character(term))
    if (tm %in% names(.TERM_ALIASES)){
      tm <- gsub("\\s+", "", .TERM_ALIASES[[tm]])
    }
    # a bare base term needs no splitting
    if (grepl("^[gs][1-9]-[gs][1-9]$", tm) || grepl("^k\\+[gs][1-9]$", tm)){
      return(data.frame(coef = 1, base = tm, stringsAsFactors = FALSE))
    }
    chars <- strsplit(tm, "")[[1]]
    depth <- 0L
    cuts <- integer(0)
    for (i in seq_along(chars)){
      ch <- chars[i]
      if (ch == "("){
        depth <- depth + 1L
      } else if (ch == ")"){
        depth <- depth - 1L
      } else if (ch %in% c("+", "-") && depth == 0L && i > 1 &&
                 !(chars[i - 1] %in% c("(", "+", "-")) &&
                 !grepl("^k$", paste(chars[max(1, i - 1)], collapse = ""))){
        cuts <- c(cuts, i)
      }
    }
    starts <- c(1L, cuts)
    ends <- c(cuts - 1L, length(chars))
    parts <- character(length(starts))
    for (j in seq_along(starts)){
      parts[j] <- paste(chars[starts[j]:ends[j]], collapse = "")
    }
    out <- vector("list", length(parts))
    for (j in seq_along(parts)){
      pt <- parts[j]
      sgn <- 1
      if (startsWith(pt, "-")){
        sgn <- -1
        pt <- substring(pt, 2)
      } else if (startsWith(pt, "+")){
        pt <- substring(pt, 2)
      }
      m <- regmatches(pt, regexec("^([0-9]*)\\*?\\((.+)\\)$", pt))[[1]]
      if (length(m) == 3){
        if (nzchar(m[2])){
          cf <- as.numeric(m[2])
        } else {
          cf <- 1
        }
        base <- m[3]
      } else {
        m2 <- regmatches(pt, regexec("^([0-9]+)\\*?(.+)$", pt))[[1]]
        if (length(m2) == 3){
          cf <- as.numeric(m2[2])
          base <- m2[3]
        } else {
          cf <- 1
          base <- pt
        }
      }
      if (base %in% names(.TERM_ALIASES)){
        base <- gsub("\\s+", "", .TERM_ALIASES[[base]])
      }
      if (!grepl("^[gs][1-9]-[gs][1-9]$", base) && !grepl("^k\\+[gs][1-9]$", base)){
        stop("Cannot read the term '", term, "': the piece '", base, "' is not ",
             "gi-gj, si-sj or k+gi. astro_terms() lists the named ones.")
      }
      out[[j]] <- data.frame(coef = sgn * cf, base = base, stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
  }

  # Frequency (arcsec/yr) of any named or combined astronomical term
  # @param term e.g. "k+g5", "g2-g5", "2(g2-g5)", "(g3-g2)-(g2-g5)", or an alias
  # @param p_arcsec precession frequency; only needed for terms containing k
  # @return list: freq, has_k, k_mult (how many times k enters), offset (the
  #   part that does not involve k), canonical
  .term_frequency <- function(term, p_arcsec = NA, g = NULL, s = NULL, data_dir = NULL){
    if (is.null(g) || is.null(s)){
      gs <- .la04_gisi(data_dir)
      if (is.null(g)){
        g <- gs$g
      }
      if (is.null(s)){
        s <- gs$s
      }
    }
    parts <- .parse_term_expr(term)
    pick <- function(let, i){
      if (let == "g"){
        g[as.integer(i)]
      } else {
        s[as.integer(i)]
      }
    }
    k_mult <- 0
    offset <- 0
    for (i in seq_len(nrow(parts))){
      b <- parts$base[i]
      cf <- parts$coef[i]
      m <- regmatches(b, regexec("^k\\+([gs])([1-9])$", b))[[1]]
      if (length(m) == 3){
        k_mult <- k_mult + cf
        offset <- offset + cf * pick(m[2], m[3])
      } else {
        m <- regmatches(b, regexec("^([gs])([1-9])-([gs])([1-9])$", b))[[1]]
        offset <- offset + cf * (pick(m[2], m[3]) - pick(m[4], m[5]))
      }
    }
    if (k_mult == 0){
      freq <- offset
    } else {
      freq <- k_mult * p_arcsec + offset
    }
    list(freq = freq,
         has_k = k_mult != 0, k_mult = k_mult,
         offset = offset, canonical = .canonical_term(term))
  }

  .canonical_term <- function(term){
    tm <- gsub("\\s+", "", as.character(term))
    if (tm %in% names(.TERM_ALIASES)){
      gsub("\\s+", "", .TERM_ALIASES[[tm]])
    } else {
      tm
    }
  }

  # --- Astronomical periods at a given age (calculate_orbit9.m) --------------
  # Nine principal periods from the La2004 secular frequencies g_i, s_i
  # (La04gisi.csv) and the precession constant k(t) of AstroGeo22
  # (AstroGeo22.txt, column 11, arcsec/yr), both bundled in Acycle
  # code/package/secular_periods and copied to acycle2R/data.
  #   1 g2-g5 (long ecc.)   2 g3-g2   3 g4-g2   4 g3-g5   5 g4-g5 (short ecc.)
  #   6 k+s3 (obliquity)    7 k+g5    8 k+g2    9 k+g4 (precession)
  # The age is rounded to the nearest integer Ma (AstroGeo22 has 1-Myr steps,
  # 0-4429 Ma). g and s are held at their present La2004 values.
  .coco_env <- new.env()

  .coco_env$src <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e){
    NULL
  })

  # --- Which terms can be named, and what they are -----------------------------

  # Alias table: friendly name -> canonical expression
  # Second-order terms (a single difference of two secular frequencies) plus the
  # higher-order combination tones that make the minor eccentricity cycles.
  # The ~200 kyr family follows Hilgen, F.J., Zeeden, C., Laskar, J. (2020)
  # Paleoclimate records reveal elusive ~200-kyr eccentricity cycle for the
  # first time. Glob. Planet. Change 194, 103296 (their Table 2).
  .TERM_ALIASES <- c(
    # main eccentricity (2nd order)
    "e2400" = "g4-g3", "m2400" = "g4-g3",
    "e405" = "g2-g5", "e130" = "g3-g2", "e125" = "g3-g2", "e124" = "g4-g2",
    "e99" = "g3-g5", "e95" = "g4-g5",
    # minor eccentricity with g1 (Mercury)
    "e973" = "g1-g5", "e696" = "g2-g1", "e110" = "g3-g1", "e105" = "g4-g1",
    # the ~200 kyr family (Hilgen et al. 2020, Table 2)
    "e203" = "2(g2-g5)", "e200" = "2(g2-g5)",
    "e193" = "(g3-g2)-(g2-g5)", "e178" = "(g4-g2)-(g2-g5)",
    "e222" = "2(g2-g5)-(g4-g3)", "e187" = "2(g2-g5)+(g4-g3)",
    # ~77 and ~56 kyr combination tones (405 or 130/124 plus a ~100 kyr term)
    "e79" = "(g2-g5)+(g3-g5)", "e77" = "(g2-g5)+(g4-g5)",
    "e56" = "(g3-g2)+(g3-g5)", "e55" = "(g3-g2)+(g4-g5)",
    "e54" = "(g4-g2)+(g4-g5)",
    # obliquity / precession
    "m1200" = "s4-s3", "m173" = "s3-s6",
    "obl" = "k+s3", "obliquity" = "k+s3",
    "prec" = "k+g5", "precession" = "k+g5")

  # ---- main code ----
  method <- match.arg(method)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (length(age_ma) == 1){
    age_ma <- c(age_ma, age_ma)
  }
  if (long <= short){
    stop("`long` must be the longer of the two cycles")
  }
  point <- isTRUE(nsim == 0)
  draw <- function(v){
    if (point){
      return(v)
    }
    if (!is.null(record_length)){
      b <- cycle_band(v, record_length = record_length, mult = mult)
      stats::runif(nsim, b[1], b[2])
    } else {
      stats::rnorm(nsim, v, rel_sd * v)
    }
  }
  T1 <- draw(long)
  T2 <- draw(short)
  gs <- .la04_gisi(data_dir)
  # each term is either a fixed eccentricity difference (no k) or k + f_term
  fterm <- function(tm){
    .term_frequency(tm, p_arcsec = 0, g = gs$g, s = gs$s, data_dir = data_dir)
  }
  t1 <- fterm(term_long)
  t2 <- fterm(term_short)
  R <- T1 / T2
  # F = k_mult * p + offset for each term, and T1/T2 = F2/F1, so
  #   p = (offset2 - R * offset1) / (R * k_mult1 - k_mult2)
  den <- R * t1$k_mult - t2$k_mult
  if (all(abs(den) < 1e-12)){
    stop("Neither ratio member constrains k (both are secular differences); ",
         "pair one with an obliquity or precession term.")
  }
  p <- (t2$offset - R * t1$offset) / den
  p[!is.finite(p) | p <= 0] <- NA
  k <- 360 * 3600 / p
  em <- ctest_earthmoon(k[is.finite(k)], age_ma = mean(age_ma), method = method,
                        data_dir = data_dir)
  if (point){
    k_sd_kyr <- NA_real_
    k_ci95_kyr <- c(NA_real_, NA_real_)
  } else {
    k_sd_kyr <- stats::sd(k / 1000, na.rm = TRUE)
    k_ci95_kyr <- stats::quantile(k / 1000, c(0.025, 0.975), na.rm = TRUE)
  }
  list(summary = em$summary, k_kyr = mean(k / 1000, na.rm = TRUE),
       k_sd_kyr = k_sd_kyr,
       k_ci95_kyr = k_ci95_kyr,
       ratio = mean(R), samples = k, terms = c(term_long, term_short),
       reference = astrogeo22(mean(age_ma), data_dir = data_dir), method = method)
}
