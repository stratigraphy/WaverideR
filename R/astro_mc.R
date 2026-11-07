#' @title Monte Carlo consistency test for tuned (time-domain) or depth-domain cycles
#'
#' @param constraints Measured ranges of the cycles you resolved (named list; each element a numeric vector
#'   \code{c(min, max)}). Names are astronomical terms: \code{"e405"}, \code{"e125"} (= e130 and e124), \code{"e95"}
#'   (= e99 and e95), \code{"k+s3"} (or \code{"obl"}, \code{"obliquity"}), \code{"k+g5"}, \code{"k+g2"},
#'   \code{"k+g4"}, \code{"k+g3"}, \code{"prec_short"} (= k+g4 and k+g3), \code{"prec_all"} (all four precession
#'   terms), any other term listed by \code{\link{astro_terms}} (e.g. \code{"g4-g3"}, \code{"s4-s3"}), or a
#'   combination such as \code{"2(g2-g5)"}. Several terms in one band: \code{"k+g4|k+g3"} (all must fit). A peak you
#'   cannot assign: \code{"any:k+g5|k+g2"} (at least one must fit). Units: kyr when \code{domain = "time"},
#'   metres when \code{domain = "depth"}. Required, no default.
#' @param age Approximate age of the record in Ma (numeric): \code{c(min, max)}, or a single value, which is
#'   widened by \code{age_tol} Ma on either side. Required, no default.
#' @param domain Domain of the measured cycles (character): \code{"time"} for a tuned record (cycles in kyr),
#'   \code{"depth"} for a depth-domain spectrum (cycles in metres; \code{scale_range} must then be given).
#'   \code{Default=c("time", "depth")}, i.e. \code{"time"}.
#' @param anchor Period in kyr that was assumed for the tuning target \code{anchor_term} when the record was tuned
#'   (numeric, > 0); with \code{domain = "time"} every predicted period is divided by (sampled period of
#'   \code{anchor_term})/\code{anchor}. \code{NULL} means that the time scale is sampled freely within
#'   \code{scale_range}. Ignored when \code{domain = "depth"}. \code{Default=405}.
#' @param anchor_term Astronomical term the record was tuned to (character, any name accepted in
#'   \code{constraints}). \code{Default="e405"} (= g2-g5).
#' @param scale_range Range sampled uniformly for the scale factor (numeric vector \code{c(min, max)}): the
#'   sedimentation rate in m/kyr when \code{domain = "depth"} (then required), or the allowed stretch of the time
#'   scale when \code{domain = "time"} and \code{anchor = NULL} (then \code{NULL} becomes \code{c(0.9, 1.1)}).
#'   \code{Default=NULL}.
#' @param g4g3_range_kyr,s4s3_range_kyr Not used by the calculation; accepted only for a common interface with
#'   other functions (numeric vector \code{c(min, max)} in kyr or \code{NULL}). \code{Default=NULL}.
#' @param age_tol Half-width in Ma of the age range used when \code{age} is a single value (numeric, >= 0).
#'   \code{Default=5}.
#' @param nsim Number of accepted (consistent) Monte Carlo samples to collect (integer, > 0). \code{Default=2000}.
#' @param ptest Range of precession constants k to sample uniformly, in years (numeric vector \code{c(min, max)}).
#'   \code{Default=c(10000, 28000)}.
#' @param sampleorb How the secular frequencies g2-g5 are sampled (character): \code{"U"} uniform within 2 sigma
#'   (as in Zeeden et al. 2023) or \code{"N"} normal with 1 sigma. Sigma is the Laskar et al. (2011) uncertainty
#'   scaled by sqrt(age/100 Myr). \code{Default=c("U", "N")}, i.e. \code{"U"}.
#' @param s3_uncertainty How the uncertainty of s3 is sampled (only when a constraint involves k+s3; character):
#'   \code{"published"} follows the published CTest code, in which the relative uncertainty is applied as an
#'   absolute one in years (s3 is effectively fixed); \code{"scaled"} applies it relative to the s3 period, as for
#'   g2-g5, giving s3 a realistic spread. \code{Default=c("published", "scaled")}, i.e. \code{"published"}.
#' @param seed Seed of the random number generator (integer), so that the Monte Carlo results can be reproduced;
#'   \code{NULL} leaves the current random state unchanged. \code{Default=NULL}.
#' @param block Number of random draws per iteration of the vectorised search (integer; values below 1000 are
#'   raised to 1000). Larger blocks are faster but use more memory. \code{Default=1e5}.
#' @param max_attempts Maximum total number of draws before the search stops with an error because fewer than
#'   \code{nsim} consistent samples were found (numeric, > 0). \code{Default=5e7}.
#' @param verbose Print a progress message after each block with accepted samples (logical, TRUE or FALSE).
#'   \code{Default=TRUE}.
#' @return object of class c("astro_mc", "ctest"), so print(), plot() and
#'   ctest_earthmoon() work on it
#' @export
astro_mc <- function(constraints, age, domain = c("time", "depth"),
    anchor = 405, anchor_term = "e405", scale_range = NULL,
    g4g3_range_kyr = NULL, s4s3_range_kyr = NULL,
    age_tol = 5, nsim = 2000, ptest = c(10000, 28000),
    sampleorb = c("U", "N"), s3_uncertainty = c("published", "scaled"),
    seed = NULL, block = 1e5, max_attempts = 5e7, verbose = TRUE){
  # Internal helper functions (local to this function)

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

  # Vectorised term frequency over Monte Carlo samples
  # @param G,S n-by-9 matrices of sampled secular frequencies (arcsec/yr)
  # @param p vector of precession frequencies (arcsec/yr), or NA for secular terms
  .term_frequency_vec <- function(term, G, S, p = NA){
    parts <- .parse_term_expr(term)
    col <- function(let, i){
      if (let == "g"){
        G[, as.integer(i)]
      } else {
        S[, as.integer(i)]
      }
    }
    n <- nrow(G)
    k_mult <- 0
    offset <- rep(0, n)
    first_parts <- NULL
    for (i in seq_len(nrow(parts))){
      b <- parts$base[i]
      cf <- parts$coef[i]
      m <- regmatches(b, regexec("^k\\+([gs])([1-9])$", b))[[1]]
      if (length(m) == 3){
        k_mult <- k_mult + cf
        offset <- offset + cf * col(m[2], m[3])
        if (is.null(first_parts)){
          first_parts <- c(m[2], m[3])
        }
      } else {
        m <- regmatches(b, regexec("^([gs])([1-9])-([gs])([1-9])$", b))[[1]]
        offset <- offset + cf * (col(m[2], m[3]) - col(m[4], m[5]))
        if (is.null(first_parts)){
          first_parts <- c(m[2], m[3], m[4], m[5])
        }
      }
    }
    if (k_mult == 0){
      freq <- offset
    } else {
      freq <- k_mult * p + offset
    }
    list(freq = freq,
        has_k = k_mult != 0, k_mult = k_mult, offset = offset,
        parts = first_parts, canonical = .canonical_term(term),
        simple = nrow(parts) == 1 && abs(parts$coef[1]) == 1)
  }

  # @title Expand a constraint name into the astronomical terms it stands for
  # @keywords internal
  # @noRd
  .expand_term <- function(nm){
    nm <- trimws(nm)
    any_of <- grepl("^any:", nm)
    nm <- sub("^any:", "", nm)
    parts <- trimws(strsplit(nm, "\\|")[[1]])
    expanded <- NULL
    for (p in parts){
      expanded <- c(expanded, switch(p,
          "e125" = c("e130", "e124"),
          "e95"  = c("e99", "e95"),
          "prec_short" = c("k+g4", "k+g3"),
          "prec_all"   = c("k+g5", "k+g2", "k+g4", "k+g3"),
          "obl" = "k+s3", "obliquity" = "k+s3",
          p))
    }
    parts <- expanded
    # validates
    for (pp in setdiff(parts, .ASTRO_TERMS)){
      .term_frequency(pp, 0)
    }
    list(terms = parts, any_of = any_of)
  }

  # Astronomical periods at a given age (calculate_orbit9.m)
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

  # Which terms can be named, and what they are

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

  # astro_mc(): the same consistency test, generalised
  #  * works on a TUNED record (cycles already in kyr) as well as on a
  #    depth-domain spectrum (cycles in m)
  #  * takes whatever terms you actually resolved -- one precession cycle is
  #    enough -- instead of requiring the full eccentricity/obliquity/precession
  #    set of the published functions
  #  * on a tuned record the time scale is not free: if you tuned assuming the
  #    405-kyr cycle equals `anchor` kyr, then every measured duration is
  #    stretched by (true g2-g5 period)/(anchor), which is recomputed for every
  #    Monte Carlo draw instead of being sampled
  # ctest() remains the faithful port of the published code; use that to
  # reproduce Zeeden et al. (2023).

  # @title Astronomical terms known to astro_mc() without extra computation
  # @keywords internal
  # @noRd
  .ASTRO_TERMS <- c("e405", "e130", "e124", "e99", "e95",
      "k+s3", "k+g5", "k+g2", "k+g4", "k+g3")

  # Main code
  domain <- match.arg(domain)
  sampleorb <- match.arg(sampleorb)
  s3_uncertainty <- match.arg(s3_uncertainty)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (length(age) == 1){
    age <- c(age - age_tol, age + age_tol)
  }
  if (!length(constraints)){
    stop("Give at least one cycle range")
  }
  cnames <- names(constraints)
  cons <- vector("list", length(cnames))
  for (j in seq_along(cnames)){
    cons[[j]] <- .expand_term(cnames[j])
  }
  rng <- vector("list", length(constraints))
  names(rng) <- names(constraints)
  for (j in seq_along(constraints)){
    rng[[j]] <- sort(as.numeric(constraints[[j]]))
  }
  needs_s3 <- FALSE
  for (j in seq_along(cons)){
    if ("k+s3" %in% cons[[j]]$terms){
      needs_s3 <- TRUE
    }
  }
  if (domain == "depth" && is.null(scale_range)){
    stop("For domain = 'depth' give scale_range = c(min, max) sedimentation rate in m/kyr")
  }
  if (domain == "time" && is.null(anchor) && is.null(scale_range)){
    scale_range <- c(0.9, 1.1)
  }

  g2 <- 173913
  g3 <- 74620
  g4 <- 72338
  g5 <- 304407
  s3 <- 68753
  dg2 <- 0.019 / 7.453
  dg3 <- 0.2 / 17.368
  dg4 <- 0.2 / 17.916
  dg5 <- 0.000030 / 4.257482
  ds3 <- 0.066 / 18.848

  keep <- vector("list", 0)
  got <- 0
  attempts <- 0
  accepted <- 0
  block <- max(1000, block)
  while (got < nsim){
    attempts <- attempts + block
    if (attempts > max_attempts){
      stop("No consistent samples found; widen the ranges, the age or ptest.")
    }
    t <- stats::runif(block, age[1] * 1e6, age[2] * 1e6)
    p <- stats::runif(block, ptest[1], ptest[2])
    sc <- sqrt(t / 1e8)
    if (sampleorb == "N"){
      G2 <- stats::rnorm(block, g2, g2 * dg2 * sc)
      G3 <- stats::rnorm(block, g3, g3 * dg3 * sc)
      G4 <- stats::rnorm(block, g4, g4 * dg4 * sc)
      G5 <- stats::rnorm(block, g5, g5 * dg5 * sc)
    } else {
      G2 <- stats::runif(block, g2 - g2 * 2 * dg2 * sc, g2 + g2 * 2 * dg2 * sc)
      G3 <- stats::runif(block, g3 - g3 * 2 * dg3 * sc, g3 + g3 * 2 * dg3 * sc)
      G4 <- stats::runif(block, g4 - g4 * 2 * dg4 * sc, g4 + g4 * 2 * dg4 * sc)
      G5 <- stats::runif(block, g5 - g5 * 2 * dg5 * sc, g5 + g5 * 2 * dg5 * sc)
    }
    if (s3_uncertainty == "scaled"){
      sd_s3 <- s3 * ds3 * sc
    } else {
      sd_s3 <- ds3 * sc
    }
    if (needs_s3){
      S3 <- stats::rnorm(block, s3, sd_s3)
    } else {
      S3 <- rep(NA_real_, block)
    }

    all_terms <- NULL
    for (j in seq_along(cons)){
      all_terms <- c(all_terms, cons[[j]]$terms)
    }
    extra <- setdiff(all_terms, .ASTRO_TERMS)
    per <- list(
        "e405" = 1 / (1 / G2 - 1 / G5) / 1000,
        "e130" = 1 / (1 / G3 - 1 / G2) / 1000,
        "e124" = 1 / (1 / G4 - 1 / G2) / 1000,
        "e99"  = 1 / (1 / G3 - 1 / G5) / 1000,
        "e95"  = 1 / (1 / G4 - 1 / G5) / 1000,
        "k+s3" = 1 / (1 / p - 1 / S3) / 1000,
        "k+g5" = 1 / (1 / p + 1 / G5) / 1000,
        "k+g2" = 1 / (1 / p + 1 / G2) / 1000,
        "k+g4" = 1 / (1 / p + 1 / G4) / 1000,
        "k+g3" = 1 / (1 / p + 1 / G3) / 1000)
    if (length(extra)){
      # any other named term (g4-g3, s4-s3, s3-s6, k+s4, ...): the four sampled
      # g's and s3 vary, every other secular frequency is held at La2004
      gsn <- .la04_gisi()
      Gm <- matrix(gsn$g, block, 9, byrow = TRUE)
      Sm <- matrix(gsn$s, block, 9, byrow = TRUE)
      Gm[, 2] <- 360 * 3600 / G2
      Gm[, 3] <- 360 * 3600 / G3
      Gm[, 4] <- 360 * 3600 / G4
      Gm[, 5] <- 360 * 3600 / G5
      Sm[, 3] <- -360 * 3600 / S3
      for (tm in extra){
        tf <- .term_frequency_vec(tm, Gm, Sm, p = 360 * 3600 / p)
        per[[tm]] <- 360 * 3600 / abs(tf$freq) / 1000
      }
    }

    if (domain == "time" && !is.null(anchor)){
      # tuning fixes the time scale: measured durations are stretched by
      # (true anchor period) / (period assumed when tuning)
      scale <- per[[.expand_term(anchor_term)$terms[1]]] / anchor
    } else {
      scale <- stats::runif(block, scale_range[1], scale_range[2])
    }

    ok <- rep(TRUE, block)
    for (j in seq_along(cons)){
      lo <- rng[[j]][1]
      hi <- rng[[j]][2]
      terms_j <- cons[[j]]$terms
      hit <- NULL
      for (m in seq_along(terms_j)){
        # predicted cycle in the units of the input: a duration stretched by the time scale (time domain), or a
        # wavelength = period (kyr) x sedimentation rate (m/kyr) (depth domain)
        obs <- if (domain == "depth") { per[[terms_j[m]]] * scale } else { per[[terms_j[m]]] / scale }
        inside <- obs >= lo & obs <= hi
        if (m == 1){
          hit <- inside
        } else if (cons[[j]]$any_of){
          hit <- hit | inside
        } else {
          hit <- hit & inside
        }
      }
      ok <- ok & (hit %in% TRUE)
    }
    idx <- which(ok)
    accepted <- accepted + length(idx)
    if (!length(idx)){
      next
    }
    idx <- idx[seq_len(min(length(idx), nsim - got))]
    keep[[length(keep) + 1]] <- data.frame(
        scale = scale[idx], k = p[idx], g2 = G2[idx], g3 = G3[idx], g4 = G4[idx],
        g5 = G5[idx], s3 = S3[idx], age_yr = t[idx],
        ecc405 = per$e405[idx], ecc130 = per$e130[idx], ecc124 = per$e124[idx],
        ecc99 = per$e99[idx], ecc95 = per$e95[idx], obl_ks3 = per$`k+s3`[idx],
        prec_kg5 = per$`k+g5`[idx], prec_kg2 = per$`k+g2`[idx],
        prec_kg4 = per$`k+g4`[idx], prec_kg3 = per$`k+g3`[idx])
    got <- got + length(idx)
    if (verbose){
      message(sprintf("astro_mc: %d of %d samples (%.4f%% of draws accepted)",
          got, nsim, 100 * accepted / attempts))
    }
  }
  samples <- do.call(rbind, keep)
  samples$SR <- samples$scale      # so ctest print/plot helpers work unchanged
  ci <- function(v){
    if (!any(is.finite(v))){
      return(rep(NA_real_, 4))
    }
    c(range(v, na.rm = TRUE), mean(v, na.rm = TRUE), stats::sd(v, na.rm = TRUE))
  }
  sum_cols <- c("ecc405", "ecc124", "ecc130", "ecc95", "ecc99",
      "obl_ks3", "prec_kg5", "prec_kg2", "prec_kg4",
      "prec_kg3", "g2", "g3", "g4", "g5", "s3", "k")
  sum_labels <- c("405 kyr ecc", "124 kyr ecc", "130 kyr ecc", "95 kyr ecc",
      "99 kyr ecc", "obliquity (k+s3)", "precession k+g5",
      "precession k+g2", "precession k+g4", "precession k+g3",
      "g2 (yr)", "g3 (yr)", "g4 (yr)", "g5 (yr)", "s3 (yr)",
      "precession constant k (yr)")
  rows <- vector("list", length(sum_cols) + 1)
  rows[[1]] <- ci(samples$scale)
  for (j in seq_along(sum_cols)){
    rows[[j + 1]] <- ci(samples[[sum_cols[j]]])
  }
  names(rows) <- c("time scale / sed. rate", sum_labels)
  summ <- as.data.frame(do.call(rbind, rows))
  names(summ) <- c("lower", "upper", "mean", "sd")
  summ <- cbind(quantity = rownames(summ), summ)
  rownames(summ) <- NULL

  out <- list(summary = summ, samples = samples, mode = paste0(domain, "-domain"),
      nsim = nsim, attempts = attempts, acceptance = accepted / attempts,
      ptest = ptest, sampleorb = sampleorb, age_ma = age,
      constraints = constraints, anchor = anchor, domain = domain)
  class(out) <- c("astro_mc", "ctest")
  out
}
