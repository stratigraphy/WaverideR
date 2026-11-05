#' @title Precession constant from one or more single-valued cycles
#'
#' @description Each cycle is read as a named astronomical term and inverted
#' directly: k = 360*3600 / (f_observed - f_term). Uncertainty is propagated by
#' Monte Carlo over the cycle's own uncertainty (frequency resolution or `tol`)
#' and over the secular frequencies at that age.
#'
#' @param cycles Measured cycle periods in kyr (named numeric vector); every name is the astronomical term the cycle
#' is read as, e.g. \code{c("k+s3" = 30.8, "k+g5" = 19.9)}. Terms with k (\code{"k+s3"}, \code{"k+g5"},
#' \code{"k+g2"}, \code{"k+g4"}, \code{"k+g3"}, ...) give the precession constant; terms without k (e.g.
#' \code{"e405"}, \code{"g2-g5"}) only constrain the secular frequencies (element \code{secular} of the output).
#' \code{\link{astro_terms}} lists the names that can be used. Required (no default).
#' @param age_ma Age of the record in Ma (numeric): one value, or \code{c(min, max)} to draw the age uniformly
#' from that range in the Monte Carlo simulation. Required (no default).
#' @param record_length_kyr Length of the analysed record in kyr (numeric, > 0); with it the uncertainty band of
#' each cycle is set by the frequency resolution \code{mult / record_length_kyr} (see \code{\link{cycle_band}}).
#' \code{NULL} uses the relative half-width \code{tol} instead. \code{Default=NULL}.
#' @param tol Relative half-width of the uncertainty band of each cycle, used only when \code{record_length_kyr}
#' is \code{NULL} (numeric, > 0; 0.05 is plus or minus 5\%). \code{NULL} gives 0.05 (see \code{\link{cycle_band}}).
#' \code{Default=NULL}.
#' @param mult Multiple of the frequency resolution \code{1 / record_length_kyr} used as half-width of the band
#' (numeric, > 0; 1 is plus or minus one bandwidth). \code{Default=1}.
#' @param nsim Number of Monte Carlo draws (integer, >= 0); within its band each cycle is drawn uniformly.
#' \code{0} gives the plain point estimate without uncertainty. \code{Default=20000}.
#' @param method How the precession constant is converted to Earth-Moon distance, length of day and obliquity
#' (character, passed to \code{\link{ctest_earthmoon}} and \code{\link{earthmoon_state}}): \code{"angular_momentum"}
#' (conservation of Earth-Moon angular momentum) or \code{"astrogeo22"} (inversion of the AstroGeo22 curve of Farhat
#' et al., 2022). \code{Default=c("angular_momentum", "astrogeo22")}, i.e. \code{"angular_momentum"}.
#' @param combine How the estimates of several k terms are combined (character): \code{"multiply"} (product of the
#' individual densities, as in Zeeden et al., 2023, via \code{\link{ctest_combine}}) or \code{"pool"} (all
#' samples together). \code{Default=c("multiply", "pool")}, i.e. \code{"multiply"}.
#' @param seed Seed of the random number generator (integer), so that the Monte Carlo results can be reproduced;
#' \code{NULL} leaves the current random state unchanged. \code{Default=NULL}.
#' @param tuned_to If the record was tuned, the period in kyr that was assumed for the tuning target (numeric, e.g.
#' 405). The true duration of that target at this age is uncertain, so every measured period is stretched by
#' (true target period)/\code{tuned_to} to propagate that error. \code{NULL} treats the input periods as true
#' durations. \code{Default=NULL}.
#' @param anchor_term Astronomical term of the tuning target (character; a secular difference without k, e.g.
#' \code{"e405"} or \code{"g2-g5"}); only used when \code{tuned_to} is given. \code{Default="e405"}.
#' @param sampleorb Shape of the distribution from which the uncertain secular frequencies g2, g3, g4, g5 and s3
#' are drawn (character): \code{"U"} uniform or \code{"N"} normal. \code{Default="U"}.
#' @param g4g3_range_kyr Range \code{c(min, max)} of the period of the g4-g3 term in kyr (numeric vector of
#' length 2), drawn uniformly to allow for a chaotic resonance transition (g3 and g4 are moved symmetrically);
#' \code{NULL} keeps g4-g3 at its La2004 value. \code{Default=NULL}.
#' @param s4s3_range_kyr Range \code{c(min, max)} of the period of the s4-s3 term in kyr (numeric vector of
#' length 2), drawn uniformly as for \code{g4g3_range_kyr}; \code{NULL} keeps s4-s3 at its La2004 value.
#' \code{Default=NULL}.
#' @param s3_uncertainty How the uncertainty of s3 is applied (character): \code{"published"} uses the relative
#' uncertainty directly as standard deviation of the s3 period in years, as in the published code;
#' \code{"scaled"} multiplies it by the s3 period (a much wider spread).
#' \code{Default=c("published", "scaled")}, i.e. \code{"published"}.
#' @param data_dir Folder that contains the data tables \code{AstroGeo22.txt} and \code{La04gisi.csv} (character
#' path). \code{NULL} uses \code{getOption("cocoR.data")} or the copies installed with WaverideR.
#' \code{Default=NULL}.
#' @return list: `per_term` (one row per cycle), `combined` (k and Earth-Moon
#'   state), `samples`, and `consistent` (do the per-term 95\% intervals overlap)
#' @export
k_from_cycles <- function(cycles, age_ma, record_length_kyr = NULL, tol = NULL,
    mult = 1, nsim = 20000, method = c("angular_momentum", "astrogeo22"),
    combine = c("multiply", "pool"), seed = NULL,
    tuned_to = NULL, anchor_term = "e405",
    sampleorb = "U", g4g3_range_kyr = NULL, s4s3_range_kyr = NULL,
    s3_uncertainty = c("published", "scaled"), data_dir = NULL){
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

  # Build n-by-9 sampled secular-frequency matrices (arcsec/yr).
  # g2, g3, g4, g5 and s3 carry the Laskar et al. (2011) uncertainty propagated
  # as sqrt(age); every other g_i / s_i is held at its La2004 value.
  # @param g4g3_range_kyr optional c(min, max) period of the g4-g3 term to sample
  #   from, instead of holding it at La2004 - use this to admit a chaotic
  #   resonance transition (g3 and g4 move symmetrically about their mean)
  # @param s4s3_range_kyr the same for the s4-s3 term
  .sample_secular <- function(n, age_yr, sampleorb = "U", s3_uncertainty = "published",
      point = FALSE, data_dir = NULL,
      g4g3_range_kyr = NULL, s4s3_range_kyr = NULL){
    gs <- .la04_gisi(data_dir)
    G <- matrix(gs$g, n, 9, byrow = TRUE)
    S <- matrix(gs$s, n, 9, byrow = TRUE)
    if (point){
      return(list(G = G, S = S))
    }
    sc <- sqrt(age_yr / 1e8)
    rel <- c(g2 = 0.019 / 7.453, g3 = 0.2 / 17.368, g4 = 0.2 / 17.916,
        g5 = 0.000030 / 4.257482, s3 = 0.066 / 18.848)
    draw <- function(mu, r){
      if (sampleorb == "N"){
        stats::rnorm(n, mu, abs(mu) * r * sc)
      } else {
        stats::runif(n, mu - abs(mu) * 2 * r * sc, mu + abs(mu) * 2 * r * sc)
      }
    }
    G[, 2] <- draw(gs$g[2], rel[["g2"]])
    G[, 3] <- draw(gs$g[3], rel[["g3"]])
    G[, 4] <- draw(gs$g[4], rel[["g4"]])
    G[, 5] <- draw(gs$g[5], rel[["g5"]])
    # s3 as in the published code: relative uncertainty applied in period-years
    s3_per <- 360 * 3600 / abs(gs$s[3])
    if (s3_uncertainty == "scaled"){
      sd_per <- s3_per * rel[["s3"]] * sc
    } else {
      sd_per <- rel[["s3"]] * sc
    }
    S[, 3] <- -360 * 3600 / stats::rnorm(n, s3_per, sd_per)
    if (!is.null(g4g3_range_kyr)){
      target <- 360 * 3600 / (stats::runif(n, min(g4g3_range_kyr), max(g4g3_range_kyr)) * 1000)
      shift <- target - (G[, 4] - G[, 3])
      G[, 4] <- G[, 4] + shift / 2
      G[, 3] <- G[, 3] - shift / 2
    }
    if (!is.null(s4s3_range_kyr)){
      target <- 360 * 3600 / (stats::runif(n, min(s4s3_range_kyr), max(s4s3_range_kyr)) * 1000)
      shift <- target - (S[, 4] - S[, 3])
      S[, 4] <- S[, 4] + shift / 2
      S[, 3] <- S[, 3] - shift / 2
    }
    list(G = G, S = S)
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

  # Main code
  method <- match.arg(method)
  combine <- match.arg(combine)
  s3_uncertainty <- match.arg(s3_uncertainty)
  if (!is.null(seed)){
    set.seed(seed)
  }
  if (length(age_ma) == 1){
    age_ma <- c(age_ma, age_ma)
  }
  terms <- names(cycles)
  if (is.null(terms) || any(!nzchar(terms))){
    stop('Name every cycle with its term, e.g. c("k+s3" = 30.8, "k+g5" = 19.9). ',
        'astro_terms() lists what you can use.')
  }
  if (any(terms == "e125")){
    message('"e125" read as the g3-g2 term (130.7 kyr today); a single unresolved ',
        '~125 kyr peak mixes g3-g2 and g4-g2. Use "e130" or "e124" to be explicit.')
  }
  info <- vector("list", length(terms))
  has_k <- logical(length(terms))
  for (j in seq_along(terms)){
    info[[j]] <- .term_frequency(terms[j], p_arcsec = 0, data_dir = data_dir)
    has_k[j] <- info[[j]]$has_k
  }
  if (!any(has_k)){
    message("None of these terms contains the precession constant, so k is not ",
        "estimated - see $secular for what they do constrain.")
  }

  point <- isTRUE(nsim == 0)
  if (point){
    n <- 1L
    t <- mean(age_ma) * 1e6
  } else {
    n <- nsim
    t <- stats::runif(n, age_ma[1] * 1e6, age_ma[2] * 1e6)
  }
  sec <- .sample_secular(n, t, sampleorb, s3_uncertainty, point, data_dir,
      g4g3_range_kyr, s4s3_range_kyr)
  nom <- .sample_secular(1, 0, point = TRUE, data_dir = data_dir)

  # tuning stretch: assumed target period vs its true duration at this age
  stretch <- rep(1, n)
  if (!is.null(tuned_to)){
    fa <- .term_frequency_vec(anchor_term, sec$G, sec$S)
    if (fa$has_k){
      stop("anchor_term must be a secular difference such as 'e405'")
    }
    stretch <- (360 * 3600 / abs(fa$freq) / 1000) / tuned_to
  }
  band_of <- function(v){
    if (point){
      return(c(v, v))
    }
    cycle_band(v, record_length_kyr, tol, mult)
  }
  draw <- function(band){
    if (point){
      return(band[1])
    }
    stats::runif(n, band[1], band[2])
  }

  # Terms without k: they pin the secular frequencies, not the precession
  sec_rows <- NULL
  if (any(!has_k)){
    idx_nok <- which(!has_k)
    sec_list <- vector("list", length(idx_nok))
    for (m in seq_along(idx_nok)){
      i <- idx_nok[m]
      tf <- .term_frequency_vec(terms[i], sec$G, sec$S)
      band <- band_of(cycles[[i]])
      f_obs <- 360 * 3600 / (draw(band) * stretch * 1000)      # observed difference
      pr <- tf$parts                                           # e.g. g2 - g5
      simple <- isTRUE(tf$simple) && length(pr) == 4
      partner <- NA_real_
      solves <- NA_character_
      ref <- NA_real_
      given <- "combination"
      if (simple){
        if (pr[3] == "g"){
          partner <- nom$G[1, as.integer(pr[4])]
        } else {
          partner <- nom$S[1, as.integer(pr[4])]
        }
        solves <- paste0(pr[1], pr[2])
        if (pr[1] == "g"){
          ref <- nom$G[1, as.integer(pr[2])]
        } else {
          ref <- nom$S[1, as.integer(pr[2])]
        }
        given <- paste0(pr[3], pr[4])
      }
      implied <- f_obs + partner
      difference_sd <- NA_real_
      if (!point){
        difference_sd <- stats::sd(f_obs)
      }
      implied_sd <- NA_real_
      if (!(point || !simple)){
        implied_sd <- stats::sd(implied)
      }
      sec_list[[m]] <- data.frame(term = tf$canonical,
          given = given,
          cycle = cycles[[i]], band_lo = band[1], band_hi = band[2],
          difference_arcsec_yr = mean(f_obs),
          difference_sd = difference_sd,
          solves = solves, implied_arcsec_yr = mean(implied),
          implied_sd = implied_sd,
          la2004_arcsec_yr = ref, offset_arcsec_yr = mean(implied) - ref,
          reference_period_kyr = 360 * 3600 / abs(mean(tf$freq)) / 1000,
          row.names = NULL)
    }
    sec_rows <- do.call(rbind, sec_list)
  }

  # Terms with k
  idx_k <- which(has_k)
  per <- vector("list", length(idx_k))
  for (m in seq_along(idx_k)){
    i <- idx_k[m]
    tf <- .term_frequency_vec(terms[i], sec$G, sec$S, p = 0)   # p = 0 -> the offset
    band <- band_of(cycles[[i]])
    f_obs <- 360 * 3600 / (draw(band) * stretch * 1000)
    # a term is k_mult * p + offset, so invert for p
    p <- (f_obs - tf$offset) / tf$k_mult
    p[p <= 0] <- NA
    per[[m]] <- list(term = tf$canonical, cycle = cycles[[i]], band = band, p = p,
        k = 360 * 3600 / p)
  }

  ci <- function(v){
    if (length(v) < 2){
      sd_v <- NA_real_
    } else {
      sd_v <- stats::sd(v, na.rm = TRUE)
    }
    c(mean = mean(v, na.rm = TRUE),
        sd = sd_v,
        stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))
  }
  if (!length(per)){
    return(list(per_term = NULL, secular = sec_rows, combined = NULL,
        samples = numeric(0), consistent = NA, method = method,
        age_ma = age_ma, point = point,
        reference = astrogeo22(mean(age_ma), data_dir = data_dir)))
  }
  per_rows <- vector("list", length(per))
  for (m in seq_along(per)){
    x <- per[[m]]
    st <- ci(x$k / 1000)
    em1 <- earthmoon_state(k_yr = st[["mean"]] * 1000, age_ma = mean(age_ma),
        method = method, data_dir = data_dir)
    per_rows[[m]] <- data.frame(term = x$term, cycle_kyr = x$cycle,
        band_lo_kyr = x$band[1], band_hi_kyr = x$band[2],
        k_kyr = st[["mean"]], k_sd_kyr = st[["sd"]],
        k_lo95 = st[[3]], k_hi95 = st[[4]],
        a_RE = em1$a_RE, lod_solar_h = em1$lod_solar_h,
        row.names = NULL)
  }
  per_term <- do.call(rbind, per_rows)

  ks <- vector("list", length(per))
  for (m in seq_along(per)){
    ks[[m]] <- per[[m]]$k[is.finite(per[[m]]$k)]
  }
  if (point){
    kcomb <- mean(unlist(ks))
  } else if (length(ks) == 1 || combine == "pool"){
    kcomb <- unlist(ks)
  } else {
    cmb <- ctest_combine(ks, n = 2000)
    kcomb <- sample(cmb$x, nsim, replace = TRUE, prob = pmax(cmb$density, 0))
  }
  em <- ctest_earthmoon(kcomb, age_ma = mean(age_ma), method = method, data_dir = data_dir)
  if (nrow(per_term) > 1){
    overlap <- max(per_term$k_lo95) <= min(per_term$k_hi95)
  } else {
    overlap <- TRUE
  }

  list(per_term = per_term, secular = sec_rows, combined = em$summary,
      samples = kcomb, consistent = overlap, method = method, age_ma = age_ma,
      point = point, reference = astrogeo22(mean(age_ma), data_dir = data_dir))
}
