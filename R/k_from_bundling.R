#' @title Precession constant from bundling: precession cycles per eccentricity bundle
#'
#' @description The scale-free way an eccentricity cycle does carry k. Count
#' how many precession couplets make up one eccentricity bundle (or measure the
#' ratio of the two periods, or the sideband spacing of the precession peak
#' relative to its carrier) and invert
#'     p = n * (g_i - g_j) - f_term.
#' No tuning, no sedimentation rate and no time scale are needed.
#'
#' @param n counted precession cycles per bundle (may be fractional)
#' @param bundle_term the eccentricity term forming the bundle ("e405" by
#'   default; also "e125"/"e130", "e124", "e99", "e95", or any gi-gj)
#' @param prec_term which precession (or obliquity) term is being bundled
#' @param age_ma approximate age, for the secular frequencies and the obliquity
#' @param n_sd 1-sigma counting uncertainty on n (default 0.5 cycles);
#'   `nsim = 0` returns the plain point estimate
#' @param nsim number of Monte Carlo draws
#' @param method "angular_momentum" or "astrogeo22"
#' @param sampleorb,g4g3_range_kyr,s4s3_range_kyr passed to the
#'   secular-frequency sampler
#' @param seed optional RNG seed
#' @param data_dir folder with the data tables
#' @return list with the Earth-Moon `summary`, `k_kyr`, `k_sd_kyr`,
#'   `k_ci95_kyr`, `n`, `terms`, `samples`, `predicted_n_model`, `reference`
#'   and `method`
#' @export
k_from_bundling <- function(n, bundle_term = "e405", prec_term = "k+g5", age_ma,
                            n_sd = 0.5, nsim = 20000,
                            method = c("angular_momentum", "astrogeo22"),
                            sampleorb = "U", g4g3_range_kyr = NULL,
                            s4s3_range_kyr = NULL, seed = NULL, data_dir = NULL){
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
  point <- isTRUE(nsim == 0)
  if (point){
    nn <- 1L
    t <- mean(age_ma) * 1e6
  } else {
    nn <- nsim
    t <- stats::runif(nn, age_ma[1] * 1e6, age_ma[2] * 1e6)
  }
  sec <- .sample_secular(nn, t, sampleorb, point = point, data_dir = data_dir,
                         g4g3_range_kyr = g4g3_range_kyr, s4s3_range_kyr = s4s3_range_kyr)
  tb <- .term_frequency_vec(bundle_term, sec$G, sec$S)
  tp <- .term_frequency_vec(prec_term, sec$G, sec$S, p = 0)
  if (!tp$has_k){
    stop("prec_term must contain k (a precession or obliquity term)")
  }
  if (point){
    nvals <- n
  } else {
    nvals <- stats::rnorm(nn, n, n_sd)
  }
  # n = (k_mult_p * p + offset_p) / (k_mult_b * p + offset_b)
  p <- (nvals * tb$offset - tp$offset) / (tp$k_mult - nvals * tb$k_mult)
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
  br <- bundling_ratios(mean(age_ma), data_dir = data_dir)
  row_b <- grep(.term_frequency(bundle_term, 0)$canonical, rownames(br), fixed = TRUE)[1]
  col_p <- grep(.term_frequency(prec_term, 0)$canonical, colnames(br), fixed = TRUE)[1]
  list(summary = em$summary,
       k_kyr = mean(k / 1000, na.rm = TRUE),
       k_sd_kyr = k_sd_kyr,
       k_ci95_kyr = k_ci95_kyr,
       n = n, terms = c(bundle_term, prec_term), samples = k,
       predicted_n_model = br[row_b, col_p],
       reference = astrogeo22(mean(age_ma), data_dir = data_dir), method = method)
}
