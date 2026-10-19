#' The menu of named astronomical terms at a given age
#'
#' @param age_ma age in Ma; the precession constant comes from AstroGeo22
#' @param k_yr optional: your own precession constant (yr) instead of the model
#' @param data_dir Directory with the data files (`NULL` = the data installed with the package).
#' @return data.frame from the longest to the shortest period
#' @export
astro_terms <- function(age_ma = 0, k_yr = NULL, data_dir = NULL){
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
  if (is.null(k_yr)){
    p <- astrogeo22(age_ma, data_dir = data_dir)$p_arcsec
  } else {
    p <- 360 * 3600 / k_yr
  }
  row <- function(term, alias, family, description, order = 1){
    data.frame(term = term, alias = alias, family = family,
               description = description, order = order, stringsAsFactors = FALSE)
  }
  tab <- rbind(
    row("g4-g3", "e2400", "eccentricity", "'2.4 Myr' eccentricity modulation"),
    row("s4-s3", "m1200", "obliquity modulation", "'1.2 Myr' obliquity modulation"),
    row("g1-g5", "e973", "eccentricity (minor)", "Mercury term"),
    row("g2-g1", "e696", "eccentricity (minor)", "Mercury term"),
    row("g2-g5", "e405", "eccentricity", "405 kyr, the metronome"),
    row("2(g2-g5)-(g4-g3)", "e222", "eccentricity (minor)", "~200 kyr family, 6th order", 3),
    row("2(g2-g5)", "e203/e200", "eccentricity (minor)", "harmonic of the 405 kyr cycle", 2),
    row("(g3-g2)-(g2-g5)", "e193", "eccentricity (minor)", "~200 kyr family, 4th order", 2),
    row("2(g2-g5)+(g4-g3)", "e187", "eccentricity (minor)", "~200 kyr family, 6th order", 3),
    row("(g4-g2)-(g2-g5)", "e178", "eccentricity (minor)", "~200 kyr family, 4th order", 2),
    row("s3-s6", "m173", "obliquity modulation", "'173 kyr' obliquity modulation"),
    row("g3-g2", "e130/e125", "eccentricity", "short eccentricity"),
    row("g4-g2", "e124", "eccentricity", "short eccentricity"),
    row("g3-g1", "e110", "eccentricity (minor)", "Mercury term"),
    row("g4-g1", "e105", "eccentricity (minor)", "Mercury term"),
    row("g3-g5", "e99", "eccentricity", "short eccentricity"),
    row("g4-g5", "e95", "eccentricity", "short eccentricity"),
    row("(g2-g5)+(g3-g5)", "e79", "eccentricity (minor)", "405 + 99 kyr combination tone", 2),
    row("(g2-g5)+(g4-g5)", "e77", "eccentricity (minor)", "405 + 95 kyr combination tone", 2),
    row("(g3-g2)+(g3-g5)", "e56", "eccentricity (minor)", "131 + 99 kyr combination tone", 2),
    row("(g3-g2)+(g4-g5)", "e55", "eccentricity (minor)", "131 + 95 kyr combination tone", 2),
    row("(g4-g2)+(g4-g5)", "e54", "eccentricity (minor)", "124 + 95 kyr combination tone", 2),
    row("k+s6", "", "obliquity", "secondary obliquity term"),
    row("k+s3", "obl", "obliquity", "main obliquity term"),
    row("k+s4", "", "obliquity", "secondary obliquity term"),
    row("k+s2", "", "obliquity (minor)", "minor obliquity term"),
    row("k+s1", "", "obliquity (minor)", "minor obliquity term"),
    row("k+g5", "prec", "precession", "longest precession term"),
    row("k+g1", "", "precession (minor)", "Mercury precession term"),
    row("k+g2", "", "precession", "precession"),
    row("k+g3", "", "precession", "precession"),
    row("k+g4", "", "precession", "shortest precession term"),
    row("2(k+g5)", "", "precession harmonic", "half-precession", 2))
  f <- numeric(nrow(tab))
  for (j in seq_len(nrow(tab))){
    f[j] <- .term_frequency(tab$term[j], p, data_dir = data_dir)$freq
  }
  names(f) <- tab$term
  tab$freq_arcsec_yr <- f
  tab$period_kyr <- 360 * 3600 / abs(f) / 1000
  tab$carries_k <- grepl("^k\\+", tab$term)
  # A secular difference is (k+g_i) - (k+g_j): k cancels in its frequency, but
  # it returns in the ratio to a precession term (bundling / AM sidebands).
  tab$k_route <- ifelse(tab$carries_k, "direct (k + this term)",
                        "only via its ratio to a precession term (bundling)")
  tab <- tab[order(-tab$period_kyr), ]
  rownames(tab) <- NULL
  attr(tab, "age_ma") <- age_ma
  attr(tab, "k_kyr") <- 360 * 3600 / p / 1000
  tab
}
