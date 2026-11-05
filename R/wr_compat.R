# Base R replacements for functions that WaverideR used to import from other packages.
#
# Replaced packages: DescTools (Closest), Hmisc (approxExtrap), Matrix (rowMeans), matrixStats (row/col
# statistics), reshape2 (melt), truncnorm (rtruncnorm), colorednoise (colored_noise, autocorrelation),
# trapezoid (rtrapezoid), fANCOVA (loess.as), DecomposeR (inst.pulse, Hilbert envelope), rlist,
# RColorBrewer, colorRamps and viridis (palettes). Only what WaverideR needs is implemented; the
# behaviour follows the original functions (see the comments). All functions are internal.

# Nearest value

# DescTools::Closest: the element(s) of x closest to a; with which = TRUE the index/indices. Ties return
# all closest elements. A list is returned if a has more than one element.
wr_closest <- function(x, a, which = FALSE, na.rm = FALSE) {
  one <- function(a1) {
    xx <- x
    if (na.rm) {
      xx <- xx[!is.na(xx)]
    }
    d <- abs(xx - a1)
    mdist <- min(d)
    if (is.na(mdist)) {
      return(NA)
    }
    idx <- abs(d - mdist) < sqrt(.Machine$double.eps)
    if (which) {
      which(idx)
    } else {
      xx[idx]
    }
  }
  if (length(a) == 1L) {
    return(one(a))
  }
  res <- vector("list", length(a))
  for (i in seq_along(a)) {
    res[[i]] <- one(a[i])
  }
  res
}

# Interpolation with linear extrapolation

# Hmisc::approxExtrap: approx() with linear extrapolation outside the range of x. Returns list(x, y).
wr_approx_extrap <- function(x, y, xout, method = "linear", n = 50, rule = 2, f = 0,
    ties = "ordered", na.rm = FALSE) {
  if (is.list(x)) {
    y <- x[[2]]
    x <- x[[1]]
  }
  if (na.rm) {
    d <- !is.na(x + y)
    x <- x[d]
    y <- y[d]
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  d <- !duplicated(x)
  x <- x[d]
  y <- y[d]
  d <- order(x)
  x <- x[d]
  y <- y[d]
  w <- stats::approx(x, y, xout = xout, method = method, n = n, rule = 2, f = f, ties = ties)$y
  r <- range(x)
  d <- xout < r[1]
  if (any(is.na(d))) {
    stop("NAs not allowed in xout")
  }
  if (any(d)) {
    w[d] <- (y[2] - y[1]) / (x[2] - x[1]) * (xout[d] - x[1]) + y[1]
  }
  d <- xout > r[2]
  n <- length(y)
  if (any(d)) {
    w[d] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1]) * (xout[d] - x[n - 1]) + y[n - 1]
  }
  list(x = xout, y = w)
}

# Row / column statistics (matrixStats)

wr_row_sds <- function(x, na.rm = FALSE) {
  x <- as.matrix(x)
  m <- rowMeans(x, na.rm = na.rm)
  n <- if (na.rm) { rowSums(!is.na(x)) } else { rep(ncol(x), nrow(x)) }
  s <- sqrt(rowSums((x - m)^2, na.rm = na.rm) / (n - 1))
  s[n < 2] <- NA_real_
  s
}

wr_row_mins <- function(x, na.rm = FALSE) {
  x <- as.matrix(x)
  out <- x[, 1]
  if (ncol(x) > 1L) {
    for (j in 2:ncol(x)) {
      out <- pmin(out, x[, j], na.rm = na.rm)
    }
  }
  as.numeric(out)
}

wr_row_maxs <- function(x, na.rm = FALSE) {
  x <- as.matrix(x)
  out <- x[, 1]
  if (ncol(x) > 1L) {
    for (j in 2:ncol(x)) {
      out <- pmax(out, x[, j], na.rm = na.rm)
    }
  }
  as.numeric(out)
}

wr_col_mins <- function(x, na.rm = FALSE) {
  x <- as.matrix(x)
  out <- numeric(ncol(x))
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    if (na.rm) {
      v <- v[!is.na(v)]
    }
    out[j] <- if (length(v)) { min(v) } else { Inf }
  }
  out
}

wr_row_quantiles <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = FALSE) {
  x <- as.matrix(x)
  out <- matrix(NA_real_, nrow(x), length(probs))
  for (i in seq_len(nrow(x))) {
    out[i, ] <- stats::quantile(x[i, ], probs = probs, na.rm = na.rm, names = FALSE, type = 7)
  }
  colnames(out) <- paste0(format(100 * probs, trim = TRUE, digits = 7), "%")
  out
}

# Reshape2::melt for a matrix

# Long format of a matrix: Var1 (row), Var2 (column), value. Dimnames are used if present, otherwise the
# integer indices (as reshape2::melt).
wr_melt <- function(m) {
  m <- as.matrix(m)
  r <- if (is.null(rownames(m))) { seq_len(nrow(m)) } else { factor(rownames(m), levels = unique(rownames(m))) }
  cc <- if (is.null(colnames(m))) { seq_len(ncol(m)) } else { factor(colnames(m), levels = unique(colnames(m))) }
  data.frame(Var1 = rep(r, times = ncol(m)), Var2 = rep(cc, each = nrow(m)), value = as.vector(m))
}

# Random numbers

# truncnorm::rtruncnorm: normal(mean, sd) truncated to [a, b] by inversion of the CDF. The tails are handled on
# the log scale / by symmetry, so far-out truncation limits do not give Inf or NaN.
wr_rtruncnorm <- function(n, a = -Inf, b = Inf, mean = 0, sd = 1) {
  if (length(n) > 1L) {
    n <- length(n)
  }
  if (n == 0L) {
    return(NULL)
  }
  a <- rep_len(as.numeric(a), n)
  b <- rep_len(as.numeric(b), n)
  mean <- rep_len(as.numeric(mean), n)
  sd <- rep_len(as.numeric(sd), n)
  alpha <- (a - mean) / sd
  beta <- (b - mean) / sd
  if (any(alpha > beta, na.rm = TRUE)) {
    stop("rtruncnorm: lower limit a is larger than upper limit b")
  }
  # work in the lower half of the distribution (flip if both limits are in the upper tail)
  flip <- alpha > 0
  lo <- ifelse(flip, -beta, alpha)
  hi <- ifelse(flip, -alpha, beta)
  u <- stats::runif(n)
  plo <- stats::pnorm(lo, log.p = TRUE)
  phi <- stats::pnorm(hi, log.p = TRUE)
  # log(p) = log(exp(plo) + u * (exp(phi) - exp(plo))), computed stably
  lp <- phi + log(u + (1 - u) * exp(plo - phi))
  z <- stats::qnorm(lp, log.p = TRUE)
  z <- pmin(pmax(z, lo), hi)
  z <- ifelse(flip, -z, z)
  mean + sd * z
}

# colorednoise::colored_noise: AR(1) noise with mean, standard deviation sd and autocorrelation phi. Uses the same
# random numbers in the same order as the original C++ code.
#' @title Colour noise and palette helpers
#' @description \code{wr_coloured_noise} generates AR(1) red noise (replaces \code{colorednoise::colored_noise});
#' \code{wr_brewer_palettes} lists the available ColorBrewer palette names and their maximum number of colours.
#' @param timesteps Number of values to generate (integer, >= 1; with 1 a single normal value is returned). Required,
#' no default.
#' @param mean Mean of the noise (numeric, in the unit of the series to be simulated). Required, no default.
#' @param sd Standard deviation of the noise (numeric, > 0); the marginal standard deviation of the AR(1) process.
#' Required, no default.
#' @param phi Lag-1 autocorrelation of the AR(1) process (numeric, between 0 and 1, below 1; 0 gives white noise,
#' values close to 1 give strongly red noise). Required, no default.
#' @return \code{wr_coloured_noise}: numeric vector. \code{wr_brewer_palettes}: data frame with the palettes.
#' @examples
#' set.seed(1)
#' x <- wr_coloured_noise(100, mean = 0, sd = 1, phi = 0.8)
#' head(wr_brewer_palettes())
#' @export
#' @rdname wr_coloured_noise
wr_coloured_noise <- function(timesteps, mean, sd, phi) {
  delta <- mean * (1 - phi)
  variance <- sd^2 * (1 - phi^2)
  first <- stats::rnorm(1, mean, sd)
  if (timesteps < 2L) {
    return(first)
  }
  eps <- delta + stats::rnorm(timesteps - 1L, 0, sqrt(variance))
  c(first, as.numeric(stats::filter(eps, phi, method = "recursive", init = first)))
}

# colorednoise::autocorrelation: bias corrected lag-1 autocorrelation
wr_autocorrelation <- function(x, biasCorrection = TRUE) {
  if (biasCorrection) {
    z <- x[!is.na(x)]
    y <- z[!is.nan(z)]
    r <- stats::acf(y, plot = FALSE)[[1]][2]
    r1 <- stats::acf(y[1:round(length(y) / 2)], plot = FALSE, na.action = stats::na.omit)[[1]][2]
    r2 <- stats::acf(y[(round(length(y) / 2) + 1):length(y)], plot = FALSE,
        na.action = stats::na.omit)[[1]][2]
    2 * r - 0.5 * (r1 + r2)
  } else {
    stats::acf(x, plot = FALSE, na.action = stats::na.omit)[[1]][2]
  }
}

# trapezoid::rtrapezoid for the standard trapezoidal distribution (n1 = n3 = 2, alpha = 1, the only case used by
# WaverideR): density rises linearly on [min, mode1], is flat on [mode1, mode2] and falls linearly on
# [mode2, max]. Sampled by inversion of the CDF.
wr_rtrapezoid <- function(n, min = 0, mode1 = 1 / 3, mode2 = 2 / 3, max = 1) {
  a <- min
  b <- mode1
  c2 <- mode2
  d <- max
  h <- 2 / (d + c2 - a - b)
  p1 <- h * (b - a) / 2
  p2 <- p1 + h * (c2 - b)
  u <- stats::runif(n)
  out <- numeric(n)
  i1 <- u < p1
  i2 <- u >= p1 & u < p2
  i3 <- u >= p2
  out[i1] <- a + sqrt(2 * u[i1] * (b - a) / h)
  out[i2] <- b + (u[i2] - p1) / h
  out[i3] <- d - sqrt(2 * (1 - u[i3]) * (d - c2) / h)
  out
}

# Loess with automatic span (fANCOVA::loess.as, one predictor)

wr_loess_as <- function(x, y, degree = 1, criterion = c("aicc", "gcv"), family = c("gaussian", "symmetric"),
    user.span = NULL, ...) {
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (anyNA(x) || anyNA(y)) {
    stop("'x' and 'y' contain missing values")
  }
  if (length(x) != length(y)) {
    stop("'x' and 'y' have different lengths")
  }
  if (length(y) < 3) {
    stop("not enough observations")
  }
  d <- data.frame(x = x, y = y)
  crit <- function(span) {
    mod <- stats::loess(y ~ x, degree = degree, family = family, span = span, data = d, ...)
    trace_l <- mod$trace.hat
    sigma2 <- sum(mod$residuals^2) / (mod$n - 1)
    if (criterion == "aicc") {
      log(sigma2) + 1 + 2 * (2 * (trace_l + 1)) / (mod$n - trace_l - 2)
    } else {
      mod$n * sigma2 / (mod$n - trace_l)^2
    }
  }
  span1 <- if (is.null(user.span)) { stats::optimize(crit, c(0.05, 0.95))$minimum } else { user.span }
  stats::loess(y ~ x, degree = degree, span = span1, family = family, data = d, ...)
}

# Hilbert envelope (DecomposeR::inst.pulse, method "HT", amplitude)

# Instantaneous amplitude = modulus of the analytic signal. Unevenly sampled records (spacing varying by more than
# 8%) are first interpolated to the median spacing and the envelope is interpolated back.
wr_hilbert_envelope <- function(x, t = NULL) {
  env <- function(v) {
    n <- length(v)
    h <- numeric(n)
    if (n %% 2 == 0) {
      h[c(1, n / 2 + 1)] <- 1
      h[2:(n / 2)] <- 2
    } else {
      h[1] <- 1
      h[2:((n + 1) / 2)] <- 2
    }
    Mod(stats::fft(h * stats::fft(v), inverse = TRUE) / n)
  }
  if (is.null(t)) {
    return(env(x))
  }
  dt <- diff(t)
  if (all(abs(dt / stats::median(dt) - 1) < 0.08)) {
    return(env(x))
  }
  tg <- seq(min(t), max(t), by = stats::median(dt))
  eg <- env(stats::approx(t, x, tg)$y)
  stats::approx(tg, eg, t)$y
}

# Palettes (RColorBrewer, colorRamps, viridis)

wr_brewer_info <- data.frame(maxcolors = unname(.wr_brewer_info), row.names = names(.wr_brewer_info))

# colours of the largest ColorBrewer set of a palette (the one WaverideR always used)
wr_brewer_pal <- function(n, name) {
  if (!(name %in% names(.wr_brewer_cols))) {
    stop(name, " is not a valid palette name for brewer.pal")
  }
  cols <- .wr_brewer_cols[[name]]
  if (n < length(cols)) {
    warning("only the largest set of a ColorBrewer palette is available; returning ", length(cols), " colours")
  }
  cols
}

wr_table_ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1) {
  x <- seq(0, 1, length.out = n)
  y <- rep(0, length(x))
  sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
  sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
  y[sill.min:sill.max] <- 1
  base.min <- round((n - 1) * (mid - base / 2)) + 1
  base.max <- round((n - 1) * (mid + base / 2)) + 1
  xi <- base.min:sill.min
  yi <- seq(0, 1, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  xi <- sill.max:base.max
  yi <- seq(1, 0, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  height * y
}

wr_rgb_tables <- function(n, red = c(0.75, 0.25, 1), green = c(0.5, 0.25, 1), blue = c(0.25, 0.25, 1)) {
  rr <- do.call(wr_table_ramp, as.list(c(n, red)))
  gr <- do.call(wr_table_ramp, as.list(c(n, green)))
  br <- do.call(wr_table_ramp, as.list(c(n, blue)))
  grDevices::rgb(rr, gr, br)
}

wr_cr_blue2green <- function(n) wr_rgb_tables(n, red = c(0, 0, 0, 0), green = c(1, 0, 2), blue = c(0, 0, 2))
wr_cr_blue2green2red <- function(n) wr_rgb_tables(n, red = c(0.8, 0.2, 1), green = c(0.5, 0.4, 0.8), blue = c(0.2, 0.2, 1))
wr_cr_blue2red <- function(n) wr_rgb_tables(n, red = c(1, 1, 1), green = c(0.5, 0, 1), blue = c(0, 1, 1))
wr_cr_blue2yellow <- function(n) wr_rgb_tables(n, red = c(1, 0, 2), green = c(1, 0, 2), blue = c(0, 0, 2))
wr_cr_cyan2yellow <- function(n) wr_rgb_tables(n, red = c(1, 0, 2), green = c(0.5, 1, 2), blue = c(0, 0, 2))
wr_cr_green2red <- function(n) wr_rgb_tables(n, red = c(1, 0, 2), green = c(0, 0, 2), blue = c(0, 0, 0, 0))
wr_cr_magenta2green <- function(n) wr_rgb_tables(n, red = c(0, 0, 2), green = c(1, 0, 2), blue = c(0, 0, 2))
wr_cr_matlab.like <- function(n) wr_rgb_tables(n)
wr_cr_matlab.like2 <- function(n) wr_rgb_tables(n, red = c(0.8, 0.2, 1), green = c(0.5, 0.4, 0.8), blue = c(0.2, 0.2, 1))
wr_cr_ygobb <- function(n) {
  rg <- stats::approx(c(0, 0.25, 0.5, 0.75, 1), c(1, 2 / 3, 2 / 3, 1 / 3, 0), n = n)$y
  b <- stats::approx(c(0, 0.25, 0.5, 0.75, 1), c(2 / 3, 2 / 3, 1 / 3, 2 / 3, 0), n = n)$y
  grDevices::rgb(rg, rg, b)
}

# viridis family (viridisLite::viridis): colours interpolated in Lab space from the 256-colour maps
wr_viridis <- function(n, option = "viridis", direction = 1, begin = 0, end = 1) {
  if (n == 0) {
    return(character(0))
  }
  if (direction == -1) {
    tmp <- begin
    begin <- end
    end <- tmp
  }
  fn_cols <- grDevices::colorRamp(.wr_viridis_maps[[option]], space = "Lab", interpolate = "spline")
  cols <- fn_cols(seq(begin, end, length.out = n)) / 255
  grDevices::rgb(cols[, 1], cols[, 2], cols[, 3])
}
wr_viridis_magma <- function(n, direction = 1) wr_viridis(n, "magma", direction)
wr_viridis_inferno <- function(n, direction = 1) wr_viridis(n, "inferno", direction)
wr_viridis_plasma <- function(n, direction = 1) wr_viridis(n, "plasma", direction)
wr_viridis_viridis <- function(n, direction = 1) wr_viridis(n, "viridis", direction)
wr_viridis_cividis <- function(n, direction = 1) wr_viridis(n, "cividis", direction)
wr_viridis_rocket <- function(n, direction = 1) wr_viridis(n, "rocket", direction)
wr_viridis_mako <- function(n, direction = 1) wr_viridis(n, "mako", direction)
wr_viridis_turbo <- function(n, direction = 1) wr_viridis(n, "turbo", direction)

# batlowK (Crameri, F., 2018. Scientific colour maps, Zenodo, doi:10.5281/zenodo.1243862; MIT licence), 256
# colours, kept in the package because it is the default palette of the scalogram plots
.wr_batlowK <- c(
  "#04050A", "#05060C", "#06070E", "#060A11", "#080C13", "#0A0C14", "#0B0E17", "#0C1019", "#0C1119", "#0D121B",
  "#0E131E", "#0E141F", "#101520", "#101722", "#111825", "#111926", "#121A28", "#121B2B", "#131E2C", "#131F2E",
  "#131F31", "#142133", "#142234", "#152537", "#172639", "#17263A", "#18283D", "#192B3F", "#192C40", "#192D43",
  "#1A2E45", "#1B3146", "#1D3247", "#1E334A", "#1F344B", "#1F374C", "#20384E", "#213950", "#223A51", "#243D52",
  "#253E53", "#263F54", "#264055", "#274357", "#284458", "#2A4559", "#2B4659", "#2C475A", "#2C485A", "#2D4A5B",
  "#2E4B5D", "#304C5D", "#314D5E", "#324E5E", "#33505E", "#33515F", "#34525F", "#35525F", "#37535F", "#38545F",
  "#39545F", "#39555F", "#3A575F", "#3B585F", "#3D595F", "#3E595F", "#3E595F", "#3F5A5E", "#3F5B5E", "#405B5E",
  "#415D5E", "#435D5D", "#435E5D", "#445F5D", "#455F5B", "#465F5B", "#465F5B", "#47605A", "#48615A", "#486159",
  "#4A6359", "#4B6359", "#4C6459", "#4C6558", "#4D6558", "#4E6657", "#506657", "#516655", "#526655", "#526754",
  "#536854", "#546853", "#556A52", "#576A52", "#586B52", "#596C51", "#596C51", "#5A6C50", "#5B6D4E", "#5E6D4E",
  "#5F6E4D", "#5F704C", "#60704C", "#61714C", "#63724B", "#65724B", "#66724A", "#667348", "#687348", "#6A7447",
  "#6B7646", "#6C7646", "#6D7746", "#6E7845", "#717845", "#727944", "#727943", "#747A43", "#767A41", "#787B40",
  "#797D40", "#7A7E3F", "#7D7E3F", "#7E7F3F", "#7F7F3E", "#80803E", "#83803D", "#85813D", "#85833B", "#87843B",
  "#8A853A", "#8C853A", "#8D8539", "#8E8639", "#918739", "#928739", "#948839", "#978A39", "#998B38", "#9A8B38",
  "#9D8C38", "#9F8C38", "#A08D38", "#A38D38", "#A58E38", "#A69038", "#A99139", "#AB9139", "#AC9239", "#AE9239",
  "#B19239", "#B2933A", "#B4943A", "#B7943B", "#B9963D", "#BB963D", "#BE973E", "#BF983F", "#C1983F", "#C49940",
  "#C59941", "#C79943", "#C99945", "#CB9A46", "#CC9A46", "#CE9B48", "#D19B4A", "#D29D4C", "#D39D4C", "#D69D4E",
  "#D89E50", "#D89E52", "#DA9F53", "#DC9F54", "#DE9F57", "#DF9F59", "#E09F5A", "#E1A05D", "#E4A05F", "#E5A05F",
  "#E5A161", "#E6A164", "#E7A366", "#EAA367", "#EBA36A", "#EBA46C", "#ECA46D", "#EDA570", "#EEA572", "#F0A573",
  "#F0A576", "#F1A578", "#F2A679", "#F2A67B", "#F2A67E", "#F3A77F", "#F4A781", "#F4A984", "#F6A985", "#F6AA87",
  "#F7AA8B", "#F7AB8C", "#F8AB8E", "#F8AC91", "#F8AC92", "#F8AC93", "#F9AC96", "#F9AD98", "#F9AD99", "#FAAE9B",
  "#FAAE9E", "#FAB09F", "#FAB0A1", "#FCB1A4", "#FCB1A5", "#FCB2A7", "#FCB2AA", "#FCB2AC", "#FCB2AD", "#FDB3B0",
  "#FDB4B2", "#FDB4B3", "#FDB6B6", "#FDB6B7", "#FDB7B8", "#FDB7BA", "#FDB8BD", "#FDB8BF", "#FDB8C0", "#FDB9C3",
  "#FDBAC5", "#FDBAC6", "#FDBBC9", "#FDBBCB", "#FDBDCC", "#FDBECE", "#FDBED1", "#FDBFD2", "#FDBFD6", "#FDC0D8",
  "#FDC0D9", "#FDC1DC", "#FDC3DE", "#FDC3DF", "#FCC4E1", "#FCC5E4", "#FCC5E6", "#FCC5E9", "#FCC6EB", "#FCC7EC",
  "#FAC7EE", "#FAC9F1", "#FACAF3", "#FACBF6", "#FACBF8", "#F9CCF9"
)

# scico palettes are an optional extra (package scico in Suggests); batlowK also works without it
wr_scico <- function(n, palette, direction = 1) {
  if (requireNamespace("scico", quietly = TRUE)) {
    return(scico::scico(n = n, palette = palette, direction = direction))
  }
  if (identical(palette, "batlowK")) {
    cl <- grDevices::colorRampPalette(.wr_batlowK)(n)
    if (direction < 0) {
      cl <- rev(cl)
    }
    return(cl)
  }
  stop("color_brewer = 'scico' with palette '", palette, "' needs the 'scico' package: install.packages('scico')",
      call. = FALSE)
}

#' @export
#' @rdname wr_coloured_noise
wr_brewer_palettes <- function() wr_brewer_info
