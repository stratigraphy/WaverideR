#' @title Prepare a depth series for COCO (cocoPrepareRegularData.m)
#' @description Removes non-finite rows, sorts by depth, averages duplicate depths and
#' linearly interpolates to the median step when the spacing is uneven.
#' @param raw two-column matrix [depth (m), value]
#' @param verbose print what was changed
#' @param min_points minimum number of points required
#' @export
coco_prepare_regular_data <- function(raw, verbose = TRUE, min_points = 4){
  # ---- internal helper functions (local to this function) ----

  # --- Interpolation ------------------------------------------------------------
  interp1_na <- function(x, y, xi){
    # Linear interpolation, NA outside range (MATLAB interp1 default).
    o <- order(x)
    x <- x[o]
    y <- y[o]
    keep <- !duplicated(x)
    stats::approx(x[keep], y[keep], xout = xi, method = "linear", rule = 1,
                  ties = "ordered")$y
  }

  # ---- main code ----
  raw <- as.matrix(raw)[, 1:2, drop = FALSE]
  raw <- raw[is.finite(raw[, 1]) & is.finite(raw[, 2]), , drop = FALSE]
  raw <- raw[order(raw[, 1]), , drop = FALSE]
  ud <- sort(unique(raw[, 1]))
  # mean of the values at every unique depth
  groups <- split(raw[, 2], factor(raw[, 1], levels = ud))
  uv <- numeric(length(groups))
  for (i in seq_along(groups)){
    uv[i] <- mean(groups[[i]])
  }
  clean <- cbind(ud, uv)
  if (nrow(clean) < min_points){
    stop("Too few unique finite observations")
  }
  sp <- diff(clean[, 1])
  med <- stats::median(sp)
  tol <- 1e-5 * med
  interp <- max(abs(sp - med)) > tol
  if (interp){
    ex <- (clean[nrow(clean), 1] - clean[1, 1]) / med
    if (abs(ex - round(ex)) <= 1e-10 * max(1, abs(ex))){
      ni <- round(ex)
    } else {
      ni <- floor(ex)
    }
    grid <- clean[1, 1] + (0:ni) * med
    vals <- interp1_na(clean[, 1], clean[, 2], grid)
    ok <- is.finite(vals)
    data <- cbind(grid[ok], vals[ok])
    if (verbose){
      message(sprintf(
        "COCO preprocessing: interpolated %d unique points to %d points at %.6g m (largest gap = %.3g x median)",
        nrow(clean), nrow(data), med, max(sp) / med))
    }
  } else {
    data <- clean
  }
  colnames(data) <- c("depth", "value")
  list(data = data, dt = stats::median(diff(data[, 1])), interpolated = interp)
}
