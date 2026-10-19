#' @title Track sedimentation-rate solutions through eCOCO peak maps
#'
#' @description Tracks sedimentation-rate paths through the per-window
#' correlation peaks of an eCOCO result, following Acycle's legacy tracking
#' routines. Peak locations (row indices of the \code{n_peaks} largest score
#' values per window that pass the correlation, p-value and orbit-count
#' criteria) are found as in eccpeaks.m. \code{mode = "slices"} reproduces
#' ecocotrack.m: the windows are divided into \code{n_slices} depth slices
#' (sharing boundary windows); in every slice a path is started from each
#' trial rate between min(sr) and max(sr) in steps of \code{sr_threshold},
#' following the nearest peak within \code{sr_threshold} (reference = median
#' of the last five tracked rates, ecocotrackn.m method 1). Candidate paths
#' are summarised with eccount.m (mean/median rate, number of windows, mean
#' rho, mean p, summed pCOCO, summed orbit count), sorted, filtered
#' (mean rho >= \code{rho_min}, mean p <= \code{p_track_max}, summed orbit
#' count >= \code{norbit_min}, near-duplicate solutions within
#' \code{sr_threshold} removed) and the solution with the largest summed
#' pCOCO is kept. \code{mode = "single"} reproduces ecocotrackn.m for one
#' path started at \code{initial_sr} (reference \code{"median5"} = method 1,
#' \code{"previous"} = method 0).
#' Translated from Acycle v3.1 ecocotrack.m, ecocotrackn.m, eccount.m and
#' eccpeaks.m (Mingsong Li).
#'
#' @param ecoco An eCOCO result list (from \code{ecoco_blocked},
#' \code{ecoco_interleaved} or cocoR \code{ecoco_adaptive}) with elements
#' sr, depth, rho, p_global, pcoco and norbit. \code{Default=NULL}
#' @param sr,depth,rho,p,pcoco,norbit,score Explicit inputs used when
#' \code{ecoco} is NULL: rate grid, window depths and rate-by-window maps
#' (\code{score} is the map searched for peaks; default pcoco).
#' \code{Default=NULL}
#' @param location Optional matrix of peak row indices (n_peaks x windows,
#' NA = no peak), replacing the eccpeaks step. \code{Default=NULL}
#' @param mode "slices" (ecocotrack) or "single" (ecocotrackn).
#' \code{Default="slices"}
#' @param n_peaks Number of peaks per window. \code{Default=3}
#' @param rho_min Minimum correlation of a peak and of a path mean.
#' \code{Default=0.3}
#' @param p_max Maximum p-value (same units as \code{p}) of a peak.
#' \code{Default=0.05}
#' @param p_track_max Maximum mean p-value of a path. \code{Default=p_max}
#' @param norbit_min Minimum orbit count of a peak and minimum summed orbit
#' count of a path. \code{Default=3}
#' @param sr_threshold Search radius and trial-rate spacing (cm/kyr).
#' \code{Default=2.5}
#' @param n_slices Number of depth slices. \code{Default=5}
#' @param initial_sr Starting rate for \code{mode = "single"} (cm/kyr).
#' \code{Default=NULL}: the best rate of the first window.
#' @param reference "median5" (median of the previous five tracked rates) or
#' "previous" (previous tracked rate) for \code{mode = "single"}.
#' \code{Default="median5"}
#' @param genplot Plot the tracked rates against depth. \code{Default=FALSE}
#'
#' @author Michiel Arts (R translation); original MATLAB code: Mingsong Li
#' @references Li, M., Kump, L.R., Hinnov, L.A., Mann, M.E., 2018. Tracking
#' variable sedimentation rates and astronomical forcing in Phanerozoic
#' paleoclimate proxy series with evolutionary correlation coefficients and
#' hypothesis testing. Earth and Planetary Science Letters 501, 165-179.
#'
#' @examples
#' \donttest{
#' sr <- seq(1, 10, by = 0.25)
#' depth <- seq(0, 20, by = 0.5)
#' true_sr <- 4 + 2 * depth / 20
#' rho <- matrix(0, length(sr), length(depth))
#' for (j in seq_along(depth)){
#'   rho[, j] <- 0.8 * exp(-(sr - true_sr[j])^2 / 0.5)
#' }
#' p <- matrix(0.01, length(sr), length(depth))
#' nor <- matrix(9, length(sr), length(depth))
#' tr <- ecoco_track(sr = sr, depth = depth, rho = rho, p = p,
#'                   pcoco = 2 * rho, norbit = nor, n_slices = 2)
#' plot(tr$sr, tr$depth)
#' }
#'
#' @return A list with \code{depth}, \code{sr} (tracked rate per window,
#' NA where untracked), \code{index} (row index into the rate grid),
#' \code{location} (peak row-index matrix); for \code{mode = "slices"} also
#' \code{slice_range} (window indices of the slice boundaries),
#' \code{paths} (rate-index paths of all trial rates, rows = trial rates) and
#' \code{summary} (per slice matrix of the eccount statistics: mean index,
#' median index, n, mean rho, mean p, summed pCOCO, summed orbit count;
#' rejected solutions NA).
#' @export
#' @importFrom stats median
#' @importFrom graphics plot
ecoco_track <- function(ecoco = NULL, sr = NULL, depth = NULL, rho = NULL,
                        p = NULL, pcoco = NULL, norbit = NULL, score = NULL,
                        location = NULL, mode = "slices", n_peaks = 3,
                        rho_min = 0.3, p_max = 0.05, p_track_max = p_max,
                        norbit_min = 3, sr_threshold = 2.5, n_slices = 5,
                        initial_sr = NULL, reference = "median5",
                        genplot = FALSE){
  # ---- internal helper functions (local to this function) ----

  # ---------------------------------------------------------------------------
  # MATLAB-compatible scalar helpers
  # ---------------------------------------------------------------------------

  ecocox_round <- function(x){
    # MATLAB round: halves away from zero (R round() rounds half to even)
    sign(x) * floor(abs(x) + 0.5)
  }

  ecocox_max_finite <- function(x){
    # MATLAB max() ignoring NaN; NaN when nothing is finite
    x <- x[is.finite(x)]
    if (length(x) < 1){
      return(NaN)
    }
    max(x)
  }

  ecocox_first_max <- function(v){
    # index of the first finite maximum (MATLAB max over finite entries)
    ok <- which(is.finite(v))
    if (length(ok) < 1){
      return(NA_integer_)
    }
    ok[which.max(v[ok])]
  }

  # ---------------------------------------------------------------------------
  # Legacy peak tracking (eccpeaks.m, ecocotrackn.m, ecocotrack.m, eccount.m)
  # ---------------------------------------------------------------------------

  ecocox_eccpeaks <- function(X, ecc, eci, norbit, rho_min = 0.3, p_max = 5,
                              norbit_min = 3, n = 1){
    # eccpeaks.m (dim = 1, fill value NaN): row indices of the n largest
    # values of every column that pass the rho / p / orbit-count criteria.
    X <- as.matrix(X)
    nc <- ncol(X)
    locat <- matrix(NaN, n, nc)
    for (i in seq_len(nc)){
      xp <- X[, i]
      xs <- sort(xp[!is.na(xp)], decreasing = TRUE)
      for (j in seq_len(n)){
        if (j > length(xs)){
          next
        }
        hit <- which(xp == xs[j])
        if (length(hit) != 1){
          next
        }
        if (isTRUE(ecc[hit, i] >= rho_min && eci[hit, i] <= p_max &&
                   norbit[hit, i] >= norbit_min)){
          locat[j, i] <- hit
        }
      }
    }
    locat
  }

  ecocox_nanmedian <- function(x){
    x <- x[!is.na(x)]
    if (length(x) < 1){
      return(NaN)
    }
    stats::median(x)
  }

  ecocox_trackn <- function(locat, srsh, iniloc, method = 1){
    # ecocotrackn.m
    locat <- as.matrix(locat)
    nrow_l <- nrow(locat)
    ncol_l <- ncol(locat)
    srn1 <- rep(NaN, ncol_l)
    locaj <- locat[, 1]
    for (j in seq_len(nrow_l)){
      if (!is.na(locaj[j]) && abs(locaj[j] - iniloc) < srsh){
        d <- abs(locaj - iniloc)
        srn1[1] <- locaj[which(d == min(d, na.rm = TRUE))][1]
      }
    }
    if (ncol_l < 2){
      return(srn1)
    }
    for (j in 2:ncol_l){
      lt <- locat[, j]
      ltm <- abs(lt - iniloc)
      if (all(is.na(ltm))){
        next
      }
      locn <- min(ltm, na.rm = TRUE)
      if (locn > srsh){
        srn1[j] <- NaN
      } else {
        loc <- which(ltm == locn)
        if (length(loc) > 1){
          ltm <- abs(lt - ecocox_nanmedian(srn1[1:j]))
          if (all(is.na(ltm))){
            next
          }
          locn <- min(ltm, na.rm = TRUE)
          loc <- which(ltm == locn)
          if (length(loc) > 1){
            loc <- loc[1]
          }
        }
        if (length(loc) < 1){
          next
        }
        srn1[j] <- lt[loc]
        if (!is.na(srn1[j - 1])){
          if (method == 0){
            iniloc <- srn1[j - 1]
          } else if (method == 1){
            if (j > 5){
              iniloc <- ecocox_nanmedian(srn1[(j - 5):(j - 1)])
            } else {
              iniloc <- ecocox_nanmedian(srn1[1:(j - 1)])
            }
          }
        }
      }
    }
    srn1
  }

  ecocox_eccount <- function(srn, ecc, eci, ecoco, enorb){
    # eccount.m
    srn <- as.matrix(srn)
    nr <- nrow(srn)
    out <- matrix(0, nr, 7)
    for (i in seq_len(nr)){
      s <- srn[i, ]
      ok <- which(!is.na(s))
      nn <- length(ok)
      m1 <- NaN
      m2 <- NaN
      if (nn > 0){
        m1 <- mean(s[ok])
        m2 <- stats::median(s[ok])
      }
      a <- 0
      b <- 0
      c3 <- 0
      d <- 0
      for (j in ok){
        a <- a + ecc[s[j], j]
        b <- b + eci[s[j], j]
        c3 <- c3 + ecoco[s[j], j]
        d <- d + enorb[s[j], j]
      }
      out[i, ] <- c(m1, m2, nn, a / nn, b / nn, c3, d)
    }
    out
  }

  ecocox_sortrows <- function(M){
    # MATLAB sortrows (ascending, NaN last), all columns as keys
    o <- order(M[, 1], M[, 2], M[, 3], M[, 4], M[, 5], M[, 6], M[, 7],
               na.last = TRUE)
    M[o, , drop = FALSE]
  }

  # ---- main code ----
  mode <- match.arg(mode, c("slices", "single"))
  reference <- match.arg(reference, c("median5", "previous"))
  if (!is.null(ecoco)){
    sr <- ecoco$sr
    depth <- ecoco$depth
    rho <- ecoco$rho
    p <- ecoco$p_global
    pcoco <- ecoco$pcoco
    norbit <- ecoco$norbit
    if (is.null(score)){
      score <- ecoco$pcoco
    }
  }
  if (is.null(score)){
    score <- pcoco
  }
  if (is.null(sr) || is.null(rho) || is.null(p) || is.null(pcoco) ||
      is.null(norbit)){
    stop("Supply an eCOCO result or sr, rho, p, pcoco and norbit")
  }
  rho <- as.matrix(rho)
  if (is.null(depth)){
    depth <- seq_len(ncol(rho))
  }
  sr <- as.numeric(sr)
  srstep <- stats::median(diff(sr))
  if (is.null(location)){
    location <- ecocox_eccpeaks(score, rho, p, norbit, rho_min, p_max,
                                norbit_min, n_peaks)
  }
  location <- as.matrix(location)
  m3 <- ncol(location)
  thr <- sr_threshold / srstep
  if (mode == "single"){
    if (is.null(initial_sr)){
      b <- ecocox_first_max(score[, 1])
      initial_sr <- sr[1]
      if (!is.na(b)){
        initial_sr <- sr[b]
      }
    }
    ini <- (initial_sr - sr[1]) / srstep + 1
    ref_method <- 1
    if (reference == "previous"){
      ref_method <- 0
    }
    idx <- ecocox_trackn(location, thr, ini, ref_method)
    out <- list(depth = depth, sr = sr[1] + srstep * (idx - 1), index = idx,
                location = location)
  } else {
    trial <- seq(sr[1], sr[length(sr)], by = sr_threshold) / srstep
    trial <- ecocox_round(trial)
    trial <- trial - trial[1] + 1
    nt <- length(trial)
    srange <- ecocox_round(seq(1, m3, length.out = n_slices + 1))
    paths <- matrix(NaN, nt, m3)
    best <- rep(NaN, m3)
    summ <- vector("list", n_slices)
    for (s in seq_len(n_slices)){
      c1 <- srange[s]
      c2 <- srange[s + 1]
      loc <- location[, c1:c2, drop = FALSE]
      srn <- matrix(NaN, nt, ncol(loc))
      for (i in seq_len(nt)){
        srn[i, ] <- ecocox_trackn(loc, thr, trial[i], 1)
      }
      paths[, c1:c2] <- srn
      ec <- ecocox_eccount(srn, rho[, c1:c2, drop = FALSE],
                           p[, c1:c2, drop = FALSE],
                           pcoco[, c1:c2, drop = FALSE],
                           norbit[, c1:c2, drop = FALSE])
      ec <- ecocox_sortrows(ec)
      keep <- rep(TRUE, nt)
      for (i in seq_len(nt)){
        if (isTRUE(ec[i, 4] < rho_min)){
          keep[i] <- FALSE
        }
        if (isTRUE(ec[i, 5] > p_track_max)){
          keep[i] <- FALSE
        }
        if (isTRUE(ec[i, 7] < norbit_min)){
          keep[i] <- FALSE
        }
        if (i > 1){
          if (isTRUE(ec[i, 1] - ec[i - 1, 1] < thr)){
            if (isTRUE(ec[i - 1, 6] > ec[i, 6])){
              keep[i] <- FALSE
            } else {
              keep[i - 1] <- FALSE
            }
          }
        }
      }
      ecn <- ec
      ecn[!keep, ] <- NaN
      bmax <- ecocox_max_finite(ecn[, 6])
      bi <- which(ecn[, 6] == bmax)
      if (length(bi) >= 1){
        # Acycle indexes the unsorted path matrix with the sorted summary row
        best[c1:c2] <- srn[bi[1], ]
      } else {
        best[c1:c2] <- NaN
      }
      summ[[s]] <- ecn
    }
    out <- list(depth = depth, sr = sr[1] + srstep * (best - 1), index = best,
                location = location, slice_range = srange, paths = paths,
                summary = summ)
  }
  if (genplot){
    graphics::plot(out$sr, out$depth, xlim = range(sr),
                   xlab = "Sedimentation rate (cm/kyr)", ylab = "Depth (m)",
                   main = "Tracked sedimentation rate")
  }
  out
}
