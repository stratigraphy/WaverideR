#' @title Fixed sinusoid-amplitude weights for a COCO target (cocoFixedTargetWeights.m)
#' @param orbit vector of target periods (kyr)
#' @export
coco_fixed_target_weights <- function(orbit){
  # ---- internal helper functions (local to this function) ----

  # Physical groups of a 9-period target by descending rank (cocoOrbitGroups.m)
  coco_orbit_groups <- function(orbit9){
    if (length(orbit9) != 9 || length(unique(orbit9)) != 9){
      stop("orbit9 must contain nine distinct periods")
    }
    o <- order(orbit9, decreasing = TRUE)
    idx <- integer(9)
    idx[o[1]] <- 1L
    idx[o[2:5]] <- 2L
    idx[o[6]] <- 3L
    idx[o[7:9]] <- 4L
    idx
  }

  # ---- main code ----
  w <- rep(0.6, length(orbit))
  if (length(orbit) == 9 && length(unique(orbit)) == 9){
    g <- coco_orbit_groups(orbit)
    w[g == 3] <- 0.8
    w[g <= 2] <- 1.0
  } else {
    w[orbit >= 30 & orbit < 80] <- 0.8
    w[orbit >= 80] <- 1.0
  }
  w
}
