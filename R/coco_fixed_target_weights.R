#' @title Fixed sinusoid-amplitude weights for a COCO target (cocoFixedTargetWeights.m)
#' @param orbit Target periods in kyr (numeric vector). With nine distinct periods the weights follow their rank
#' (longest five, i.e. eccentricity: 1.0; sixth, obliquity: 0.8; shortest three, precession: 0.6); otherwise
#' periods of 80 kyr and longer get 1.0, 30 to 80 kyr 0.8 and shorter ones 0.6. Required, no default.
#' @export
coco_fixed_target_weights <- function(orbit){
  # Internal helper functions (local to this function)

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

  # Main code
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
