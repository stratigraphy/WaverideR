test_that("synchrosqueezing follows a sedimentation rate change more sharply than the CWT", {
  set.seed(1)
  z <- seq(0, 1600, by = 1)
  Pt <- ifelse(z < 800, 40, 20)
  d <- data.frame(z, sin(cumsum(1 / Pt) * 2 * pi) + rnorm(length(z), 0, 0.3))
  conc <- function(w) {
    j <- which(z > 150 & z < 650)
    band <- abs(w$Period / 40 - 1) <= 0.1
    mean(colSums(w$Power[band, j, drop = FALSE]) / colSums(w$Power[, j, drop = FALSE]))
  }
  cw <- analyse_wavelet(d, lowerPeriod = 12, upperPeriod = 100, dj = 1 / 40, omega_nr = 12)
  ss <- analyse_wavelet_sst(d, lowerPeriod = 12, upperPeriod = 100, dj = 1 / 40, omega_nr = 12)
  expect_gt(conc(ss), conc(cw))
  # ridge period after the change
  j <- which(z > 1000 & z < 1500)
  rid <- ss$Period[max.col(t(ss$Power[, j]), ties.method = "first")]
  expect_lt(abs(median(rid) / 20 - 1), 0.05)
})

test_that("synchrosqueezing with a high omega_nr separates close cycles", {
  z <- 0:1499
  d <- data.frame(z, sin(2 * pi * z / 40) + sin(2 * pi * z / 44 + 1))
  mid <- which(z > 250 & z < 1250)
  ss <- analyse_wavelet_sst(d, lowerPeriod = 20, upperPeriod = 100, dj = 1 / 100, omega_nr = 20)
  avg <- rowMeans(ss$Power[, mid])
  pk <- which(diff(sign(diff(avg))) == -2) + 1
  expect_true(any(abs(ss$Period[pk] / 40 - 1) < 0.03) && any(abs(ss$Period[pk] / 44 - 1) < 0.03))
})
