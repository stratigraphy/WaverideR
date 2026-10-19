# Internal helper (padded Hilbert envelope). Not exported: astrochron::hilbert() provides this functionality.
wr_hilbert_padded <- function(data = NULL,
                              demean = TRUE,
                              nr_pad = 100) {

  mean_dat <- mean(data[, 2])
  data[, 2] <- data[, 2] - mean_dat

  step <- (data[2, 1] - data[1, 1])
  top <- seq(from = data[1, 1] - (nr_pad * step),
             by = step ,
             length.out = nr_pad)
  top_dat <- rep(data[1, 2], times = nr_pad)
  top <- cbind(top, top_dat)
  colnames(top) <- colnames(data)

  bottom <- seq(from = data[nrow(data), 1] + step ,
                by = step ,
                length.out = nr_pad)
  bottom_dat <- rep(data[nrow(data), 2], times = nr_pad)
  bottom <- cbind(bottom, bottom_dat)
  colnames(bottom) <- colnames(data)

  my.data <- rbind(top, data, bottom)
  my.hilbert <- as.data.frame(my.data)


  # instantaneous amplitude (modulus of the analytic signal), base R replacement of DecomposeR::inst.pulse
  my.hilbert[, 2] <- wr_hilbert_envelope(my.data[, 2], my.data[, 1])

  if (demean != TRUE) {
    my.hilbert[, 2] <- my.hilbert[, 2]+mean_dat
  }

  dat_2 <- my.hilbert[c((nr_pad + 1):(nrow(data) + nr_pad)),]

  return(dat_2)
}
