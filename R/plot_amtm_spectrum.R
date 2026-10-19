# amtm_spectrum.R
#
# Adaptive multitaper (aMTM) PSD background estimation and periodic signal
# detection (gamma-test + harmonic F-test), R translation of the Python package
# amtm_signaldetect v1.0.6 (https://github.com/hsalinasGIT/amtm_signaldetect)
# and of the parts of Prieto's 'multitaper' Python package (v1.2.0,
# https://github.com/gaprieto/multitaper) that amtm_signaldetect calls
# (MTSpec with iadapt = 0: dpss normalisation, eigenspec, adaptspec, rspec, ftest).
#
# ---------------------------------------------------------------------------
# Original licences (both MIT):
#
# amtm_signaldetect -- MIT License, Copyright (c) 2024 Hector A Salinas
# multitaper        -- MIT License, Copyright (c) 2022 German A. Prieto
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ---------------------------------------------------------------------------
# The AR(1) background model, the AIC / MERIT model-selection criteria and the
# AR(1) starting value follow the IDL package SPD_MTM (c) Simone Di Matteo 2020,
# Apache License 2.0 (SPEDAS: general/science/spd_mtm/models/ar1/*.pro,
# spd_mtm_modelgof.pro, spd_mtm_fitmodel.pro). They are not part of the Python
# package. The prewhitening option follows the "Python PreWhiten Filter Example"
# distributed with amtm_signaldetect.
#
# Uses from helpers_matlab_equivalents.R: acy_dpss(), acy_check_two_column().
# ---------------------------------------------------------------------------



#' @title Plot an adaptive multitaper spectrum with background and detections
#'
#' @description The \code{\link{plot_amtm_spectrum}} function plots the output of
#' \code{\link{amtm_spectrum}}: (top) the aMTM PSD with all fitted background models,
#' the selected model and its gamma-test threshold; (middle) the harmonic F-test with
#' its threshold; (bottom) the gamma statistic S/B with its threshold and the gamma
#' bands. Detected periodic signals are marked with vertical lines. The frequency axis
#' is in cycles per unit of the input depth/time column. Translated/adapted from the
#' plotting code of the amtm_signaldetect tutorial notebook (Hector A. Salinas).
#'
#' @param amtm Output list of \code{\link{amtm_spectrum}}
#' @param xlim Frequency limits \code{Default=NULL} (full range)
#' @param log_freq Use a logarithmic frequency axis \code{Default=FALSE}
#' @param freq_lab Label of the frequency axis \code{Default="Frequency (cycles/unit)"}
#' @param show_all_models Draw all fitted background models \code{Default=TRUE}
#'
#' @author Michiel Arts (R translation); original Python code: Hector A. Salinas
#' @references Di Matteo, S., Viall, N. M., & Kepko, L. (2021). Journal of Geophysical
#' Research: Space Physics, 126, e2020JA028748. \doi{10.1029/2020JA028748}
#' @examples
#' \donttest{
#' set.seed(2)
#' n <- 512
#' t <- 0:(n - 1)
#' x <- as.numeric(stats::filter(rnorm(n), 0.9, method = "recursive")) +
#'   0.5 * sin(2 * pi * 0.1 * t) + 0.3 * sin(2 * pi * 0.32 * t)
#' res <- amtm_spectrum(data = cbind(t, x))
#' plot_amtm_spectrum(res)
#' }
#' @return No return value; a plot is drawn.
#' @export
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics rect
#' @importFrom graphics points
#' @importFrom grDevices adjustcolor
plot_amtm_spectrum <- function(amtm = NULL,
                               xlim = NULL,
                               log_freq = FALSE,
                               freq_lab = "Frequency (cycles/unit)",
                               show_all_models = TRUE){
  if (is.null(amtm) || is.null(amtm$psd)){
    stop("amtm must be the output of amtm_spectrum()")
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(3, 1), mar = c(4, 4.5, 1.5, 1))
  f <- amtm$freq
  keep <- f > 0
  if (is.null(xlim)){
    xlim <- range(f[keep])
  }
  logarg <- "y"
  if (log_freq == TRUE){
    logarg <- "xy"
  }
  ff <- amtm$fit_freq
  cols <- c(WHT = "grey50", PL = "darkorange", AR1 = "forestgreen", BPL = "purple")
  graphics::plot(f[keep], amtm$psd[keep], type = "l", log = logarg, xlim = xlim,
                 xlab = freq_lab, ylab = "aMTM PSD", col = "black")
  mnames <- colnames(amtm$background)
  leg <- "aMTM PSD"
  lcol <- "black"
  llwd <- 1
  llty <- 1
  if (show_all_models == TRUE){
    for (j in seq_along(mnames)){
      graphics::lines(ff, amtm$background[, j], col = cols[[mnames[j]]], lwd = 1)
      leg <- c(leg, mnames[j])
      lcol <- c(lcol, cols[[mnames[j]]])
      llwd <- c(llwd, 1)
      llty <- c(llty, 1)
    }
  }
  graphics::lines(ff, amtm$background_selected, col = "red", lwd = 2)
  graphics::lines(ff, amtm$background_selected * amtm$gamma_threshold, col = "red",
                  lty = 2)
  leg <- c(leg, paste0("selected: ", amtm$selected_model),
           paste0(100 * amtm$settings$conf_level, "% gamma-test"))
  lcol <- c(lcol, "red", "red")
  llwd <- c(llwd, 2, 1)
  llty <- c(llty, 1, 2)
  if (nrow(amtm$detected) > 0){
    graphics::abline(v = amtm$detected$frequency, col = "blue", lty = 3)
  }
  graphics::legend("topright", legend = leg, col = lcol, lwd = llwd, lty = llty,
                   bty = "n", cex = 0.8)

  logx <- ""
  if (log_freq == TRUE){
    logx <- "x"
  }
  graphics::plot(f[keep], amtm$ftest[keep], type = "l", log = logx, xlim = xlim,
                 xlab = freq_lab, ylab = "Harmonic F-test")
  graphics::abline(h = amtm$f_threshold, col = "red", lty = 2)
  if (length(amtm$f_peaks) > 0){
    graphics::points(ff[amtm$f_peaks], amtm$ftest[amtm$fit_index][amtm$f_peaks],
                     pch = 1, col = "red")
  }
  if (nrow(amtm$detected) > 0){
    graphics::abline(v = amtm$detected$frequency, col = "blue", lty = 3)
  }

  graphics::plot(ff, amtm$gamma, type = "l", log = logx, xlim = xlim,
                 xlab = freq_lab, ylab = "gamma = S / B")
  if (nrow(amtm$gamma_bands) > 0){
    yl <- graphics::par("usr")[3:4]
    for (j in seq_len(nrow(amtm$gamma_bands))){
      x1 <- amtm$gamma_bands$f_start[j]
      x2 <- amtm$gamma_bands$f_end[j]
      if (log_freq == TRUE){
        x1 <- log10(x1)
        x2 <- log10(x2)
      }
      graphics::rect(x1, yl[1], x2, yl[2], border = NA,
                     col = grDevices::adjustcolor("orange", alpha.f = 0.3))
    }
    graphics::lines(ff, amtm$gamma)
  }
  graphics::abline(h = amtm$gamma_threshold, col = "red", lty = 2)
  if (nrow(amtm$detected) > 0){
    graphics::abline(v = amtm$detected$frequency, col = "blue", lty = 3)
  }
  invisible(NULL)
}
