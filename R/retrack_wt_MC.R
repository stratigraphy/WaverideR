#'@title Re-track cycles using a Monte-Carlo simulation
#'
#'@description When analysing multi-proxy records an age-model can be created
#'for each proxy. These age-models can be in general agreement but might also
#'indicate conflicting deposition rates. Picking one age-model out of the all
#'multi-proxy age-models and stating that, that age-model is the best overlooks
#'the information contained within the other proxies and hence a degree of error
#'remains the age-model exists.To combine the multiple age-models all the age
#'models can be averaged out and the uncertainty can be calculated by means
#'of the standard deviation. The result is an age-model which takes into
#'account all the age-models from the proxy records. The averaged out age-model
#'does not take into account any small user errors during  the creation of the
#'individual age-models nor does the averaging take into account the differences
#'between the age-models and how the initial age-model of a certain proxy might
#'be off in certain intervals. the \link[WaverideR]{retrack_wt_MC} mitigates
#'these problems by re-tracking periods of astronomical cycles in the wavelet
#' spectra. The re-tracking is based on the information provided by the
#' age-models constructed from the different proxy records.
#'First a synthetic tracked curve is created by adding up fractions (0-1) of
#' the tracked periods of the different proxy records. This synthetic curve is
#' then used to re-track the period/spectral peaks of an astronomical cycle
#' in a randomly select wavelet scalogram. This process is repeated x times.
#' The result x tracked curves which take into account all the original
#' age-models. From the re-tracked curves one can calculate the mean period and
#' the standard deviation. The resulting  standard deviation is a good
#' indicator of the quality of the imprint of of astronomical cycles in the
#' proxy records. A small standard deviation indicates that given the input
#' of the different tracked cycles similar periods keep on being tracked
#' indicating the an astronomical is well recorded in the proxy records and as
#' such the age-model is very reliable in set interval. A high standard
#' deviation on the other hand means that the tracking results in vastly
#' different periods of the tracked astronomical cycle, as such the quality of
#' the imprint of the astronomical cycle proxy records is poor and hence the
#' age-model is less-reliable in this interval.
#'
#'@param wt_list Wavelet scalograms of the different proxy records (required): a \code{list} of objects created with
#'\code{\link{analyse_wavelet}}, one per proxy, in the same order as the columns of \code{data_track}. In every
#'simulation one of them is picked at random and the cycle is re-tracked in it. \code{Default=NULL}.
#'@param data_track Tracked periods of the cycle in each proxy record (required): numeric matrix with one column per
#'proxy (same order as \code{wt_list}) and one row per value of \code{x_axis}; contains only the period values, not the
#'depth axis. Build it with \code{\link[base]{cbind}} from tracked curves that were interpolated to the same depths.
#'\code{Default=NULL}.
#'@param x_axis Depth/time values of the rows of \code{data_track} (numeric vector, evenly spaced, same length as
#'\code{nrow(data_track)}). Required. \code{Default=NULL}.
#'@param smoothing Smoothing of every re-tracked curve (character vector of length 1 or 2):
#'\itemize{
#'\item \code{"auto"}: loess smoothing with an automatically selected span (\code{\link{loess_auto}}),
#'\item \code{c("loess", span)}: loess smoothing with the given span (0 to 1; \code{smooth} of \code{astrochron::noLow}),
#'\item \code{c("window", length)}: moving average with the given window length in depth/time units (\code{win} of
#'\code{astrochron::mwStats}).
#'}
#'Any other value applies no smoothing. \code{Default=c("auto")}.
#'@param nr_simulations Number of Monte Carlo simulations (re-tracked curves) (integer, > 0). \code{Default=50}.
#'@param seed_nr Seed of the random number generator (integer), so that the random choice of scalogram and of the
#'weights of the synthetic curve can be reproduced. \code{Default=1337}.
#'@param verbose Show the progress of the simulations (logical, TRUE or FALSE). \code{Default=FALSE}.
#'@param genplot Plot the mean re-tracked period (green) with the mean minus (red) and plus (blue) one standard
#'deviation of the frequency (logical). \code{Default=FALSE}.
#'@param keep_editable Keep the graphical parameters of the plot after the function returns, so that extra elements
#'can be added (logical); only used with \code{genplot = TRUE}. \code{Default=FALSE}.
#'@param create_GIF Save one image per simulation, showing the re-tracked curve on the selected scalogram (drawn with
#'\code{\link{plot_wavelet}}), in the folder \code{file_name} and join them into an animated GIF (logical; needs the
#'\pkg{magick} package). \code{Default=FALSE}.
#'@param plot_GIF Display the GIF after it is created (logical); only used with \code{create_GIF = TRUE}.
#'\code{Default=FALSE}.
#'@param width_plt Width of the GIF images in pixels (integer). \code{Default=600}.
#'@param height_plt Height of the GIF images in pixels (integer). \code{Default=450}.
#'@param period_up Factor (numeric, > 1) by which the synthetic tracked curve is multiplied to get the upper boundary of
#'the search window for the nearest spectral peak (passed to \code{\link{completed_series}}). If no peak is found
#'within the window the interpolated value is used. \code{Default=1.5}.
#'@param period_down Factor (numeric, between 0 and 1) by which the synthetic tracked curve is multiplied to get the
#'lower boundary of the search window for the nearest spectral peak (passed to \code{\link{completed_series}}).
#'\code{Default=0.5}.
#'@param plot.COI Draw the cone of influence on the scalograms of the GIF images (logical); only used with
#'\code{create_GIF = TRUE}. \code{Default=TRUE}.
#'@param n.levels Number of colour levels of the scalograms of the GIF images (integer, > 1). \code{Default=100}.
#' @param palette_name Name of the colour palette of the power (colour) scale, a character string. Which names are
#' valid depends on \code{color_brewer}:
#' \itemize{
#' \item \code{color_brewer = "grDevices"}: \code{"rainbow"}, \code{"heat.colors"}, \code{"terrain.colors"}, \code{"topo.colors"}, \code{"cm.colors"},
#' or any of the 115 \code{\link[grDevices]{hcl.colors}} palettes (list them with \code{grDevices::hcl.pals()}), for
#' example \code{"viridis"}, \code{"Plasma"}, \code{"Inferno"}, \code{"YlOrRd"}, \code{"Spectral"}, \code{"Blue-Red"}, \code{"Zissou 1"}.
#' \item \code{color_brewer = "viridis"}: \code{"viridis"}, \code{"magma"}, \code{"inferno"}, \code{"plasma"}, \code{"cividis"}, \code{"rocket"}, \code{"mako"}, \code{"turbo"}.
#' \item \code{color_brewer = "RColorBrewer"}: \code{"BrBG"}, \code{"PiYG"}, \code{"PRGn"}, \code{"PuOr"}, \code{"RdBu"}, \code{"RdGy"}, \code{"RdYlBu"}, \code{"RdYlGn"}, \code{"Spectral"}, \code{"Accent"}, \code{"Dark2"}, \code{"Paired"}, \code{"Pastel1"}, \code{"Pastel2"}, \code{"Set1"}, \code{"Set2"}, \code{"Set3"}, \code{"Blues"}, \code{"BuGn"}, \code{"BuPu"}, \code{"GnBu"}, \code{"Greens"}, \code{"Greys"}, \code{"Oranges"}, \code{"OrRd"}, \code{"PuBu"}, \code{"PuBuGn"}, \code{"PuRd"}, \code{"Purples"}, \code{"RdPu"}, \code{"Reds"}, \code{"YlGn"}, \code{"YlGnBu"}, \code{"YlOrBr"}, \code{"YlOrRd"} (also listed by \code{wr_brewer_palettes()}).
#' \item \code{color_brewer = "colorRamps"}: \code{"blue2green"}, \code{"blue2green2red"}, \code{"blue2red"}, \code{"blue2yellow"}, \code{"cyan2yellow"}, \code{"green2red"}, \code{"magenta2green"}, \code{"matlab.like"}, \code{"matlab.like2"}, \code{"ygobb"}.
#' \item \code{color_brewer = "scico"} (all except \code{"batlowK"} need the \pkg{scico} package): \code{"acton"}, \code{"bam"}, \code{"bamako"}, \code{"bamO"}, \code{"batlow"}, \code{"batlowK"}, \code{"batlowW"}, \code{"berlin"}, \code{"bilbao"}, \code{"broc"}, \code{"brocO"}, \code{"buda"}, \code{"bukavu"}, \code{"cork"}, \code{"corkO"}, \code{"davos"}, \code{"devon"}, \code{"fes"}, \code{"glasgow"}, \code{"grayC"}, \code{"hawaii"}, \code{"imola"}, \code{"lajolla"}, \code{"lapaz"}, \code{"lipari"}, \code{"lisbon"}, \code{"managua"}, \code{"navia"}, \code{"nuuk"}, \code{"oleron"}, \code{"oslo"}, \code{"roma"}, \code{"romaO"}, \code{"tofino"}, \code{"tokyo"}, \code{"turku"}, \code{"vanimo"}, \code{"vik"}, \code{"vikO"}.
#' }
#' Low power is drawn with the first colour of the palette and high power with the last. \code{Default="viridis"}.
#' @param color_brewer Palette family from which \code{palette_name} is taken, a character string: \code{"grDevices"},
#' \code{"viridis"}, \code{"RColorBrewer"}, \code{"colorRamps"} or \code{"scico"} (the last one, except for \code{"batlowK"}, needs the
#' \pkg{scico} package). See \code{palette_name} for the palette names of each family. \code{Default="viridis"}.
#'@param periodlab Label of the period axis of the plots (character). \code{Default="Period (metres)"}.
#'@param x_lab Label of the depth/time axis of the plots (character). \code{Default="depth (metres)"}.
#'@param add_avg Draw the average wavelet power next to the scalograms of the GIF images (logical; passed to
#'\code{\link{plot_wavelet}}); only used with \code{create_GIF = TRUE}. \code{Default=FALSE}.
#'@param time_dir Direction of the depth/time axis in the GIF images (logical): \code{TRUE} if time increases with
#'increasing depth/time values (e.g. a core that gets older downwards), \code{FALSE} if it decreases (e.g. a section
#'log where 0 m is the base); only used with \code{create_GIF = TRUE}. \code{Default=TRUE}.
#'@param file_name Name (character) of the folder created in the working directory for the GIF images, also used as the
#'base name of the images (\code{file_name_<simulation number>.jpeg}) and of the GIF (\code{file_name.gif}). Required
#'when \code{create_GIF = TRUE}. \code{Default=NULL}.
#'@param run_multicore Run the Monte Carlo simulations in parallel on several cores (logical). \code{Default=FALSE}.
#'@param n_cores Number of cores used when \code{run_multicore = TRUE} (integer, or \code{NULL}). \code{NULL} uses an
#'automatic, safe number: physical cores minus one, at most 4, and limited by the free memory; set
#'\code{options(waverider.n_cores = ...)} to change it. A value > 1 also switches on parallel computing. Results do
#'not depend on the number of cores. \code{Default=NULL}.
#'@param output What is returned (integer, 1 to 7):
#'\itemize{
#'\item 1: list of (1) a matrix with the x-axis, the mean tracked frequency and its standard deviation (columns
#'\code{depth}, \code{mean_freq}, \code{sd}), (2) a matrix with all the re-tracked periods (one column per simulation)
#'and (3) the GIF (\code{NULL} without \code{create_GIF = TRUE}),
#'\item 2: list of the mean/sd matrix and the matrix of re-tracked periods,
#'\item 3: list of the mean/sd matrix and the GIF,
#'\item 4: list of the matrix of re-tracked periods and the GIF,
#'\item 5: the mean/sd matrix only,
#'\item 6: the matrix of re-tracked periods only,
#'\item 7: the GIF only.
#'}
#'\code{Default=1}.
#'@param n_imgs Number of images (simulations) put in the GIF (integer, > 0; at most \code{nr_simulations}). Many
#'images take long to process and give a large GIF. \code{Default=50}.
#'@param plot_horizontal Orientation of the scalograms of the GIF images (logical): \code{TRUE} puts depth/time on the
#'x-axis, \code{FALSE} on the y-axis. \code{Default=TRUE}.
#'@param empty_folder Delete all files in the folder \code{file_name} before the new images are written (logical);
#'only used with \code{create_GIF = TRUE}. \code{Default=FALSE}.
#'
#'
#'@return
#'The output depends on the output setting
#'If genplot = TRUE a plot will be generated in which the mean period and
#'standard deviation is plotted
#'if plot_GIF = TRUE a GIF with n number of n_imgs will be plotted in which the
#'retraced curve is plotted in a wavelet scalogram
#'If output = 1, output is a list which contain 3 objects. object 1 is a matrix
#'with the x-axis and the mean tracked frequency and standard deviation.
#'object 2 is a matrix with all the tracked periods. Object 3 is a GIF in which
#'all the tracked periods are plotted.
#'If output = 2, output is a list which contain 2 objects. object 1 is a matrix
#'with the x-axis and the mean tracked frequency and standard deviation.
#'object 2 is a matrix with all the tracked periods.
#'If output = 3, output is a list which contain 2 objects. object 1 is a matrix
#'with the x-axis and the mean tracked frequency and standard deviation.
#'Object 2 is a GIF in which
#'all the tracked periods are plotted.
#'If output = 4, output is a list which contain 3 objects. Object 1 is a matrix
#'with all the tracked periods. Object 2 is a GIF in which  all the tracked
#'periods are plotted.
#'If output = 4 output is a list which contain 3 objects. Object 1 is a matrix
#'with all the tracked periods. Object 2 is a GIF in which  all the tracked
#'periods are plotted.
#'If output = 5 a matrix with the x-axis and the mean tracked period and
#'standard deviation is returned.
#'If output = 6, a matrix with all the tracked periods is returned.
#'If output = 7, a GIF in which all the tracked periods are plotted is returned
#'
#'@examples
#'\donttest{
#'# Re-track the 110kyr eccentricity cycle in the wavelet scalogram
#'# from the XRF record of the Bisciaro data set of Arts (2014)
#'
#'Bisciaro_al <- Bisciaro_XRF[, c(1, 61)]
#'Bisciaro_al <- astrochron::sortNave(Bisciaro_al,verbose=FALSE,genplot=FALSE)
#'Bisciaro_al <- astrochron::linterp(Bisciaro_al, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_al <- Bisciaro_al[Bisciaro_al[, 1] > 2, ]
#'
#'Bisciaro_al_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_al,
#'    dj = 1 /50 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_al_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_al_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'# Bisciaro_al_wt_track <- completed_series(
#'#   wavelet = Bisciaro_al_wt,
#'#   tracked_curve = Bisciaro_al_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_al_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_al_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'Bisciaro_ca <- Bisciaro_XRF[, c(1, 55)]
#'Bisciaro_ca <- astrochron::sortNave(Bisciaro_ca,verbose=FALSE,genplot=FALSE)
#'Bisciaro_ca <- astrochron::linterp(Bisciaro_ca, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_ca <- Bisciaro_ca[Bisciaro_ca[, 1] > 2, ]
#'
#'Bisciaro_ca_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_ca,
#'    dj = 1 /50 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_ca_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_ca_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'# Bisciaro_ca_wt_track <- completed_series(
#'#   wavelet = Bisciaro_ca_wt,
#'#   tracked_curve = Bisciaro_ca_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_ca_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_ca_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE)
#'
#'Bisciaro_sial <- Bisciaro_XRF[,c(1,64)]
#'Bisciaro_sial <- astrochron::sortNave(Bisciaro_sial,verbose=FALSE,genplot=FALSE)
#'Bisciaro_sial <- astrochron::linterp(Bisciaro_sial, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_sial <- Bisciaro_sial[Bisciaro_sial[, 1] > 2, ]
#'
#'Bisciaro_sial_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_sial,
#'    dj = 1 /50 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_sial_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_sial_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'#
#'# Bisciaro_sial_wt_track <- completed_series(
#'#   wavelet = Bisciaro_sial_wt,
#'#   tracked_curve = Bisciaro_sial_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_sial_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_sial_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'
#'Bisciaro_Mn <- Bisciaro_XRF[,c(1,46)]
#'Bisciaro_Mn <- astrochron::sortNave(Bisciaro_Mn,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mn <- astrochron::linterp(Bisciaro_Mn, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mn <- Bisciaro_Mn[Bisciaro_Mn[, 1] > 2, ]
#'
#'Bisciaro_Mn_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_Mn,
#'    dj = 1 /50 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_Mn_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_Mn_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'#
#'# Bisciaro_Mn_wt_track <- completed_series(
#'#   wavelet = Bisciaro_Mn_wt,
#'#   tracked_curve = Bisciaro_Mn_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'# Bisciaro_Mn_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_Mn_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE
#'#   )
#'
#'Bisciaro_Mg <- Bisciaro_XRF[,c(1,71)]
#'Bisciaro_Mg <- astrochron::sortNave(Bisciaro_Mg,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mg <- astrochron::linterp(Bisciaro_Mg, dt = 0.01,verbose=FALSE,genplot=FALSE)
#'Bisciaro_Mg <- Bisciaro_Mg[Bisciaro_Mg[, 1] > 2, ]
#'
#'Bisciaro_Mg_wt <-
#'  analyse_wavelet(
#'    data = Bisciaro_Mg,
#'    dj = 1 /50 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'# Bisciaro_Mg_wt_track <-
#'#   track_period_wavelet(
#'#     astro_cycle = 110,
#'#     wavelet = Bisciaro_Mg_wt,
#'#     n.levels = 100,
#'#     periodlab = "Period (metres)",
#'#     x_lab = "depth (metres)"
#'#   )
#'#
#'#
#'# Bisciaro_Mg_wt_track <- completed_series(
#'#   wavelet = Bisciaro_Mg_wt,
#'#   tracked_curve = Bisciaro_Mg_wt_track,
#'#   period_up = 1.2,
#'#   period_down = 0.8,
#'#   extrapolate = TRUE,
#'#   genplot = FALSE,
#'#   keep_editable = FALSE
#'# )
#'#
#'# Bisciaro_Mg_wt_track <-
#'#   loess_auto(
#'#     time_series = Bisciaro_Mg_wt_track,
#'#     genplot = FALSE,
#'#     print_span = FALSE,
#'#     keep_editable = FALSE)
#'
#'
#'
#'
#'wt_list_bisc <- list(Bisciaro_al_wt,
#'                Bisciaro_ca_wt,
#'                Bisciaro_sial_wt,
#'                Bisciaro_Mn_wt,
#'                Bisciaro_Mg_wt)
#'
#'#Instead of tracking, the tracked solution data sets Bisciaro_al_wt_track,
#'#Bisciaro_ca_wt_track, Bisciaro_sial_wt_track, Bisciaro_Mn_wt_track,
#'# Bisciaro_Mn_wt_track and Bisciaro_Mg_wt_track are used
#'
#'data_track_bisc <- cbind(Bisciaro_al_wt_track[,2],
#'                      Bisciaro_ca_wt_track[,2],
#'                      Bisciaro_sial_wt_track[,2],
#'                      Bisciaro_Mn_wt_track[,2],
#'                      Bisciaro_Mg_wt_track[,2])
#'
#'x_axis_bisc <- Bisciaro_al_wt_track[,1]
#'
#'
#'bisc_retrack <- retrack_wt_MC(wt_list = wt_list_bisc,
#'              data_track = data_track_bisc,
#'              x_axis = x_axis_bisc,
#'              nr_simulations = 20,
#'              seed_nr = 1337,
#'              verbose = FALSE,
#'              genplot = FALSE,
#'              keep_editable = FALSE,
#'              create_GIF = FALSE,
#'              plot_GIF = FALSE,
#'              width_plt =  600,
#'              height_plt = 450,
#'             period_up  =  1.5,
#'              period_down = 0.5,
#'              plot.COI = TRUE,
#'              n.levels = 100,
#'              palette_name = "rainbow",
#'              color_brewer = "grDevices",
#'              periodlab = "Period (metres)",
#'              x_lab = "depth (metres)",
#'              add_avg = FALSE,
#'              time_dir = TRUE,
#'              file_name = NULL,
#'              run_multicore = FALSE,
#'              output = 1,
#'              n_imgs = 50,
#'              plot_horizontal = TRUE,
#'              empty_folder = FALSE)
#'
#'
#'}
#'
#' @param colour_brewer British spelling of \code{color_brewer}, with the same options; if given (not \code{NULL}) it overrides
#' \code{color_brewer}. \code{Default=NULL}.
#' @export
#' @importFrom astrochron linterp
#' @importFrom stats quantile
#' @importFrom graphics par
#' @importFrom graphics image
#' @importFrom graphics axis
#' @importFrom graphics mtext
#' @importFrom graphics text
#' @importFrom graphics box
#' @importFrom graphics polygon
#' @importFrom graphics title
#' @importFrom grDevices rgb
#' @importFrom astrochron mtm
#' @importFrom graphics abline
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom astrochron sortNave
#' @importFrom astrochron noLow
#' @importFrom astrochron mwStats


retrack_wt_MC <- function(wt_list = NULL,
    data_track = NULL,
    x_axis = NULL,
    smoothing = c("auto"),
    nr_simulations = 50,
    seed_nr = 1337,
    verbose = FALSE,
    genplot = FALSE,
    keep_editable = FALSE,
    create_GIF = FALSE,
    plot_GIF = FALSE,
    width_plt =  600,
    height_plt = 450,
    period_up  =  1.5,
    period_down = 0.5,
    plot.COI = TRUE,
    n.levels = 100,
    palette_name = "viridis",
    color_brewer = "viridis",
    periodlab = "Period (metres)",
    x_lab = "depth (metres)",
    add_avg = FALSE,
    time_dir = TRUE,
    file_name = NULL,
    run_multicore = FALSE,
    n_cores = NULL,
    output = 1,
    n_imgs = 50,
    plot_horizontal = TRUE,
    empty_folder = FALSE,
    colour_brewer = NULL){
  if (!is.null(colour_brewer)) {
    color_brewer <- colour_brewer
  }  # British spelling synonym

  simulations = nr_simulations
  img_animated <- NULL





  if (empty_folder==TRUE & create_GIF==TRUE){
    f <- list.files(file_name, include.dirs = F, full.names = T, recursive = T)
    file.remove(f)
  }



  if (create_GIF==TRUE){
    dir.create(file_name, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }

  n_curves <- ncol(data_track)

  set.seed(seed_nr)
  fit <-  wr_lapply(seq_len(simulations), function(j) {

    sel_curve <- sample(1:n_curves, 1, replace=F)
    n <- n_curves
    x <- runif(n, 0, 1)
    y <- x / sum(x)


    vals <- matrix(rep(t(y),nrow(data_track)),ncol=ncol(data_track),byrow=TRUE)
    fractions <- data_track*vals
    fractions <- rowSums(fractions)
    fractions <- cbind(x_axis,fractions)

    wt_sel <- wt_list[[sel_curve]]

    completed_curve <- WaverideR::completed_series(
        wavelet =  wt_sel,
        tracked_curve =  fractions[,c(1,2)],
        period_up  = period_up,
        period_down  = period_down,
        extrapolate = TRUE,
        genplot = FALSE
    )

    # on the positions of the scalogram (linterp() on a grid of step x_axis[2] - x_axis[1] could end just short of
    # the last position by rounding, which left the last row NA)
    completed_curve <- as.matrix(completed_curve)
    completed_curve <- data.frame(x_axis, stats::approx(completed_curve[, 1], completed_curve[, 2], xout = x_axis,
        rule = 2, ties = mean)$y)

    if(smoothing[1]=="auto"){
      completed_curve <- WaverideR::loess_auto(completed_curve)
    }

    if(smoothing[1]=="loess"){
      completed_curve <- astrochron::noLow(completed_curve,smooth=smoothing[2],output=2,genplot=F,verbose=F)

    }
    if(smoothing[1]=="window"){
      completed_curve <- astrochron::mwStats(completed_curve,win=smoothing[2],conv=1,ends=T,CI=0,output=T,genplot=F,verbose=F)
    }







    interp <- approx(completed_curve[, 1],
        completed_curve[, 2], x_axis, method = "linear")
    completed_curve <- cbind(interp$x, interp$y)
    completed_curve <- completed_curve[,c(1,2)]


    if (create_GIF==TRUE){
      png(filename=paste0(file_name,"/",file_name,"_",j,".jpeg"),width=width_plt,height=height_plt)

      WaverideR::plot_wavelet(wavelet = wt_sel,
          plot.COI = plot.COI,
          n.levels = n.levels,
          useRaster = TRUE,
          palette_name = palette_name,
          color_brewer = color_brewer,
          periodlab = periodlab,
          x_lab = x_lab,
          add_lines=completed_curve,
          add_avg= add_avg,
          dev_new = FALSE,
          time_dir = time_dir,
          plot_horizontal = plot_horizontal)
      dev.off()
    }
    completed_curve <- completed_curve[,2]

  }, n_cores = n_cores, run_multicore = run_multicore, combine = "cbind", progress = verbose)






  sims_2 <- 1/fit
  sims_mean_2 <- rowMeans(sims_2,na.rm = TRUE)
  sims_2 <-  as.matrix(sims_2)
  sims_sd_2 <- wr_row_sds(sims_2,na.rm = TRUE)

  sed_run <- cbind(x_axis,sims_mean_2,sims_sd_2)

  if (genplot == TRUE) {
    if (keep_editable == FALSE) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
    }
    layout.matrix <- matrix(c(1), nrow = 1, ncol = 1)
    graphics::layout(mat = layout.matrix,
        heights = c(1),
        # Heights of the two rows
        widths = c(1))
    par(mar = c(4, 4, 1, 1))
    plot(x=sed_run[,1],y=1/sed_run[,2],type="l",ylim=c(min(1/(sed_run[,2]-sed_run[,3])),
        max(1/(sed_run[,2]+sed_run[,3]))),col="green",lwd=2,
        xlab=x_lab,ylab=periodlab)
    lines(x=sed_run[,1],y=1/(sed_run[,2]-sed_run[,3]),col="red",lwd=2)
    lines(x=sed_run[,1],y=1/(sed_run[,2]+sed_run[,3]),col="blue",lwd=2)
  }

  if (create_GIF==TRUE){
    if (!requireNamespace("magick", quietly = TRUE)) {
      stop("create_GIF = TRUE needs the 'magick' package: install.packages('magick')", call. = FALSE)
    }
    imgs <- list.files(file_name, full.names = TRUE)

    if (n_imgs > nr_simulations){
      n_imgs <- nr_simulations
    }

    imgs <- imgs[1:n_imgs]

    img_list <- vector("list", length(imgs))
    for (k in seq_along(imgs)) {
      img_list[[k]] <- magick::image_read(imgs[[k]])
    }
    img_joined <- magick::image_join(img_list)
    img_animated <- magick::image_animate(img_joined, fps = 5)

    if (plot_GIF==TRUE){
      img_animated
    }

    magick::image_write(image = img_animated,
        path =paste0(file_name,"/",file_name,".gif"))
  }

  colnames(sed_run) <- c("depth","mean_freq","sd")

  if (output == 1) {
    res <- list(sed_run,fit,img_animated)

  }

  if (output == 2) {
    res <- list(sed_run,fit)

  }

  if (output == 3) {
    res <- list(sed_run,img_animated)

  }


  if (output == 4) {
    res <- list(fit,img_animated)

  }


  if (output == 5) {
    res <- sed_run

  }

  if (output == 6) {
    res <- fit

  }

  if (output == 7) {
    res <- img_animated

  }


  return(res)
}
