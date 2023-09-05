#'@title Re-track cycles using a Monte-Carlo simulation
#'
#'@description When analyzing multi-proxy records an age-model can be created
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
#'be off in certain intervals. the link[WaverideR]{retrack_wt_MC} mitigates
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
#'@param wt_list a list containing all the wavelet objects created using the
#'link[WaverideR]{analyze_wavelet} wavelet function
#'To create a list use the link[base]{list} function
#'@param data_track a matrix containing all the tracked period values.
#'To create the matrix use the link[base]{cbind} function and only add the
#'tracked period values so do not add the depth axis. When combining the
#'tracked period values make sure that all curves have a similar depth
#'spacing.
#'@param x_axis The x-axis of the tracked period values
#'@param nr_simulations The number of Monte-Carlo simulations which are to be
#' conducted\code{Default=50}
#'@param seed_nr The seed number of the Monte-Carlo simulations.
#' \code{Default=1337}
#'@param verbose Print text when running the function  \code{Default=FALSE}
#'@param genplot Plot a plot with the mean period and + and - standard
#'deviation \code{Default=FALSE}
#'@param keep_editable Keep option to add extra features after plotting
#' \code{Default=FALSE}
#'@param create_GIF Create a GIF with the re-tracked lines in the wavelet
#' scalograms \code{Default=FALSE}
#'@param plot_GIF Plot a GIF with the re-tracked lines in the wavelet
#' scalograms\code{Default=FALSE}
#'@param width_plt width of the re-tracked plot \code{Default=600}
#'@param height_plt width of the re-tracked plot \code{Default=450}
#'@param period_up The period_up parameter is the factor with which the linear interpolated tracked_curve
#'curve is multiplied by. This linear interpolated tracked_curve multiplied by the period_up factor is
#'the upper boundary which is used  for detecting the spectral peak nearest to the linear interpolated tracked_curve
#'curve. If no spectral peak is detected within the specified boundary the interpolated
#'value is used instead. between spectral peaks \code{Default=1.5},
#'@param period_down  The period_down parameter is the factor with which the linear interpolated tracked_curve
#'curve is multiplied by. This linear interpolated tracked_curve multiplied by the period_down factor is
#'the lower boundary which is used  for detecting the spectral peak nearest to the linear interpolated tracked_curve
#'curve. If no spectral peak is detected within the specified boundary the interpolated
#'value is used instead. between spectral peaks \code{Default=0.5},
#'@param plot.COI Option to plot the cone of influence \code{Default=TRUE}.
#'@param n.levels  Number of color levels \code{Default=100}.
#'@param palette_name Name of the color palette which is used for plotting.
#'The color palettes than can be chosen depends on which the R package is specified in
#'the color_brewer parameter. The included R packages from which palettes can be chosen
#'from are; the 'RColorBrewer', 'grDevices', 'ColorRamps' and 'Viridis' R packages.
#'There are many options to choose from so please
#'read the documentation of these packages \code{Default=rainbow}.
#'The R package 'viridis' has the color palette options: “magma”, “plasma”,
#'“inferno”, “viridis”, “mako”, and “rocket”  and “turbo”
#'To see the color palette options of the The R pacakge 'RColorBrewer' run
#'the RColorBrewer::brewer.pal.info() function
#'The R package 'colorRamps' has the color palette options:"blue2green",
#'"blue2green2red", "blue2red",	"blue2yellow", "colorRamps",	"cyan2yellow",
#'"green2red", "magenta2green", "matlab.like", "matlab.like2" and	"ygobb"
#'The R package 'grDevices' has the built in  palette options:"rainbow",
#'"heat.colors", "terrain.colors","topo.colors" and "cm.colors"
#'To see even more color palette options of the The R pacakge 'grDevices' run
#'the grDevices::hcl.pals() function
#'@param color_brewer Name of the R package from which the color palette is chosen from.
#'The included R packages from which palettes can be chosen
#'are; the RColorBrewer, grDevices, ColorRamps and Viridis R packages.
#'There are many options to choose from so please
#'read the documentation of these packages. "\code{Default=grDevices}
#'@param periodlab Label for the y-axis \code{Default="Period (metres)"}.
#'@param x_lab Label for the x-axis \code{Default="depth (metres)"}.
#'@param add_avg Plot the average wavelet spectral power to the side of the wavelet \code{Default=FALSE}
#'@param time_dir The direction of the proxy record which is assumed for tuning if time increases with increasing depth/time values
#'(e.g. bore hole data which gets older with increasing depth ) then time_dir should be set to TRUE
#'if time decreases with depth/time values (eg stratospheric logs where 0m is the bottom of the section)
#'then time_dir should be set to FALSE \code{time_dir=TRUE}
#'@param file_name Name of the images created using this function. Each file
#'gets a number added to it which corresponds to which number of simulation it was
#'the files are saved in a folder with a similar name created in the current
#'directory
#'@param run_multicore Run function using multiple cores \code{Default="FALSE"}
#'@param output #'If output = 1, output is a list which contain 3 objects.
#'object 1 is a matrix with the x-axis and the mean tracked frequency and
#'standard deviation. #'object 2 is a matrix with all the tracked periods.
#'Object 3 is a GIF in which #'all the tracked periods are plotted.
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
#'If output = 5 a matrix with the x-axis and the mean tracked frequency and
#'standard deviation is returned.
#'If output = 6, a matrix with all the tracked periods is returned.
#'If output = 7, a GIF in which all the tracked periods are plotted is returned.
#'\code{Default=1}
#'@param n_imgs Number images used in creating the GIF a high number of images
#'is computationally intensive and will create a large sized GIF \code{Default=50}
#'@param plot_horizontal plot the wavelet horizontal or vertical eg y axis
#' is depth or y axis power \code{Default=TRUE}
#'@param empty_folder Empty the folder in which the images created using
#'this function are saved \code{Default=FALSE}
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
#'  analyze_wavelet(
#'    data = Bisciaro_al,
#'    dj = 1 /200 ,
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
#'  analyze_wavelet(
#'    data = Bisciaro_ca,
#'    dj = 1 /200 ,
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
#'  analyze_wavelet(
#'    data = Bisciaro_sial,
#'    dj = 1 /200 ,
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
#'  analyze_wavelet(
#'    data = Bisciaro_Mn,
#'    dj = 1 /200 ,
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
#'  analyze_wavelet(
#'    data = Bisciaro_Mg,
#'    dj = 1 /200 ,
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
#' @export
#' @importFrom magick image_write
#' @importFrom magick image_read
#' @importFrom magick image_join
#' @importFrom magick image_animate
#' @importFrom rlist list.extract
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
#' @importFrom WaveletComp analyze.wavelet
#' @importFrom WaveletComp wt.image
#' @importFrom biwavelet wt
#' @importFrom astrochron mtm
#' @importFrom DescTools Closest
#' @importFrom graphics abline
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom colorRamps blue2green
#' @importFrom colorRamps blue2green2red
#' @importFrom colorRamps blue2red
#' @importFrom colorRamps blue2yellow
#' @importFrom colorRamps cyan2yellow
#' @importFrom colorRamps green2red
#' @importFrom colorRamps magenta2green
#' @importFrom colorRamps matlab.like
#' @importFrom colorRamps matlab.like2
#' @importFrom colorRamps ygobb
#' @importFrom viridis viridis
#' @importFrom viridis magma
#' @importFrom viridis plasma
#' @importFrom viridis inferno
#' @importFrom viridis cividis
#' @importFrom viridis mako
#' @importFrom viridis rocket
#' @importFrom viridis turbo
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom colorRamps blue2green
#' @importFrom colorRamps blue2green2red
#' @importFrom colorRamps blue2red
#' @importFrom colorRamps blue2yellow
#' @importFrom colorRamps cyan2yellow
#' @importFrom colorRamps green2red
#' @importFrom colorRamps magenta2green
#' @importFrom colorRamps matlab.like
#' @importFrom colorRamps matlab.like2
#' @importFrom colorRamps ygobb
#' @importFrom viridis viridis
#' @importFrom viridis magma
#' @importFrom viridis plasma
#' @importFrom viridis inferno
#' @importFrom viridis cividis
#' @importFrom viridis mako
#' @importFrom viridis rocket
#' @importFrom viridis turbo
#' @importFrom grDevices rainbow
#' @importFrom grDevices heat.colors
#' @importFrom grDevices terrain.colors
#' @importFrom grDevices topo.colors
#' @importFrom grDevices cm.colors
#' @importFrom grDevices hcl.colors
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom astrochron sortNave
#' @importFrom parallel stopCluster

retrack_wt_MC <- function(wt_list = NULL,
                          data_track = NULL,
                          x_axis = NULL,
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
                          palette_name = "rainbow",
                          color_brewer = "grDevices",
                          periodlab = "Period (metres)",
                          x_lab = "depth (metres)",
                          add_avg = FALSE,
                          time_dir = TRUE,
                          file_name = NULL,
                          run_multicore = FALSE,
                          output = 1,
                          n_imgs = 50,
                          plot_horizontal = TRUE,
                          empty_folder = FALSE){

  simulations = nr_simulations
  img_animated <- NULL

  if (run_multicore == TRUE) {
    numCores <- detectCores()
    cl <- parallel::makeCluster(numCores - 2)
    registerDoSNOW(cl)
  }

  if (empty_folder==TRUE & create_GIF==TRUE){
    f <- list.files(file_name, include.dirs = F, full.names = T, recursive = T)
    file.remove(f)}


  if (verbose==TRUE){
    pb <- txtProgressBar(max = simulations, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)}else{opts=NULL}

  if (create_GIF==TRUE){
    dir.create(file_name, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }

  n_curves <- ncol(data_track)



  if (run_multicore == TRUE) {
    j <- 1 # needed to assign 1 to ijk to avoid note
    set.seed(seed_nr)
    fit <-  foreach (j = 1:simulations, .options.snow = opts, .combine = 'cbind') %dopar% {

      sel_curve <- sample(1:n_curves, 1, replace=F)
      n <- n_curves
      x <- runif(n, 0, 1)
      y <- x / sum(x)


      vals <- matrix(rep(t(y),nrow(data_track)),ncol=ncol(data_track),byrow=TRUE)
      fractions <- data_track*vals
      fractions <- rowSums(fractions)
      fractions <- cbind(x_axis,fractions)

      wt_sel <- rlist::list.extract(wt_list, sel_curve)

      completed_curve <- WaverideR::completed_series(
        wavelet =  wt_sel,
        tracked_curve =  fractions[,c(1,2)],
        period_up  = period_up,
        period_down  = period_down,
        extrapolate = TRUE,
        genplot = FALSE
      )

      completed_curve <- astrochron::linterp(completed_curve,dt=x_axis[2]-x_axis[1],start=x_axis[1],genplot = FALSE,verbose=FALSE)
      completed_curve <- completed_curve[1:length(x_axis),]
      completed_curve <- WaverideR::loess_auto(completed_curve)
      completed_curve <- completed_curve[,c(1,2)]


      if (create_GIF==TRUE){
        png(filename=paste0(file_name,"/",file_name,"_",j,".jpeg"),type="cairo",width=width_plt,height=height_plt)

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

    }

    stopCluster(cl)
  }

  if (run_multicore == FALSE) {

  fit <- matrix(data=NA,nrow=length(x_axis),ncol=simulations)

  for (j in 1:simulations){
  sel_curve <- sample(1:n_curves, 1, replace=F)
  n <- n_curves
  x <- runif(n, 0, 1)
  y <- x / sum(x)


  vals <- matrix(rep(t(y),nrow(data_track)),ncol=ncol(data_track),byrow=TRUE)
  fractions <- data_track*vals
  fractions <- rowSums(fractions)
  fractions <- cbind(x_axis,fractions)

  wt_sel <- rlist::list.extract(wt_list, sel_curve)

  completed_curve <- WaverideR::completed_series(
    wavelet =  wt_sel,
    tracked_curve =  fractions[,c(1,2)],
    period_up  = period_up,
    period_down  = period_down,
    extrapolate = TRUE,
    genplot = FALSE
  )

  completed_curve <- astrochron::linterp(completed_curve,dt=x_axis[2]-x_axis[1],start=x_axis[1],genplot = FALSE,verbose=FALSE)
  completed_curve <- completed_curve[1:length(x_axis),]
  completed_curve <- WaverideR::loess_auto(completed_curve)
  completed_curve <- completed_curve[,c(1,2)]


  if (create_GIF==TRUE){
    png(filename=paste0(file_name,"/",file_name,"_",j,".jpeg"),type="cairo",width=width_plt,height=height_plt)

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
  fit[,j] <- completed_curve[,2]

}
}


  sims_2 <- 1/fit
  sims_mean_2 <- rowMeans(sims_2)
  sims_2 <-  as.matrix(sims_2)
  sims_sd_2 <- matrixStats::rowSds(sims_2)

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
    imgs <- list.files(file_name, full.names = TRUE)

    if (n_imgs > nr_simulations){
      n_imgs <- nr_simulations
    }

    imgs <- imgs[1:n_imgs]

    img_list <- lapply(imgs, image_read)
    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 5)

    if (plot_GIF==TRUE){
      img_animated
    }

    image_write(image = img_animated,
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


  return(res)}
