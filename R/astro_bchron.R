#' @title Use astrochronological input in a Bchronology run
#'
#' @description Run a modified version of the Bchronology algorithm from
#' the R package "Bchron and modifiedBChron " (Haslett and Parnell, 2008) and (Trayler et al., 2020)
#'
#'
#' @param age_constraints A data frame object with age constraints
#' The data frame should consist of 6 columns
#' column 1: ID code of the age constraint
#' column 2: Age of the age constraint
#' column 3: Standard deviation in age (should be the SI number as the age)
#' column 4: Position of the age constraints in the depth domain
#' column 5: Thickness of the interval from which the the age constraints
#' originate from
#' column 6: calCurve values assign as "G" for a Gaussian distribution
#' @param n_simulations number of simulation runs
#' @param track_period_incl_sd A matrix of 3 columns in which the first column
#' is depth/height.The second column is the period of the tracked cycle.
#' The thirds column is uncertainty given as 1 standard deviation for the
#' period of the tracked cycle
#'@param tracked_cycle_period period in time of the tracked cycle
#'@param tracked_cycle_period_unc uncertainty in the period of the tracked cycle
#'@param tracked_cycle_period_unc_dist distribution of the uncertainty of the
#'tracked cycle value need to be either "u" for uniform distribution or
#'"n" for normal distribution  \code{Default="u"}
#' @param output If output = 1 a matrix with the predicted ages given the input
#'If output = 2 a matrix with 5 columns is generated, the first column is
#'depth/height, the second column are the quantile
#'(0.025,0.373,0.5,0.6827,0.975) age values of the runs
#'@param run_multicore Run function using multiple cores \code{Default="FALSE"}
#'
#'@references Jennifer Kasbohm, Blair Schoene, Alessandro Montanari,
#'Rodolfo Coccioni,High-precision U-Pb zircon geochronology of the Miocene
#'Bisciaro Formation, Contessa Section, Italy: A case study for requisite
#'radioisotopic calibration of bio- and magnetostratigraphy, Palaeogeography,
#'Palaeoclimatology, Palaeoecology, Volume 576, 2021, 110487, ISSN 0031-0182,
#'<\doi{doi:10.1016/j.palaeo.2021.110487}>
#'
#'
#'@examples
#'\donttest{
#'
#'# U/Pb ages of Kasbohm et al., (2021)
#'id <- c("CCT18_322", "CCT18_315", "CCT18_311")
#' ages <- c(20158, 20575, 20857)
#' ageSds <- c(28, 40, 34)
#' position <- c(10.2, 6.1, 3.2)
#' thickness <- c(0.2, 0.1, 0.1)
#' calCurves <- rep("G", length(ages))
#' ash_Bisc <- as.data.frame(cbind(id, ages, ageSds, position, thickness, calCurves))
#'
#'
#'Bisciaro_al <- Bisciaro_XRF[, c(1, 61)]
#'Bisciaro_al <-
#'  astrochron::sortNave(Bisciaro_al, verbose = FALSE, genplot = FALSE)
#'Bisciaro_al <-
#'  astrochron::linterp(Bisciaro_al,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_al <- Bisciaro_al[Bisciaro_al[, 1] > 2.8,]
#'
#'Bisciaro_al_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_al,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'
#'Bisciaro_ca <- Bisciaro_XRF[, c(1, 55)]
#'Bisciaro_ca <-
#'  astrochron::sortNave(Bisciaro_ca, verbose = FALSE, genplot = FALSE)
#'Bisciaro_ca <-
#'  astrochron::linterp(Bisciaro_ca,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_ca <- Bisciaro_ca[Bisciaro_ca[, 1] > 2.8,]
#'
#'Bisciaro_ca_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_ca,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'Bisciaro_sial <- Bisciaro_XRF[, c(1, 64)]
#'Bisciaro_sial <-
#'  astrochron::sortNave(Bisciaro_sial, verbose = FALSE, genplot = FALSE)
#'Bisciaro_sial <-
#'  astrochron::linterp(Bisciaro_sial,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_sial <- Bisciaro_sial[Bisciaro_sial[, 1] > 2.8,]
#'
#'Bisciaro_sial_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_sial,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'
#'Bisciaro_Mn <- Bisciaro_XRF[, c(1, 46)]
#'Bisciaro_Mn <-
#'  astrochron::sortNave(Bisciaro_Mn, verbose = FALSE, genplot = FALSE)
#'Bisciaro_Mn <-
#'  astrochron::linterp(Bisciaro_Mn,
#'                      dt = 0.01,
#'                      verbose = FALSE,
#'                      genplot = FALSE)
#'Bisciaro_Mn <- Bisciaro_Mn[Bisciaro_Mn[, 1] > 2.8,]
#'
#'Bisciaro_Mn_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_Mn,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'
#'Bisciaro_Mg <- Bisciaro_XRF[, c(1, 71)]
#'Bisciaro_Mg <-
#'  astrochron::sortNave(Bisciaro_Mg, verbose = FALSE, genplot = FALSE)
#'Bisciaro_Mg <-astrochron::linterp(Bisciaro_Mg,
#'                                  dt = 0.01,
#'                                  verbose = FALSE,
#'                                  genplot = FALSE)
#'Bisciaro_Mg <- Bisciaro_Mg[Bisciaro_Mg[, 1] > 2.8,]
#'
#'Bisciaro_Mg_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_Mg,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 8
#'  )
#'
#'Bisciaro_ca_wt <-
#'  analyze_wavelet(
#'    data = Bisciaro_ca,
#'    dj = 1 / 200 ,
#'    lowerPeriod = 0.01,
#'    upperPeriod = 50,
#'    verbose = FALSE,
#'    omega_nr = 6
#'  )
#'
#'wt_list_bisc <- list(Bisciaro_al_wt,
#'                     Bisciaro_ca_wt,
#'                     Bisciaro_sial_wt,
#'                     Bisciaro_Mn_wt,
#'                     Bisciaro_Mg_wt)
#'
#'data_track_bisc <- cbind(
#'  Bisciaro_al_wt_track[, 2],
#'  Bisciaro_ca_wt_track[, 2],
#'  Bisciaro_sial_wt_track[, 2],
#'  Bisciaro_Mn_wt_track[, 2],
#'  Bisciaro_Mg_wt_track[, 2]
#')
#'
#'x_axis_bisc <- Bisciaro_al_wt_track[, 1]
#'
#'retrack <- retrack_wt_MC(
#'  wt_list = wt_list_bisc,
#'  data_track = data_track_bisc,
#'  x_axis = x_axis_bisc,
#'  nr_simulations = 1000,
#'  seed_nr = 1337,
#'  verbose = TRUE,
#'  genplot = TRUE,
#'  keep_editable = FALSE,
#'  create_GIF = FALSE,
#'  plot_GIF = FALSE,
#'  width_plt =  600,
#'  height_plt = 450,
#'  period_up  =  1.5,
#'  period_down = 0.5,
#'  plot.COI = TRUE,
#'  n.levels = 100,
#'  palette_name = "rainbow",
#'  color_brewer = "grDevices",
#'  periodlab = "Period (metres)",
#'  x_lab = "depth (metres)",
#'  add_avg = FALSE,
#'  time_dir = TRUE,
#'  file_name = NULL,
#'  run_multicore = FALSE,
#'  output = 1,
#'  n_imgs = 50,
#'  plot_horizontal = TRUE,
#'  empty_folder = TRUE)
#'
#'retrack <- retrack[[1]]
#'
#'
#'astro_bchron_res <- astro_bchron(age_constraints = ash_Bisc,
#'n_simulations = 10,
#'track_period_incl_sd = retrack,
#'tracked_cycle_period = 110,
#'tracked_cycle_period_unc = 21,
#'tracked_cycle_period_unc_dist = "u",
#'output=1,
#'run_multicore = FALSE)}
#'
#'
#'@return
#'If output = 1 a matrix with the predicted ages given the input for each run
#'If output = 2 a matrix with 6 columns is generated, the first column is
#'depth/height, the other columns are the quantile
#'(0.025,0.373,0.5,0.6827,0.975) age values of the runs
#'
#' @useDynLib WaverideR, .registration = TRUE
#' @export
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom stats dunif
#' @importFrom stats var
#' @importFrom stats quantile
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom tcltk setTkProgressBar
#' @importFrom tcltk setTkProgressBar


astro_bchron <- function(age_constraints = NULL,
                              n_simulations = NULL,
                              track_period_incl_sd = NULL,
                              tracked_cycle_period = NULL,
                         tracked_cycle_period_unc =NULL,
                         tracked_cycle_period_unc_dist = "u",
                              output=1,
                              run_multicore = FALSE) {
  ash <- age_constraints

  multi_tracked <- track_period_incl_sd

  age_curves <- matrix(data = NA,
                       nrow = nrow(multi_tracked),
                       ncol = n_simulations)
  new_curve <- multi_tracked[, c(1, 2)]

  for (i in 1:ncol(age_curves)) {
    val <-rnorm(1, mean = multi_tracked[1, 2], sd = multi_tracked[1, 3])
    pnorm_val <- pnorm(val, mean = multi_tracked[1, 2], sd = multi_tracked[1, 3])
    for (j in 1:nrow(new_curve)) {
      new_curve[j, 2] <- qnorm(pnorm_val, mean = multi_tracked[j, 2], sd = multi_tracked[j, 3])
    }

    if (tracked_cycle_period_unc_dist == "u"){
    tracked_cycle_period_new <- runif(1, min = tracked_cycle_period-tracked_cycle_period_unc, max = tracked_cycle_period+tracked_cycle_period_unc)
    }

    if (tracked_cycle_period_unc_dist == "n"){
    tracked_cycle_period_new <- rnorm(1, mean = tracked_cycle_period, sd = tracked_cycle_period_unc)
    }

    time_curve <- WaverideR::curve2time(
      tracked_cycle_curve = new_curve,
      tracked_cycle_period = tracked_cycle_period_new,
      genplot = FALSE,
      keep_editable = FALSE
    )

    age_curves[, i] <-
      time_curve[, 2]
  }

  ash$postime <- as.numeric(ash$ages) - as.numeric(ash$ages[1])

  retrack_1_time <- curve2time(
    tracked_cycle_curve = multi_tracked[, c(1, 2)],
    tracked_cycle_period = tracked_cycle_period,
    genplot = TRUE,
    keep_editable = FALSE
  )

  ash$pos_time_rand <- NA
  ash_2 <- ash

  pos_thick_rand <-
    matrix(data = NA, nrow(ash_2), ncol = n_simulations)
  pos_time_rand <-
    matrix(data = NA, nrow(ash_2), ncol = n_simulations)
  pos_rownr_rand <-
    matrix(data = NA, nrow(ash_2), ncol = n_simulations)
  pos_rownr_min <-
    matrix(data = NA, nrow(ash_2), ncol = n_simulations)
  pos_rownr_plus <-
    matrix(data = NA, nrow(ash_2), ncol = n_simulations)


  for (i in 1:ncol(pos_time_rand)) {
    for (k in 1:nrow(ash_2)) {
      unc <- as.numeric(ash_2[k, 5])
      position <-  as.numeric(ash_2[k, 4])
      rand_pos_plus <-
        DescTools::Closest(multi_tracked[, 1], position + (unc / 2), which = TRUE)
      rand_pos_min <-
        DescTools::Closest(multi_tracked[, 1], position - (unc / 2), which = TRUE)
      rand_pos_plus <- rand_pos_plus[1]
      rand_pos_min <- rand_pos_min[1]
      pos_rownr_min[k, i] <- rand_pos_min
      pos_rownr_plus[k, i] <- rand_pos_plus
      pos_thick_rand[k, i] <-
        age_curves[rand_pos_plus, i] - age_curves[rand_pos_min, i]
      rown_nr <-
        DescTools::Closest(multi_tracked[, 1], position, which = TRUE)
      rown_nr <- rown_nr[1]
      pos_rownr_rand[k, i] <- rown_nr
      pos_time_rand[k, i] <- age_curves[rown_nr, i]
    }
  }




  #load C functions
  truncatedWalk = function(old, sd, low, high) {
    if (isTRUE(all.equal(low, high, tolerance = 1e-10)))
      return(list(new = low, rat = 1))
    new = .C(
      "truncatedWalk",
      as.double(old),
      as.double(sd),
      as.double(low),
      as.double(high),
      as.double(0),
      PACKAGE = "WaverideR"
    )[5][[1]]
    rat = .C(
      "truncatedRat",
      as.double(old),
      as.double(sd),
      as.double(low),
      as.double(high),
      as.double(new),
      as.double(0),
      PACKAGE = "WaverideR"
    )[6][[1]]
    return(list(new = new, rat = rat))
  }
  dtweediep1 = Vectorize(function(y, p, mu, phi) {
    return(
      .C(
        "dtweediep1",
        as.double(y),
        as.double(p),
        as.double(mu),
        as.double(phi),
        as.double(0),
        PACKAGE = "WaverideR"
      )[5][[1]]
    )
  })

  predictInterp = function(alpha,
                           lambda,
                           beta,
                           predictPositions,
                           diffPositionj,
                           currPositionsj,
                           currPositionsjp1,
                           thetaj,
                           thetajp1) {
    return(
      .C(
        "predictInterp",
        as.double(alpha),
        as.double(lambda),
        as.double(beta),
        as.double(predictPositions),
        as.integer(length(predictPositions)),
        as.double(diffPositionj),
        as.double(currPositionsj),
        as.double(currPositionsjp1),
        as.double(thetaj),
        as.double(thetajp1),
        as.double(rep(0, length(
          predictPositions
        ))),
        PACKAGE = "WaverideR"
      )[11][[1]]
    )
  }
  predictExtrapUp = function(alpha,
                             lambda,
                             beta,
                             predictPositions,
                             currPositions1,
                             theta1,
                             maxExtrap,
                             extractDate) {
    return(
      .C(
        "predictExtrapUp",
        as.double(alpha),
        as.double(lambda),
        as.double(beta),
        as.double(predictPositions),
        as.integer(length(predictPositions)),
        as.double(currPositions1),
        as.double(theta1),
        as.integer(maxExtrap),
        as.double(extractDate),
        as.double(rep(0, length(
          predictPositions
        ))),
        PACKAGE = "WaverideR"
      )[10][[1]]
    )
  }
  predictExtrapDown = function(alpha,
                               lambda,
                               beta,
                               predictPositions,
                               currPositions1,
                               theta1,
                               maxExtrap,
                               extractDate) {
    return(
      .C(
        "predictExtrapDown",
        as.double(alpha),
        as.double(lambda),
        as.double(beta),
        as.double(predictPositions),
        as.integer(length(predictPositions)),
        as.double(currPositions1),
        as.double(theta1),
        as.integer(maxExtrap),
        as.double(extractDate),
        as.double(rep(0, length(
          predictPositions
        ))),
        PACKAGE = "WaverideR"
      )[10][[1]]
    )
  }

  compoundProb <- function(ages, sigs, distType, x) {
    interval <- matrix(0, nrow = length(x), ncol = length(ages))
    for (i in 1:length(ages)) {
      if (distType[i] == "G") {
        interval[, i] <- dnorm(x, ages[i], sigs[i]) / length(ages) *
          mean(diff(x))
      }
      else if (distType[i] == "U") {
        interval[, i] <- dunif(x, ages[i] - sigs[i],
                               ages[i] + sigs[i]) / length(ages) * mean(diff(x))
      }
    }
    G <- apply(interval, 1, sum)
    return(G)
  }



  if (run_multicore == TRUE) {
    j <- 1
    numCores <- detectCores()
    cl <- makeCluster(numCores - 2)
    registerDoSNOW(cl)

    pb <- txtProgressBar(max = n_simulations, style = 3)
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    age_res <- foreach::foreach(
      j = 1:(n_simulations),
      .options.snow = opts,
      .errorhandling = "remove",
      .combine = "cbind"
    ) %dopar% {
      ages = as.numeric(ash_2[, 2])
      ageSds = as.numeric(ash_2[, 3])
      positions = as.numeric(pos_time_rand[, j])
      positionThicknesses = as.numeric(pos_thick_rand[, j])
      ids = ash_2[, 1]
      distTypes = ash[, 6]
      iterations = 5000
      burn = 1000
      probability = 0.95
      predictPositions = age_curves[, j]
      truncateUp = 0
      extrapUp = 1000
      truncateDown = 1e+9
      extrapDown = 1000


      #prepair data
      o = order(positions)
      if (any(positions[o] != positions)) {
        #warning("positions not given in order - re-ordering")
        ages <- ages[o]
        ageSds <- ageSds[o]
        positions <- positions[o]
        distTypes <- distTypes[o]
        positionThicknesses <- positionThicknesses[o]
        ids = ids[o]
      }
      nSamples <- length(unique(ids))
      masterPositions <- vector()
      nNames <- unique(ids)

      for (i in 1:nSamples) {
        masterPositions[i] <- positions[ids == nNames[i]][1]
      }

      nThicknesses <- vector(length = nSamples)

      for (i in 1:nSamples) {
        nThicknesses[i] <- unique(positionThicknesses[positions ==
                                                        masterPositions[i]])
      }
      nAges <- rep(NA, length(nSamples))

      for (i in 1:nSamples) {
        nAges[i] <- length(ids[ids == nNames[i]])
      }
      prob <- matrix(0, nrow = 1e+05, ncol = nSamples)
      ageGrid <- seq(min(ages - ageSds * 10), max(ages + ageSds *
                                                    10), length.out = 1e+05)
      for (j in 1:nSamples) {
        prob[, j] <- compoundProb(ages[ids == nNames[j]], ageSds[ids ==
                                                                   nNames[j]], distType = distTypes[ids == nNames[j]],
                                  x = ageGrid)
      }
      rm(j)
      thetaStore <- matrix(NA, nrow = iterations, ncol = nSamples)
      muStore <- vector(length = iterations)
      psiStore <- muStore
      modelStore <-
        matrix(NA, ncol = iterations, nrow = length(predictPositions))
      positionStore <-
        matrix(NA, nrow = iterations, ncol = nSamples)
      predictStore <-
        matrix(NA, ncol = iterations, nrow = length(predictPositions))
      psiSDStore <- psiStore
      muSDStore <- psiStore
      mhSDStore <- thetaStore
      thetas <- vector(length = nSamples)


      for (i in 1:nSamples) {
        thetas[i] <- ageGrid[which.max(prob[, i])] + rnorm(1,
                                                           0, 1e-05)
      }
      currPositions <- masterPositions
      mhSD <- runif(length(unique(ids)))
      psiSD <- runif(1)
      muSD <- runif(1)
      mu <-
        abs(rnorm(1, mean = mean(diff(thetas)) / mean(diff(
          currPositions
        )),
        sd = muSD))
      psi <- abs(rnorm(1, 1, sd = psiSD))
      p = 1.2
      alpha <- (2 - p) / (p - 1)

      #run Bchronology algorithm
      for (n in 1:iterations) {
        lambda <- (mu ^ (2 - p)) / (psi * (2 - p))
        beta <- 1 / (psi * (p - 1) * mu ^ (p - 1))
        currPositions <-
          runif(nSamples,
                masterPositions - nThicknesses,
                masterPositions + nThicknesses)
        do <- order(currPositions)
        diffPositions <- diff(currPositions[do])
        thetas[do] <- sort(thetas, decreasing = T)
        positionStore[n, ] <- currPositions
        for (j in 1:(nSamples - 1)) {
          inRange <- predictPositions >= currPositions[do[j]] &
            predictPositions <= currPositions[do[j + 1]]
          predval <- predictInterp(
            alpha = alpha,
            lambda = lambda,
            beta = beta,
            predictPositions = predictPositions[inRange],
            diffPositionj = diffPositions[j],
            currPositionsj = currPositions[do[j]],
            currPositionsjp1 = currPositions[do[j + 1]],
            thetaj = thetas[do[j]],
            thetajp1 = thetas[do[j +
                                   1]]
          )
          modelStore[inRange, n] <- predval
        }
        if (any(predictPositions >= currPositions[do[nSamples]])) {
          inRange <- predictPositions >= currPositions[do[nSamples]]
          predval <- predictExtrapUp(
            alpha = alpha,
            lambda = lambda,
            beta = beta,
            predictPositions = predictPositions[inRange],
            theta1 = thetas[do[nSamples]],
            currPositions1 = currPositions[do[nSamples]],
            maxExtrap = extrapUp,
            extractDate = truncateUp
          )
          modelStore[inRange, n] <- predval
        }
        if (any(predictPositions <= currPositions[do[1]])) {
          inRange <- predictPositions <= currPositions[do[1]]
          predval <- predictExtrapDown(
            alpha = alpha,
            lambda = lambda,
            beta = beta,
            predictPositions = predictPositions[inRange],
            theta1 = thetas[do[1]],
            currPositions1 = currPositions[do[1]],
            maxExtrap = extrapDown,
            extractDate = truncateDown
          )
          modelStore[inRange, n] <- predval
        }

        for (i in 1:nSamples) {
          thetaCurrent <- thetas[do[i]]
          thetaProposed <- truncatedWalk(
            old = thetaCurrent,
            sd = mhSD[do[i]],
            low = ifelse(i == nSamples,
                         truncateUp, thetas[do[i + 1]] - 1e-10),
            high = ifelse(i ==
                            1, 1e+10, thetas[do[i - 1]] + 1e-10)
          )
          pProposed <- log(prob[, do[i]][which.min(abs(ageGrid -
                                                         thetaProposed$new))])
          pProposed <- max(pProposed, -1e+06, na.rm = T)
          pCurrent <- log(prob[, do[i]][which.min(abs(ageGrid -
                                                        thetaCurrent))])
          pCurrent <- max(pCurrent, -1e+06, na.rm = T)
          priorProposed <-
            ifelse(i == 1, 0, log(
              dtweediep1(
                thetas[do[i -
                            1]] - thetaProposed$new,
                p = p,
                mu = mu * (diffPositions[i -
                                           1]),
                phi = psi * (diffPositions[i - 1]) ^ (p -
                                                        1)
              )
            )) + ifelse(i == nSamples, 0, log(
              dtweediep1(
                thetaProposed$new -
                  thetas[do[i + 1]],
                p = p,
                mu = mu * (diffPositions[i]),
                phi = psi * (diffPositions[i]) ^
                  (p - 1)
              )
            ))
          priorCurrent <-
            ifelse(i == 1, 0, log(
              dtweediep1(
                thetas[do[i -
                            1]] - thetaCurrent,
                p = p,
                mu = mu * (diffPositions[i -
                                           1]),
                phi = psi * (diffPositions[i - 1]) ^ (p -
                                                        1)
              )
            )) + ifelse(i == nSamples, 0, log(
              dtweediep1(
                thetaCurrent -
                  thetas[do[i + 1]],
                p = p,
                mu = mu * (diffPositions[i]),
                phi = psi * (diffPositions[i]) ^
                  (p - 1)
              )
            ))
          priorCurrent <- max(priorCurrent, -1e+06, na.rm = T)
          priorProposed <- max(priorProposed, -1e+06, na.rm = T)
          logRtheta <- priorProposed - priorCurrent + pProposed -
            pCurrent + log(thetaProposed$rat)
          logRtheta <- max(logRtheta, -1e+06, na.rm = T)
          if (runif(1, 0, 1) < min(1, exp(logRtheta))) {
            thetas[do[i]] <- thetaProposed$new
          }
        }
        thetaStore[n, ] <- thetas
        muCurrent <- mu
        muProposed <- truncatedWalk(
          old = mu,
          sd = muSD,
          low = 1e-10,
          high = 1e+10
        )
        priorMu <- sum(log(
          dtweediep1(
            abs(diff(thetas[do])),
            p = p,
            mu = muProposed$new * diffPositions,
            phi = psi / diffPositions ^ (p -
                                           1)
          )
        )) - sum(log(
          dtweediep1(
            abs(diff(thetas[do])),
            p = p,
            mu = mu * diffPositions,
            phi = psi / diffPositions ^ (p -
                                           1)
          )
        )) + log(muProposed$rat)
        mu <- ifelse(runif(1, 0, 1) < min(1, exp(priorMu)),
                     muProposed$new, muCurrent)
        muStore[n] <- mu
        psiCurrent <- psi
        psiProposed <- truncatedWalk(
          old = psi,
          sd = psiSD,
          low = 1e-10,
          high = 1e+10
        )
        priorPsi <- sum(log(
          dtweediep1(
            abs(diff(thetas[do])),
            p = p,
            mu = mu * diffPositions,
            phi = psiProposed$new / diffPositions ^ (p -
                                                       1)
          )
        )) - sum(log(dtweediep1(
          abs(diff(thetas[do])),
          p = p,
          mu = mu * diffPositions,
          phi = psi / (diffPositions ^ (p -
                                          1))
        ))) + log(psiProposed$rat)
        psi <- ifelse(runif(1, 0, 1) < min(1, exp(priorPsi)),
                      psiProposed$new,
                      psiCurrent)
        psiStore[n] <- psi
        h = 200
        if (n %% h == 0) {
          cd <- 2.4 / sqrt(1)
          psiK <- psiStore[(n - h + 1):n] - mean(psiStore[(n -
                                                             h + 1):n])
          psiRt <- var(psiK)
          psiSD <- sqrt((cd ^ 2) * psiRt)
          psiSD <- ifelse(is.na(psiSD), cd * sd(psiStore,
                                                na.rm = T), psiSD)
          psiSD <- ifelse(psiSD == 0, cd * sd(psiStore, na.rm = T),
                          psiSD)
          muK <- muStore[(n - h + 1):n] - mean(muStore[(n -
                                                          h + 1):n])
          muRt <- var(muK)
          muSD <- sqrt(cd ^ 2 * muRt)
          muSD <- ifelse(is.na(muSD), cd * sd(muStore, na.rm = T),
                         muSD)
          muSD <- ifelse(muSD == 0, cd * sd(muStore, na.rm = T),
                         muSD)
          for (q in 1:nSamples) {
            thetaK <- thetaStore[(n - h + 1):n, q] - mean(thetaStore[(n -
                                                                        h + 1):n, q])
            thetaRt <- var(thetaK)
            mhSD[q] <- ifelse(is.na(sqrt(cd ^ 2 * thetaRt)),
                              cd * sd(thetaStore[, q], na.rm = T),
                              sqrt(cd ^ 2 *
                                     thetaRt))
          }
        }
        psiSDStore[n] <- psiSD
        muSDStore[n] <- muSD
        mhSDStore[n, ] <- mhSD
      }


      HDI <-
        as.vector(cbind(t(apply(
          modelStore[, burn:iterations], 1, quantile, c(0.5)
        ))))
    }

  } else{
    res_ages <- matrix(data = NA, nrow(age_curves), ncol = 0)

    for (j in 1:n_simulations) {
      ages = as.numeric(ash_2[, 2])
      ageSds = as.numeric(ash_2[, 3])
      positions = as.numeric(pos_time_rand[, j])
      positionThicknesses = as.numeric(pos_thick_rand[, j])
      ids = ash_2[, 1]
      distTypes = ash[, 6]
      iterations = 5000
      burn = 1000
      probability = 0.95
      predictPositions = age_curves[, j]
      truncateUp = 0
      extrapUp = 1000
      truncateDown = 1e+9
      extrapDown = 1000


      #prepair data
      o = order(positions)
      if (any(positions[o] != positions)) {
        #warning("positions not given in order - re-ordering")
        ages <- ages[o]
        ageSds <- ageSds[o]
        positions <- positions[o]
        distTypes <- distTypes[o]
        positionThicknesses <- positionThicknesses[o]
        ids = ids[o]
      }
      nSamples <- length(unique(ids))
      masterPositions <- vector()
      nNames <- unique(ids)

      for (i in 1:nSamples) {
        masterPositions[i] <- positions[ids == nNames[i]][1]
      }

      nThicknesses <- vector(length = nSamples)

      for (i in 1:nSamples) {
        nThicknesses[i] <- unique(positionThicknesses[positions ==
                                                        masterPositions[i]])
      }
      nAges <- rep(NA, length(nSamples))

      for (i in 1:nSamples) {
        nAges[i] <- length(ids[ids == nNames[i]])
      }
      prob <- matrix(0, nrow = 1e+05, ncol = nSamples)
      ageGrid <- seq(min(ages - ageSds * 10), max(ages + ageSds *
                                                    10), length.out = 1e+05)
      for (j in 1:nSamples) {
        prob[, j] <- compoundProb(ages[ids == nNames[j]], ageSds[ids ==
                                                                   nNames[j]], distType = distTypes[ids == nNames[j]],
                                  x = ageGrid)
      }
      rm(j)
      thetaStore <- matrix(NA, nrow = iterations, ncol = nSamples)
      muStore <- vector(length = iterations)
      psiStore <- muStore
      modelStore <-
        matrix(NA, ncol = iterations, nrow = length(predictPositions))
      positionStore <-
        matrix(NA, nrow = iterations, ncol = nSamples)
      predictStore <-
        matrix(NA, ncol = iterations, nrow = length(predictPositions))
      psiSDStore <- psiStore
      muSDStore <- psiStore
      mhSDStore <- thetaStore
      thetas <- vector(length = nSamples)


      for (i in 1:nSamples) {
        thetas[i] <- ageGrid[which.max(prob[, i])] + rnorm(1,
                                                           0, 1e-05)
      }
      currPositions <- masterPositions
      mhSD <- runif(length(unique(ids)))
      psiSD <- runif(1)
      muSD <- runif(1)
      mu <-
        abs(rnorm(1, mean = mean(diff(thetas)) / mean(diff(
          currPositions
        )),
        sd = muSD))
      psi <- abs(rnorm(1, 1, sd = psiSD))
      p = 1.2
      alpha <- (2 - p) / (p - 1)

      #run Bchronology algorithm
      for (n in 1:iterations) {
        lambda <- (mu ^ (2 - p)) / (psi * (2 - p))
        beta <- 1 / (psi * (p - 1) * mu ^ (p - 1))
        currPositions <-
          runif(nSamples,
                masterPositions - nThicknesses,
                masterPositions + nThicknesses)
        do <- order(currPositions)
        diffPositions <- diff(currPositions[do])
        thetas[do] <- sort(thetas, decreasing = T)
        positionStore[n, ] <- currPositions
        for (j in 1:(nSamples - 1)) {
          inRange <- predictPositions >= currPositions[do[j]] &
            predictPositions <= currPositions[do[j + 1]]
          predval <- predictInterp(
            alpha = alpha,
            lambda = lambda,
            beta = beta,
            predictPositions = predictPositions[inRange],
            diffPositionj = diffPositions[j],
            currPositionsj = currPositions[do[j]],
            currPositionsjp1 = currPositions[do[j + 1]],
            thetaj = thetas[do[j]],
            thetajp1 = thetas[do[j +
                                   1]]
          )
          modelStore[inRange, n] <- predval
        }
        if (any(predictPositions >= currPositions[do[nSamples]])) {
          inRange <- predictPositions >= currPositions[do[nSamples]]
          predval <- predictExtrapUp(
            alpha = alpha,
            lambda = lambda,
            beta = beta,
            predictPositions = predictPositions[inRange],
            theta1 = thetas[do[nSamples]],
            currPositions1 = currPositions[do[nSamples]],
            maxExtrap = extrapUp,
            extractDate = truncateUp
          )
          modelStore[inRange, n] <- predval
        }
        if (any(predictPositions <= currPositions[do[1]])) {
          inRange <- predictPositions <= currPositions[do[1]]
          predval <- predictExtrapDown(
            alpha = alpha,
            lambda = lambda,
            beta = beta,
            predictPositions = predictPositions[inRange],
            theta1 = thetas[do[1]],
            currPositions1 = currPositions[do[1]],
            maxExtrap = extrapDown,
            extractDate = truncateDown
          )
          modelStore[inRange, n] <- predval
        }

        for (i in 1:nSamples) {
          thetaCurrent <- thetas[do[i]]
          thetaProposed <- truncatedWalk(
            old = thetaCurrent,
            sd = mhSD[do[i]],
            low = ifelse(i == nSamples,
                         truncateUp, thetas[do[i + 1]] - 1e-10),
            high = ifelse(i ==
                            1, 1e+10, thetas[do[i - 1]] + 1e-10)
          )
          pProposed <- log(prob[, do[i]][which.min(abs(ageGrid -
                                                         thetaProposed$new))])
          pProposed <- max(pProposed, -1e+06, na.rm = T)
          pCurrent <- log(prob[, do[i]][which.min(abs(ageGrid -
                                                        thetaCurrent))])
          pCurrent <- max(pCurrent, -1e+06, na.rm = T)
          priorProposed <-
            ifelse(i == 1, 0, log(
              dtweediep1(
                thetas[do[i -
                            1]] - thetaProposed$new,
                p = p,
                mu = mu * (diffPositions[i -
                                           1]),
                phi = psi * (diffPositions[i - 1]) ^ (p -
                                                        1)
              )
            )) + ifelse(i == nSamples, 0, log(
              dtweediep1(
                thetaProposed$new -
                  thetas[do[i + 1]],
                p = p,
                mu = mu * (diffPositions[i]),
                phi = psi * (diffPositions[i]) ^
                  (p - 1)
              )
            ))
          priorCurrent <-
            ifelse(i == 1, 0, log(
              dtweediep1(
                thetas[do[i -
                            1]] - thetaCurrent,
                p = p,
                mu = mu * (diffPositions[i -
                                           1]),
                phi = psi * (diffPositions[i - 1]) ^ (p -
                                                        1)
              )
            )) + ifelse(i == nSamples, 0, log(
              dtweediep1(
                thetaCurrent -
                  thetas[do[i + 1]],
                p = p,
                mu = mu * (diffPositions[i]),
                phi = psi * (diffPositions[i]) ^
                  (p - 1)
              )
            ))
          priorCurrent <- max(priorCurrent, -1e+06, na.rm = T)
          priorProposed <- max(priorProposed, -1e+06, na.rm = T)
          logRtheta <- priorProposed - priorCurrent + pProposed -
            pCurrent + log(thetaProposed$rat)
          logRtheta <- max(logRtheta, -1e+06, na.rm = T)
          if (runif(1, 0, 1) < min(1, exp(logRtheta))) {
            thetas[do[i]] <- thetaProposed$new
          }
        }
        thetaStore[n, ] <- thetas
        muCurrent <- mu
        muProposed <- truncatedWalk(
          old = mu,
          sd = muSD,
          low = 1e-10,
          high = 1e+10
        )
        priorMu <- sum(log(
          dtweediep1(
            abs(diff(thetas[do])),
            p = p,
            mu = muProposed$new * diffPositions,
            phi = psi / diffPositions ^ (p -
                                           1)
          )
        )) - sum(log(
          dtweediep1(
            abs(diff(thetas[do])),
            p = p,
            mu = mu * diffPositions,
            phi = psi / diffPositions ^ (p -
                                           1)
          )
        )) + log(muProposed$rat)
        mu <- ifelse(runif(1, 0, 1) < min(1, exp(priorMu)),
                     muProposed$new, muCurrent)
        muStore[n] <- mu
        psiCurrent <- psi
        psiProposed <- truncatedWalk(
          old = psi,
          sd = psiSD,
          low = 1e-10,
          high = 1e+10
        )
        priorPsi <- sum(log(
          dtweediep1(
            abs(diff(thetas[do])),
            p = p,
            mu = mu * diffPositions,
            phi = psiProposed$new / diffPositions ^ (p -
                                                       1)
          )
        )) - sum(log(dtweediep1(
          abs(diff(thetas[do])),
          p = p,
          mu = mu * diffPositions,
          phi = psi / (diffPositions ^ (p -
                                          1))
        ))) + log(psiProposed$rat)
        psi <- ifelse(runif(1, 0, 1) < min(1, exp(priorPsi)),
                      psiProposed$new,
                      psiCurrent)
        psiStore[n] <- psi
        h = 200
        if (n %% h == 0) {
          cd <- 2.4 / sqrt(1)
          psiK <- psiStore[(n - h + 1):n] - mean(psiStore[(n -
                                                             h + 1):n])
          psiRt <- var(psiK)
          psiSD <- sqrt((cd ^ 2) * psiRt)
          psiSD <- ifelse(is.na(psiSD), cd * sd(psiStore,
                                                na.rm = T), psiSD)
          psiSD <- ifelse(psiSD == 0, cd * sd(psiStore, na.rm = T),
                          psiSD)
          muK <- muStore[(n - h + 1):n] - mean(muStore[(n -
                                                          h + 1):n])
          muRt <- var(muK)
          muSD <- sqrt(cd ^ 2 * muRt)
          muSD <- ifelse(is.na(muSD), cd * sd(muStore, na.rm = T),
                         muSD)
          muSD <- ifelse(muSD == 0, cd * sd(muStore, na.rm = T),
                         muSD)
          for (q in 1:nSamples) {
            thetaK <- thetaStore[(n - h + 1):n, q] - mean(thetaStore[(n -
                                                                        h + 1):n, q])
            thetaRt <- var(thetaK)
            mhSD[q] <- ifelse(is.na(sqrt(cd ^ 2 * thetaRt)),
                              cd * sd(thetaStore[, q], na.rm = T),
                              sqrt(cd ^ 2 *
                                     thetaRt))
          }
        }
        psiSDStore[n] <- psiSD
        muSDStore[n] <- muSD
        mhSDStore[n, ] <- mhSD
      }


      HDI <-
        as.vector(cbind(t(apply(
          modelStore[, burn:iterations], 1, quantile, c(0.5)
        ))))

      res_ages <- cbind(res_ages, HDI)

    }
  }

  if (output == 1) {
    res <- res_ages

  }

  if (output == 2) {
    res <- run_multicore[,1]
    res <- cbind(res, quantile(res_ages,probs=c(0.025,0.373,0.5,0.6827,0.975),na.rm=TRUE))
  }



  return(res)

}
