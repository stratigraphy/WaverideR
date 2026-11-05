# Internal definitions of the mother wavelets used by analyse_wavelet() and
# analyse_wavelet_sst(). All wavelets are defined in the frequency domain
# (Torrence and Compo 1998): the daughter wavelet at scale s is
# sqrt(2 * pi * s / dt) * psi_hat(s * omega), with psi_hat normalised so that
# the integral of |psi_hat|^2 over omega is 1.

.wr_wavelet_names <- c("morlet", "paul", "dog", "morse", "bump", "cmor", "cgau", "shannon",
    "fbsp", "meyer", "loggabor", "morlet0", "morl", "hhat", "poisson", "gammatone",
    "hann", "haar", "db", "sym", "battle", "custom")

# real (non-analytic) wavelets: the transform is real, gives no phase/amplitude envelope and cannot
# be synchrosqueezed
.wr_real_wavelets <- c("dog", "morl", "hhat", "poisson", "haar", "db", "sym", "battle")

# Alternative names: name -> c(canonical name, order or NA to use the default order)
.wr_wavelet_aliases <- list(
    gabor = list("morlet", NULL), mexh = list("dog", 2), mexican_hat = list("dog", 2),
    ricker = list("dog", 2), marr = list("dog", 2), mhat = list("dog", 2),
    gmw = list("morse", NULL), cauchy = list("paul", NULL), lognormal = list("loggabor", NULL),
    complex_morlet = list("cmor", NULL), complex_gaussian = list("cgau", NULL),
    spline = list("fbsp", NULL), sinc = list("shannon", NULL), `log-gabor` = list("loggabor", NULL),
    amor = list("morlet", NULL), morlex = list("morlet0", NULL), exact_morlet = list("morlet0", NULL),
    morlet_exact = list("morlet0", NULL), real_morlet = list("morl", NULL), hermitian_hat = list("hhat", NULL),
    daubechies = list("db", NULL), symlet = list("sym", NULL), symlets = list("sym", NULL),
    battle_lemarie = list("battle", NULL), lemarie = list("battle", NULL), bl = list("battle", NULL),
    raised_cosine = list("hann", NULL), hanning = list("hann", NULL), cmhat = list("cgau", 2),
    complex_mexican_hat = list("cgau", 2), sombrero = list("dog", 2), gaus = list("dog", NULL),
    gaussian = list("dog", NULL), airy = list("morse", NULL), shan = list("shannon", NULL),
    morl0 = list("morlet0", NULL), gammachirp = list("gammatone", NULL)
)

# Resolve a wavelet name (aliases, PyWavelets style names with the order appended such as "db4",
# "sym8", "cgau3", "gaus2", "cmor1.5-1", "shan2-1.5", "fbsp2-1-1.5") or a user supplied function
# psi_hat(omega) into list(name, order)
.wr_resolve <- function(wavelet, order = NULL){
  if (is.function(wavelet)) {
    return(list(name = "custom", order = wavelet))
  }
  if (!is.character(wavelet) || length(wavelet) != 1) {
    stop("wavelet must be a name or a function psi_hat(omega)")
  }
  name <- gsub("[ ]+", "_", tolower(wavelet))
  if (!(name %in% c(.wr_wavelet_names, names(.wr_wavelet_aliases)))) {
    m <- regmatches(name, regexec("^([a-z_]+?)([0-9][0-9.]*(-[0-9.]+)*)$", name))[[1]]
    if (length(m) && nchar(m[3]) > 0) {
      name <- m[2]
      if (is.null(order)) {
        order <- as.numeric(strsplit(m[3], "-")[[1]])
      }
    }
  }
  if (name %in% names(.wr_wavelet_aliases)){
    al <- .wr_wavelet_aliases[[name]]
    name <- al[[1]]
    if (is.null(order) && !is.null(al[[2]])) {
      order <- al[[2]]
    }
  }
  if (!(name %in% .wr_wavelet_names)) {
    stop("unknown wavelet '", wavelet, "'; available: ", paste(.wr_wavelet_names[-length(.wr_wavelet_names)], collapse = ", "),
        " (or a function psi_hat(omega))")
  }
  if (name == "custom" && !is.function(order)) {
    stop("for wavelet = \"custom\" give wavelet_order = function(omega)")
  }
  list(name = name, order = order)
}

# cardinal B-spline of order m (support [0, m])
.wr_bspline <- function(x, m){
  out <- numeric(length(x))
  for (k in 0:m){
    out <- out + (-1)^k * choose(m, k) * pmax(x - k, 0)^(m - 1)
  }
  out / factorial(m - 1)
}

# psi_hat(omega) of the unit-scale wavelet for a vector of angular frequencies
.wr_psi_hat <- function(name, order, w){
  switch(name,
      morlet = {
        w0 <- order[1]
        ifelse(w > 0, pi^(-1/4) * exp(-(w - w0)^2 / 2), 0)
      },
      paul = {
        m <- order[1]
        ifelse(w > 0, exp(m * log(2) - 0.5 * (log(m) + lfactorial(2 * m - 1)) + m * log(pmax(w, 1e-300)) - w), 0)
      },
      dog = {
        m <- order[1]
        # (-1)^(m + 1) i^m / sqrt(Gamma(m + 1/2)) w^m exp(-w^2/2); the modulus is used
        # for even m (real wavelet) and the phase factor i^m is kept as complex otherwise
        lg <- -0.5 * lgamma(m + 1/2) + m * log(pmax(abs(w), 1e-300)) - w^2 / 2
        sgn <- ifelse(w >= 0, 1, cos(pi * m))
        exp(1i * pi * (m + 1)) * (1i)^m * sgn * exp(lg)
      },
      morse = {
        gam <- order[1]
        beta <- order[2]
        wp <- (beta / gam)^(1 / gam)
        loga <- log(2) - (beta / gam) * log(beta / gam) + beta / gam
        ifelse(w > 0, exp(loga + beta * log(pmax(w, 1e-300)) - w^gam), 0)
      },
      bump = {
        mu <- order[1]
        sg <- order[2]
        u <- (w - mu) / sg
        ifelse(abs(u) < 1, exp(1 - 1 / (1 - u^2)), 0)
      },
      cmor = {
        # complex Morlet with bandwidth fb and centre frequency fc (in cycles), as in Matlab/PyWavelets
        fb <- order[1]
        fc <- order[2]
        ifelse(w > 0, exp(-fb * (w - 2 * pi * fc)^2 / 4), 0)
      },
      cgau = {
        # analytic part of the m-th derivative of a Gaussian
        m <- order[1]
        ifelse(w > 0, exp(m * log(pmax(w, 1e-300)) - w^2 / 2 - (m / 2) * (log(m) - 1)), 0)
      },
      shannon = {
        # complex Shannon (sinc) wavelet: box in the frequency domain
        fb <- order[1]
        fc <- order[2]
        ifelse(abs(w - 2 * pi * fc) <= pi * fb, 1, 0)
      },
      fbsp = {
        # frequency B-spline wavelet of order m, bandwidth fb and centre frequency fc
        m <- order[1]
        fb <- order[2]
        fc <- order[3]
        a <- fb / m
        .wr_bspline((w - 2 * pi * fc) / (2 * pi * a) + m / 2, m)
      },
      meyer = {
        # Meyer wavelet, positive frequencies only
        nu <- function(x) {
          x <- pmin(pmax(x, 0), 1)
          x^4 * (35 - 84 * x + 70 * x^2 - 20 * x^3)
        }
        a <- abs(w)
        out <- ifelse(a >= 2 * pi / 3 & a < 4 * pi / 3, sin(pi / 2 * nu(3 * a / (2 * pi) - 1)),
            ifelse(a >= 4 * pi / 3 & a <= 8 * pi / 3, cos(pi / 2 * nu(3 * a / (4 * pi) - 1)), 0))
        ifelse(w > 0, out, 0)
      },
      loggabor = {
        # log-Gabor: Gaussian on a logarithmic frequency axis, centre frequency w0, log-width sg
        w0 <- order[1]
        sg <- order[2]
        ifelse(w > 0, exp(-(log(pmax(w, 1e-300) / w0))^2 / (2 * sg^2)), 0)
      },
      morlet0 = {
        # Morlet wavelet with the correction term that makes it exactly zero mean (admissible also for a
        # small omega_0)
        w0 <- order[1]
        ifelse(w > 0, pi^(-1/4) * (exp(-(w - w0)^2 / 2) - exp(-(w^2 + w0^2) / 2)), 0)
      },
      morl = {
        # real Morlet wavelet exp(-t^2/2) cos(omega_0 t) (PyWavelets "morl", omega_0 = 5)
        w0 <- order[1]
        exp(-(w - w0)^2 / 2) + exp(-(w + w0)^2 / 2)
      },
      hhat = {
        # Hermitian hat (1 - t^2 + i t) exp(-t^2/2) (Szu et al. 1992)
        (w^2 + w) * exp(-w^2 / 2)
      },
      poisson = {
        # Poisson wavelet of order n: |omega|^n exp(-|omega|); n = 1 is (1 - t^2) / (1 + t^2)^2
        n <- order[1]
        exp(n * log(pmax(abs(w), 1e-300)) - abs(w))
      },
      gammatone = {
        # analytic gammatone (auditory filter) wavelet of order n, centre frequency w0 and bandwidth b;
        # asymmetric in time (sharp onset, slow decay) like the Paul wavelet
        n <- order[1]
        w0 <- order[2]
        b <- order[3]
        # the factor 1 - exp(-(4 w / w0)^2) makes it exactly zero mean (admissible)
        ifelse(w > 0, (1 + 1i * (w - w0) / b)^(-n) * (1 - exp(-(4 * w / w0)^2)), 0)
      },
      hann = {
        # raised cosine (Hann) band in the frequency domain, centre w0 and half width d
        w0 <- order[1]
        d <- order[2]
        ifelse(abs(w - w0) < d & w > 0, cos(pi * (w - w0) / (2 * d))^2, 0)
      },
      haar = .wr_table_psi("db", 1, w),
      db = .wr_table_psi("db", order[1], w),
      sym = .wr_table_psi("sym", order[1], w),
      battle = .wr_table_psi("battle", order[1], w),
      custom = order(w),
      stop("unknown wavelet")
  )
}

# Peak angular frequency of the unit-scale wavelet
.wr_peak_omega <- function(name, order){
  switch(name,
      morlet = order[1],
      paul = (2 * order[1] + 1) / 2,
      dog = sqrt(order[1] + 1/2),
      morse = (order[2] / order[1])^(1 / order[1]),
      bump = order[1],
      cmor = 2 * pi * order[2],
      cgau = sqrt(order[1]),
      shannon = 2 * pi * order[2],
      fbsp = 2 * pi * order[3],
      meyer = {
        w <- seq(0.5, 12, length.out = 4000)
        w[which.max(.wr_psi_hat("meyer", order, w))]
      },
      loggabor = order[1] * exp(-order[2]^2),
      morlet0 = ,
      morl = {
        w <- seq(0.01, order[1] + 6, length.out = 6000)
        w[which.max(Mod(.wr_psi_hat(name, order, w)))]
      },
      gammatone = order[2],
      hann = order[1],
      {
        # numerical peak of |psi_hat| over the positive frequencies
        w <- if (name %in% c("haar", "db", "sym", "battle")) { .wr_table(name, order)$w } else { seq(1e-4, 100, length.out = 20001) }
        w <- w[w > 0]
        w[which.max(Mod(.wr_psi_hat(name, order, w)))]
      })
}

.wr_default_order <- function(name){
  switch(name, morlet = 8, paul = 4, dog = 2, morse = c(3, 20), bump = c(5, 0.6),
      cmor = c(1.5, 1), cgau = 4, shannon = c(2, 1.5), fbsp = c(2, 2, 1.5), meyer = NA_real_,
      loggabor = c(6, 0.25), morlet0 = 6, morl = 5, hhat = NA_real_, poisson = 2,
      gammatone = c(4, 6, 1.5), hann = c(6, 3), haar = 1, db = 4, sym = 8, battle = 3, custom = NULL)
}

# Full definition of a wavelet: normalisation, Fourier factor, reconstruction constant psi0
# and the factor that scales the cone of influence relative to the Morlet wavelet
.wr_wavelet <- function(wavelet = "morlet", order = NULL, omega_nr = 8, ngrid = 40001){
  rs <- .wr_resolve(wavelet, order)
  name <- rs$name
  order <- rs$order
  if (is.null(order)){
    order <- if (name %in% c("morlet", "morlet0")) { omega_nr } else { .wr_default_order(name) }
  }
  if (name == "gammatone" && length(order) == 1) {
    order <- c(order, 6, 1.5)
  }
  if (name == "gammatone" && (length(order) != 3 || any(order <= 0))) {
    stop("for the gammatone wavelet give wavelet_order = c(n, centre angular frequency, bandwidth)")
  }
  if (name == "hann" && length(order) == 1) {
    order <- c(6, order)
  }
  if (name == "hann" && (length(order) != 2 || order[2] <= 0 || order[2] >= order[1])) {
    stop("for the hann wavelet give wavelet_order = c(centre angular frequency w0, half width d) with 0 < d < w0")
  }
  if (name == "poisson" && order[1] < 1) {
    stop("the order of the Poisson wavelet must be >= 1")
  }
  if (name == "db" && !(order[1] %in% 1:20)) {
    stop("Daubechies wavelets are available for orders 1 to 20")
  }
  if (name == "sym" && !(order[1] %in% 2:20)) {
    stop("Symlets are available for orders 2 to 20")
  }
  if (name == "battle" && !(order[1] %in% c(1, 3, 5, 7))) {
    stop("Battle-Lemarie wavelets are available for orders 1, 3, 5 and 7")
  }
  if (name == "cmor" && length(order) == 1) {
    order <- c(order, 1)
  }
  if (name == "shannon" && length(order) == 1) {
    order <- c(order, 1.5)
  }
  if (name == "fbsp" && length(order) == 1) {
    order <- c(order, 2, 1.5)
  }
  if (name == "fbsp" && length(order) == 2) {
    order <- c(order[1], order[2], 1.5)
  }
  if (name == "loggabor" && length(order) == 1) {
    order <- c(6, order)
  }
  if (name == "bump" && length(order) == 1) {
    order <- c(order, 0.6)
  }
  if (name == "bump" && (length(order) != 2 || any(is.na(order)) || order[2] <= 0 || order[1] <= order[2])){
    stop("for the bump wavelet give wavelet_order = c(mu, sigma) with mu > sigma > 0 (default c(5, 0.6))")
  }
  if (name == "cmor" && (length(order) != 2 || any(order <= 0))) {
    stop("for the cmor wavelet give wavelet_order = c(bandwidth, centre frequency)")
  }
  if (name == "cgau" && order[1] < 1) {
    stop("the order of the cgau wavelet must be >= 1")
  }
  if (name == "shannon" && (length(order) != 2 || any(order <= 0))) {
    stop("for the shannon wavelet give wavelet_order = c(bandwidth, centre frequency) with fc > fb/2")
  }
  if (name == "shannon" && order[2] <= order[1] / 2) {
    stop("for the shannon wavelet the centre frequency must exceed half the bandwidth")
  }
  if (name == "fbsp" && (length(order) != 3 || any(order <= 0) || order[1] != round(order[1]))) {
    stop("for the fbsp wavelet give wavelet_order = c(m, bandwidth, centre frequency), m an integer")
  }
  if (name == "fbsp" && order[3] <= order[2] / 2) {
    stop("for the fbsp wavelet the centre frequency must exceed half the bandwidth")
  }
  if (name == "loggabor" && (length(order) != 2 || any(order <= 0))) {
    stop("for the loggabor wavelet give wavelet_order = c(centre angular frequency, log-width)")
  }
  if (name == "morse" && length(order) == 1){
    order <- c(3, order)
  }
  if (name == "paul" && order[1] < 1){
    stop("the order of the Paul wavelet must be >= 1")
  }
  if (name == "dog" && order[1] < 1){
    stop("the order of the DOG wavelet must be >= 1")
  }
  if (name == "morse" && (length(order) != 2 || any(order <= 0))){
    stop("for the Morse wavelet give wavelet_order = c(gamma, beta) or beta alone (gamma = 3)")
  }
  wp <- .wr_peak_omega(name, order)
  # numerical normalisation on a fine frequency grid
  wmax <- if (name == "morlet") { max(order[1] * 3, 30) } else { max(wp * 6, 60) }
  if (name %in% c("haar", "db", "sym", "battle")) {
    wmax <- max(.wr_table(name, order)$w)
  }
  w <- seq(-wmax, wmax, length.out = ngrid)
  ph <- .wr_psi_hat(name, order, w)
  dw <- w[2] - w[1]
  nrm <- sqrt(sum(Mod(ph)^2) * dw)
  if (name == "morlet") {
    nrm <- 1
  }
  psi0 <- Re(sum(ph) * dw / sqrt(2 * pi)) / nrm
  if (name == "dog") {
    psi0 <- Mod(sum(ph) * dw / sqrt(2 * pi)) / nrm
  }
  analytic <- !(name %in% .wr_real_wavelets)
  if (name == "custom") {
    analytic <- sum(Mod(ph[w < 0])^2) < 1e-3 * sum(Mod(ph)^2)
  }
  legacy <- name %in% c("morlet", "paul", "dog", "morse", "bump", "cmor", "cgau", "shannon", "fbsp", "meyer", "loggabor")
  if (!legacy) {
    # reconstruction constant of the new wavelets: scaled so that the reconstruction of a sinusoid
    # (sum over scales of Re(W) / sqrt(s)) behaves as for the Morlet wavelet with omega_0 = 6
    rec <- function(nm, od, nr, an) {
      u <- seq(1e-4, max(wmax, 60), length.out = 40001)
      (if (an) { 1 } else { 2 }) * Mod(sum(.wr_psi_hat(nm, od, u) / u) * (u[2] - u[1])) / nr
    }
    psi0 <- pi^(-1/4) * rec(name, order, nrm, analytic) / rec("morlet", 6, 1, TRUE)
  }
  list(name = name, order = order, analytic = analytic, norm = nrm,
      peak_omega = wp, fourier_factor = 2 * pi / wp,
      # shape normalised to one at the peak and restricted to positive frequencies: used for the
      # amplitude calibrated transforms (superlets)
      shape = function(scale, omega_k){
        v <- .wr_psi_hat(name, order, scale * omega_k)
        pk <- Mod(.wr_psi_hat(name, order, wp))
        v * (omega_k > 0) / pk
      },
      psi0 = if (name == "morlet") { pi^(-1/4) } else { abs(psi0) },
      real = !analytic,
      daughter = function(scale, omega_k, dt){
        sqrt(2 * pi * scale / dt) * .wr_psi_hat(name, order, scale * omega_k) / nrm
      })
}

# Decay time of the wavelet at scale 1: the time at which the energy of the wavelet envelope
# beyond it has dropped to exp(-2) of the energy of one half. It is a measure of how far from
# the edge of the record a wavelet is affected by it. Only ratios to the Morlet wavelet are used
# (see .wr_coi_factor), which keeps the cone of influence of the Morlet wavelet unchanged.
.wr_decay_time <- function(def){
  wmax <- max(def$peak_omega * 6, 60)
  w <- seq(-wmax, wmax, length.out = 8001)
  dw <- w[2] - w[1]
  ph <- .wr_psi_hat(def$name, def$order, w) / def$norm
  t_all <- seq(0, 15, by = 0.02)
  e <- numeric(length(t_all))
  for (ti in seq_along(t_all)) {
    e[ti] <- Mod(sum(ph * exp(1i * w * t_all[ti])) * dw / sqrt(2 * pi))^2
  }
  tail_e <- rev(cumsum(rev(e)))
  r <- tail_e / tail_e[1]
  t_all[which(r <= exp(-2))[1]]
}

# Factor of the cone of influence: Fourier factor times coi_factor times the distance to the
# edge gives the period below which the wavelet is influenced by the edge. For the Morlet
# wavelet this is sqrt(2), as before.
.wr_coi_factor <- function(def){
  if (def$name == "morlet"){
    return(sqrt(2))
  }
  ref <- .wr_wavelet("morlet", omega_nr = def$peak_omega)
  sqrt(2) * .wr_decay_time(ref) / .wr_decay_time(def)
}

# Value of the mother wavelet at time zero used in the reconstruction of a cycle from the
# wavelet transform (Torrence and Compo 1998, eq. 11). Objects made by older versions have no
# psi0 element and use the Morlet value.
.wr_psi0 <- function(w){
  if (is.null(w$psi0)){
    pi^(-1/4)
  } else {
    w$psi0
  }
}

# Resolve the name of a wavelet (aliases) without building it
.wr_canonical_name <- function(wavelet){
  .wr_resolve(wavelet)$name
}

# Quality factor of a wavelet: peak frequency divided by the root mean square width of the
# energy spectrum |psi_hat|^2 over the positive frequencies. For a Morlet wavelet it equals
# sqrt(2) * omega_0, which is how "number of cycles" of other wavelets is expressed in a Morlet
# equivalent (used by the adaptive wavelet and the superlets).
.wr_quality <- function(name, order){
  wp <- .wr_peak_omega(name, order)
  w <- seq(1e-6, max(wp * 6, 60), length.out = 12001)
  e <- Mod(.wr_psi_hat(name, order, w))^2
  e <- e / sum(e)
  m1 <- sum(e * w)
  wp / sqrt(sum(e * (w - m1)^2))
}

# Order of a wavelet that has the same relative bandwidth as a Morlet wavelet with the given
# omega_0. Families without a free bandwidth parameter (Meyer) return their fixed order.
.wr_order_for_omega <- function(name, omega, base = NULL){
  if (is.null(base)) {
    base <- .wr_default_order(name)
  }
  if (name %in% c("morlet", "morlet0", "morl")) {
    return(omega)
  }
  if (name %in% c("meyer", "hhat", "haar", "db", "sym", "battle", "custom")) {
    return(base)
  }
  omega_q <- omega * sqrt(2)      # quality factor of a Morlet wavelet with this omega_0
  clamp_root <- function(f, lo, hi){
    flo <- f(lo)
    fhi <- f(hi)
    if (!is.finite(flo) || !is.finite(fhi) || sign(flo) == sign(fhi)){
      return(if (abs(flo) < abs(fhi)) { lo } else { hi })
    }
    stats::uniroot(f, c(lo, hi), tol = 1e-6)$root
  }
  switch(name,
      paul = clamp_root(function(m) .wr_quality("paul", m) - omega_q, 1, 400),
      dog = clamp_root(function(m) .wr_quality("dog", m) - omega_q, 1, 400),
      cgau = clamp_root(function(m) .wr_quality("cgau", m) - omega_q, 1, 400),
      morse = c(base[1], clamp_root(function(b) .wr_quality("morse", c(base[1], b)) - omega_q, 0.5, 2000)),
      bump = c(base[1], clamp_root(function(sg) .wr_quality("bump", c(base[1], sg)) - omega_q, 0.02, base[1] * 0.98)),
      loggabor = c(base[1], clamp_root(function(sg) .wr_quality("loggabor", c(base[1], sg)) - omega_q, 0.01, 3)),
      cmor = c(base[1], omega_q / (2 * pi * sqrt(base[1]))),
      poisson = clamp_root(function(m) .wr_quality("poisson", m) - omega_q, 1, 400),
      gammatone = c(base[1], base[2], clamp_root(function(b) .wr_quality("gammatone", c(base[1], base[2], b)) - omega_q, 0.01, 20)),
      hann = c(base[1], clamp_root(function(d) .wr_quality("hann", c(base[1], d)) - omega_q, 0.01, base[1] * 0.99)),
      shannon = {
        q1 <- .wr_quality("shannon", c(base[1], 1.5))
        fc <- max(1.5 * omega_q / q1, base[1] / 2 * 1.0001)
        c(base[1], fc)
      },
      fbsp = {
        q1 <- .wr_quality("fbsp", base)
        fc <- max(base[3] * omega_q / q1, base[2] / 2 * 1.0001)
        c(base[1], base[2], fc)
      },
      stop("unknown wavelet"))
}

# Wavelet definition with the bandwidth of a Morlet wavelet of the given omega_0 (cycles-based
# interface of the adaptive wavelet and the superlets)
.wr_wavelet_equiv <- function(wavelet, omega, wavelet_order = NULL, ngrid = 40001){
  rs <- .wr_resolve(wavelet, wavelet_order)
  name <- rs$name
  if (!is.null(rs$order)) {
    wavelet_order <- rs$order
  }
  if (!is.null(wavelet_order) || name %in% c("morlet", "morlet0", "morl")){
    return(.wr_wavelet(wavelet, wavelet_order, omega_nr = omega, ngrid = ngrid))
  }
  ord <- .wr_order_for_omega(name, omega, .wr_wavelet(name, ngrid = 2001)$order)
  .wr_wavelet(name, ord, ngrid = ngrid)
}

## Orthogonal wavelets (Haar, Daubechies, Symlets, Battle-Lemarie) in the frequency domain
.wr_cache <- new.env(parent = emptyenv())

# Low-pass filter (normalised to sum 1) of the Daubechies (extremal phase) or Symlet (least
# asymmetric) wavelet with N vanishing moments, by spectral factorisation of the Daubechies
# polynomial (Daubechies 1992, ch. 6)
.wr_daub_filter <- function(N, type = c("db", "sym")){
  type <- match.arg(type)
  if (N == 1) {
    return(c(0.5, 0.5))
  }
  # z^(N-1) P(y(z)), y = (2 - z - 1/z) / 4 = (-1 + 2 z - z^2) / (4 z)
  pc <- numeric(2 * N - 1)
  base <- c(-1, 2, -1) / 4
  yk <- 1
  for (k in 0:(N - 1)) {
    term <- choose(N - 1 + k, k) * yk               # coefficients of (-1 + 2z - z^2)^k / 4^k
    idx <- (N - 1 - k) + seq_along(term)            # times z^(N-1-k)
    pc[idx] <- pc[idx] + term
    yk <- stats::convolve(yk, rev(base), type = "open")
  }
  r <- polyroot(pc)
  inside <- r[Mod(r) < 1]
  # group conjugate pairs; each group can be replaced by its reciprocal
  grp <- list()
  used <- rep(FALSE, length(inside))
  for (i in seq_along(inside)) {
    if (used[i]) {
      next
    }
    if (abs(Im(inside[i])) < 1e-9) {
      grp[[length(grp) + 1]] <- i
      used[i] <- TRUE
      next
    }
    j <- which(!used & seq_along(inside) != i & Mod(inside - Conj(inside[i])) < 1e-6)[1]
    grp[[length(grp) + 1]] <- c(i, j)
    used[c(i, j)] <- TRUE
  }
  make <- function(sel) {
    roots <- inside
    for (g in seq_along(grp)) {
      if (sel[g]) {
        roots[grp[[g]]] <- 1 / roots[grp[[g]]]
      }
    }
    a <- 1
    for (k in seq_len(N)) {
      a <- stats::convolve(a, rev(c(0.5, 0.5)), type = "open")
    }
    for (rt in roots) {
      a <- c(0, a) - rt * c(a, 0)
    }    # multiply by (z - root)
    a <- Re(a)
    a / sum(a)
  }
  if (type == "db") {
    return(make(rep(FALSE, length(grp))))
  }
  # Symlet: the root selection with the most linear phase
  best <- NULL
  bestv <- Inf
  wq <- seq(0.01, pi - 0.01, length.out = 200)
  ng <- length(grp)
  for (code in 0:(2^ng - 1)) {
    sel <- as.logical(bitwAnd(code, 2^(0:(ng - 1))))
    a <- make(sel)
    m0 <- colSums(a * exp(-1i * outer(seq_along(a) - 1, wq)))
    ph <- .wr_unwrap(Arg(m0))
    v <- sum(stats::resid(stats::lm(ph ~ wq))^2)
    if (v < bestv) {
      bestv <- v
      best <- a
    }
  }
  best
}

.wr_unwrap <- function(p){
  d <- diff(p)
  d <- d - 2 * pi * round(d / (2 * pi))
  c(p[1], p[1] + cumsum(d))
}

# psi_hat on a fine grid of positive angular frequencies (cached), centred in time
.wr_table <- function(type, order){
  if (type == "haar") {
    type <- "db"
    order <- 1
  }
  key <- paste(type, order[1], sep = "_")
  if (!is.null(.wr_cache[[key]])) {
    return(.wr_cache[[key]])
  }
  w <- seq(0, 2 * pi * 48, length.out = 2^15 + 1)
  if (type == "battle") {
    m <- order[1]
    A <- function(om) {
      k <- -60:60
      x <- outer(om, 2 * pi * k, "+") / 2
      v <- ifelse(abs(x) < 1e-12, 1, sin(x) / x)
      rowSums(v^(2 * m + 2))
    }
    sincp <- function(om) ifelse(abs(om) < 1e-12, 1, sin(om / 2) / (om / 2))
    phi <- function(om) sincp(om)^(m + 1) / sqrt(A(om))
    m0 <- function(om) cos(om / 2)^(m + 1) * sqrt(A(om) / A(2 * om))
    psi <- exp(-1i * w / 2) * Conj(m0(w / 2 + pi)) * phi(w / 2)
  } else {
    a <- .wr_daub_filter(order[1], if (type == "sym") { "sym" } else { "db" })
    kk <- seq_along(a) - 1
    m0 <- function(om) as.vector(exp(-1i * outer(om, kk)) %*% a)
    phi <- rep(1 + 0i, length(w))
    for (j in 1:30) {
      phi <- phi * m0(w / 2^(j + 1))
    }
    # phi_hat(w / 2) = prod_{j >= 1} m0(w / 2^(j + 1))
    psi <- exp(-1i * w / 2) * Conj(m0(w / 2 + pi)) * phi
  }
  # remove the linear phase (time shift) so that the wavelet is centred at t = 0
  e <- Mod(psi)^2
  ph <- .wr_unwrap(Arg(psi))
  ok <- e > max(e) * 1e-3
  c0 <- stats::coef(stats::lm(ph[ok] ~ w[ok], weights = e[ok]))[2]
  psi <- psi * exp(-1i * c0 * w)
  out <- list(w = w, psi = psi)
  .wr_cache[[key]] <- out
  out
}

# Interpolate a tabulated real wavelet: psi_hat(-w) = Conj(psi_hat(w))
.wr_table_psi <- function(type, order, w){
  tb <- .wr_table(type, order)
  aw <- abs(w)
  re <- stats::approx(tb$w, Re(tb$psi), aw, rule = 1)$y
  im <- stats::approx(tb$w, Im(tb$psi), aw, rule = 1)$y
  re[is.na(re)] <- 0
  im[is.na(im)] <- 0
  ifelse(w >= 0, complex(real = re, imaginary = im), complex(real = re, imaginary = -im))
}
