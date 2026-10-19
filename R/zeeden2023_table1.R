#' Devonian reference values of Zeeden et al. (2023), EPSL 621, 118348, Table 1
#' @export
zeeden2023_table1 <- function(){
  data.frame(
    source = c("combined Devonian (CTest)", "La2004 model", "Farhat et al. 2022 model",
               "Berger et al. 1992", "La2004 recent (0 Ma)"),
    age_ma = c(NA, 385, 385, 385, 0),
    p_kyr = c(20.70, 20.29, 21.2, 21.98, 25.68),
    p_sd_kyr = c(0.22, NA, NA, NA, NA),
    `k+g5_kyr` = c(19.35, 19.02, 19.82, 20.5, 23.68),
    `k+g2_kyr` = c(18.47, 18.17, 18.89, 19.51, 22.37),
    `k+g4_kyr` = c(16.13, 15.95, 16.39, 16.86, 18.95),
    `k+g3_kyr` = c(16.18, 15.85, 16.51, 16.98, 19.10),
    `k+s3_kyr` = c(29.50, 27.78, 30.64, 32.31, 40.98),
    check.names = FALSE)
}
