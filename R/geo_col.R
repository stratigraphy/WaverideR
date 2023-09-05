#' @title Generate standard color codes for the Geological Time Scale
#'
#' @description Generates the R color code which corresponds its respective
#' geological subdivision
#'
#' @param name Name of the geologchronological subdivision
#'
#'
#'@references
#'Ogg, Gabi & Ogg, James & Gradstein, Felix. (2021).
#'Recommended color coding of stages - Appendix 1
#'from Geologic Time Scale 2020.
#'
#' @examples
#'#generate the Silurian part of the GTS
#'plot.new()
#'plot(
#'  x = c(0, 1),
#'  y = c(419.2, 443.8),
#'  col = "white",
#'  xlab = "",
#'  ylab = "Time (Ma)",
#'  xaxt = "n",
#'  xaxs = "i",
#'  yaxs = "i",
#'  ylim = rev(c(419, 444))
#')            # Draw empty plot
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Rhuddanian"),
#'  col =geo_col("Rhuddanian")
#')
#'
#'text(
#'  0.85,geo_mid("Rhuddanian"),
#'  "Rhuddanian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Aeronian"),
#'  col =geo_col("Aeronian")
#')
#'
#'text(
#'  0.85,geo_mid("Aeronian"),
#'  "Aeronian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Telychian"),
#'  col =geo_col("Telychian")
#')
#'
#'text(
#'  0.85,geo_mid("Telychian"),
#'  "Telychian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Sheinwoodian"),
#'  col =geo_col("Sheinwoodian")
#')
#'
#'text(
#'  0.85,geo_mid("Sheinwoodian"),
#'  "Sheinwoodian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Homerian"),
#'  col =geo_col("Homerian")
#')
#'
#'text(
#'  0.85,geo_mid("Homerian"),
#'  "Homerian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Gorstian"),
#'  col =geo_col("Gorstian")
#')
#'
#'text(
#'  0.85,geo_mid("Gorstian"),
#'  "Gorstian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Ludfordian"),
#'  col =geo_col("Ludfordian")
#')
#'
#'text(
#'  0.85,geo_mid("Ludfordian"),
#'  "Ludfordian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0.66, 1, 1, 0.66),
#'  y = geo_loc("Pridoli_Age"),
#'  col =geo_col("Pridoli_Age")
#')
#'
#'
#'
#'polygon(
#'  x = c(0.33, 0.66, 0.66, 0.33),
#'  y = geo_loc("Pridoli"),
#'  col =geo_col("Pridoli")
#')
#'
#'text(
#'  0.5,geo_mid("Pridoli"),
#'  "Pridoli",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'
#'polygon(
#'  x = c(0.33, 0.66, 0.66, 0.33),
#'  y = geo_loc("Ludlow"),
#'  col =geo_col("Ludlow")
#')
#'
#'text(
#'  0.5,geo_mid("Ludlow"),
#'  "Ludlow",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0.33, 0.66, 0.66, 0.33),
#'  y = geo_loc("Wenlock"),
#'  col =geo_col("Wenlock")
#')
#'
#'text(
#'  0.5,geo_mid("Wenlock"),
#'  "Wenlock",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0.33, 0.66, 0.66, 0.33),
#'  y = geo_loc("Llandovery"),
#'  col =geo_col("Llandovery")
#')
#'
#'text(
#'  0.5,geo_mid("Llandovery"),
#'  "Llandovery",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'polygon(
#'  x = c(0, 0.33, 0.33, 0),
#'  y = geo_loc("Silurian"),
#'  col =geo_col("Silurian")
#')
#'
#'text(
#'  0.165,geo_mid("Silurian"),
#'  "Silurian",
#'  cex = 1,
#'  col = "black",
#'  srt = 0
#')
#'
#'
#'@return
#'Returns the color code of the geological subdivision
#' @export

geo_col <- function(name=NULL){
  GTS_info <- WaverideR::GTS_info
  cols_geo_sel <- rgb(as.data.frame(c(GTS_info[GTS_info[,1]==name,11:13],255))/255)
  return(cols_geo_sel)
}
