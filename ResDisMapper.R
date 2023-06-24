# Landscape genetic analysis with ResDisMapper
rm(list = ls())

# Read sample coordinate file and define coordinate system
Geo_raw = "./coords_sorted.txt"
sample.points <- read.table(file = Geo_raw, sep = '\t', header = T, row.names = 1, fill = T)
str(sample.points)
d <- data.frame(lon=sample.points[,2], lat=sample.points[,1]) #lon= column Y, #lat= column X
library(sp)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326")
CRS.new <- CRS("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )

library(rgdal)
d.new <- spTransform(d, CRS.new)
sample.points.new <- sample.points
sample.points.new[,1] = d.new$lon
sample.points.new[,2] = d.new$lat

par(mfrow=c(1,3))
plot.default(d$lon,d$lat, main="Raw data", cex.axis=.95)
plot(d, axes=TRUE, main="Original lat-lon", cex.axis=.95)
plot(d.new, axes=TRUE, main="Projected", cex.axis=.95)
dev.off()

unclass(d.new)
write.table(sample.points.new, file = './coords_converted.txt', sep = '\t')

##############################
# ResDisMapper
##############################

library(ResDisMapper)

# Ensure that rows of both files match, if not refer to sorting section in https://github.com/takfung/ResDisMapper/blob/master/Documentation/ResDisMapper_manual_1.1.pdf

# Read genepop and coordinate files
# Must add the word "Genetic" to the first cell (first row + first column) of genepop file, if not it would not work
Gen_raw = "./albo_genepop.gen" 
Geo_raw = "./coords_converted.txt"

# Calculate residuals
IBD.res <- rdm_IBD(Gen_raw, Geo_raw, Dist_method = 1, IBD_method = 1) 

# Testing rdm_residual
# carry out spatial autocorrelation 
res_SLDF <- rdm_residual(IBD.res, Geo_raw, min.dist = 1, max.dist = 20000, n_resolution = 50, proj=sp::CRS("+init=epsg:4326"))

# Testing rdm_resistance
F.df<-rdm_resistance(IBD.res, res_SLDF, nrows = 50, ncols = 100, 
                     conf_intervals = 0.95, random_rep = 1000, 
                     outputfile = "./resistance_map_20km.csv")

tiff("./resistance_plot_20km.tiff", width = 15, height = 6, units = 'in', res = 300)
rdm_mapper(F.df, Geo_raw, r_size = 5, p_signf = 0.05, p_size = 2, p_col = "yellow", disp_all_cells = 0, disp_contours = 1)
dev.off()

# To plot layer with different colours, specify the colours in the function below and run the function again
rdm_mapper2 <- function(F.df, Geo_raw, p_signf = 0.05, p_size = 2, p_col = "yellow", disp_all_cells = 0, disp_contours = 1){
  
  sample.points <- read.table(file = Geo_raw, sep = '\t', header = T, row.names = 1, fill = T)
  
  names(sample.points) <- c('x', 'y')
  
  p_signf_u <- 1 - p_signf
  
  resistance2 = F.df$resistance
  if(disp_all_cells==0){
    sign.check = F.df$sign.check
    sign.check_eqminus1_ind = which(sign.check==-1)
    resistance2[sign.check_eqminus1_ind] = 0
    F.df$resistance = resistance2
  }else{
    F.df$resistance = resistance2
  }
  
  if(disp_contours==0){
    ggplot2::ggplot(F.df, ggplot2::aes(x = x, y = y, colour=resistance)) + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))+
      ggplot2::geom_raster(ggplot2::aes(fill=resistance)) +
      ggplot2::scale_fill_gradient2(low="#9e9e9e", high="#FF8A65", mid="white", midpoint=0)+
      ggplot2::geom_point(data = sample.points, cex = p_size, shape = 21, color = "black", fill = p_col,stroke = 1)
  }else{
    ggplot2::ggplot(F.df, ggplot2::aes(x = x, y = y, colour=resistance)) + 
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))+
      ggplot2::geom_raster(ggplot2::aes(fill=resistance)) +
      ggplot2::scale_fill_gradient2(low="#9e9e9e", high="#FF8A65", mid="white", midpoint=0)+
      ggplot2::stat_contour(ggplot2::aes(z = Prob), bins=4, size=1, col="#FF8A65", breaks = c(p_signf_u))+
      ggplot2::stat_contour(ggplot2::aes(z = Prob), bins=4, size=1, col="#9e9e9e", breaks = c(p_signf))+
      ggplot2::geom_point(data = sample.points, cex = p_size, shape = 21, color = "black", fill = p_col,stroke = 1)
  }
  
}