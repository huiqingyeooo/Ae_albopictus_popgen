# Landscape genetic optimization with ResistanceGA

rm(list = ls())
library(raster)
library(rgdal)

###################################
# Converting shape file to raster
###################################
# Make sure that all raster layers intended for ResistanceGA have the same resolution
# Start with generating a base raster in the resolution that you need. Here I import an existing raster layer for consistency.

base <-raster("surface.tif")
plot(base)
base[!is.na(base)] <- 0
plot(base)

# Read in shape file
aegPop<-readOGR(dsn="./shape_files",layer="areas_with_high_aedes_population")
CRS.new <- CRS("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" )
aegPop.new= spTransform(aegPop, CRS.new)

# Select features
aegPop.new@data$tessellate
aegPop.sel = aegPop.new[grepl(c("-1"), aegPop.new@data$tessellate), ]

# Rasterize with base raster
aegPop.raster <- rasterize(aegPop.sel, base)
aegPop.raster[!is.na(aegPop.raster)] <- 1
aegPop.final<-merge(aegPop.raster,base)

# Check plot and resolution
plot(aegPop.final)
aegPop.final

# Write to .asc format
writeRaster(aegPop.final, "aegPop.asc", format="ascii",overwrite=TRUE)

###################################
# ResistanceGA Optimization
###################################
rm(list = ls())

# Write all the GIS layers (rasters) in the same folder as .asc files
library(raster)
library(rgdal)

# Calculate genetic distance with a genepop file using poppr
library(poppr)
library(ResistanceGA)

gen<-read.genepop("albo_genepop.gen", ncode=3) #check first line, dataset title should be present
gen_dist <- diss.dist(gen)
gen_dist<-as.vector(gen_dist)
length(gen_dist)

# Single layer optimization using Circuitscape (download circuitscape beforehand)
GA.inputs <- GA.prep(ASCII.dir = "./asc")
CS.program <- paste("./Circuitscape/cs_run.exe")
CS.inputs <- CS.prep(n.Pops = 105, response = gen_dist, 
                     CS_Point.File = "./samples.txt", 
                     CS.program = CS.program)
SS_results_vegRC <- SS_optim(CS.inputs = CS.inputs, GA.inputs = GA.inputs)

# Codes for quick plotting to check asc layers and resistance maps
result.map <- raster("./aegPop.asc")
plot(result.map, legend=TRUE)
summary(result.map)
unique(result.map)