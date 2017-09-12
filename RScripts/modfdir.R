# Functions for converting flow direction coordinate system among different
# softwares' conventions

library(raster)
library(rgdal)

setwd("/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/GIS")
fdir <- raster("tuolumne_flowdir.tif")

arc2vic <- function(fdir) {
  fdir.copy <- fdir
  
  fdir.copy[fdir == 1] <- 3
  fdir.copy[fdir == 2] <- 4
  fdir.copy[fdir == 4] <- 5
  fdir.copy[fdir == 8] <- 6
  fdir.copy[fdir == 16] <- 7
  fdir.copy[fdir == 32] <- 8
  fdir.copy[fdir == 64] <- 1
  fdir.copy[fdir == 128] <- 2
  
  fdir <- fdir.copy
}

vic2grass <- function(fdir) {
  fdir.copy <- fdir
  
  fdir.copy[fdir == 1] <- 2
  fdir.copy[fdir == 2] <- 1
  fdir.copy[fdir == 3] <- 8
  fdir.copy[fdir == 4] <- 7
  fdir.copy[fdir == 5] <- 6
  fdir.copy[fdir == 6] <- 5
  fdir.copy[fdir == 7] <- 4
  fdir.copy[fdir == 8] <- 3
  
  fdir <- fdir.copy
}

grass2vic <- function(fdir) {
  # My version
  fdir.copy <- abs(fdir)
  
  fdir.copy[abs(fdir) == 2] <- 1
  fdir.copy[abs(fdir) == 1] <- 2
  fdir.copy[abs(fdir) == 8] <- 3
  fdir.copy[abs(fdir) == 7] <- 4
  fdir.copy[abs(fdir) == 6] <- 5
  fdir.copy[abs(fdir) == 5] <- 6
  fdir.copy[abs(fdir) == 4] <- 7
  fdir.copy[abs(fdir) == 3] <- 8
  
  fdir <- fdir.copy
}

grass2vic_ka <- function(fdir) {
  # Kostas' version
  fdir.copy <- abs(fdir)
  
  fdir.copy[abs(fdir) == 2] <- 1
  fdir.copy[fdir == 1] <- 2
  fdir.copy[fdir == -1] <- 1
  fdir.copy[abs(fdir) == 8] <- 3
  fdir.copy[abs(fdir) == 7] <- 4
  fdir.copy[abs(fdir) == 6] <- 5
  fdir.copy[abs(fdir) == 5] <- 6
  fdir.copy[abs(fdir) == 4] <- 7
  fdir.copy[abs(fdir) == 3] <- 8
  
  fdir <- fdir.copy
}

fdir <- grass2vic_ka(fdir)
writeRaster(fdir, "tuolumne_flowdir_ka.asc", format = "ascii", datatype = "INT2S")
writeRaster(fdir, "flowdir_vic_input_vic_convention.tif", format = "GTiff", datatype = "INT2S")
