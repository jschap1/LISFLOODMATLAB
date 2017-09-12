# Calculate the number of 1000 m pixels with one 1/8 degree pixel in my study region
# For a large region, this number will not be constant
#
# INCOMPLETE

library(raster)

coarse <- raster("tuolumne.fract")
fine <- raster("basin.tif")
values(coarse) <- 1
values(fine) <- 1

intersect(fine, coarse[1,1])
