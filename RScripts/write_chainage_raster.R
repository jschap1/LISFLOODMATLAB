setwd("/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/GIS")

# Save chainage as text file using numpy.savetxt() in Python

chainage <- read.table("chainage.txt")
basin <- raster("basin.tif")
chainage.rast <- raster(basin)

values(chainage.rast) <- data.matrix(chainage)

chainage.rast[chainage.rast==0] <- NA

writeRaster(chainage.rast, "chainage.tif", format="GTiff", overwrite = T)
