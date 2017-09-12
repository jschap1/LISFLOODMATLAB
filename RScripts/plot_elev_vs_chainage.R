# Plot elevation vs. chainage

library(raster)

setwd("/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/GIS")

felev <- raster("tuolumne_elev.tif")
labels <- raster('labels.tif')
chainage <- raster("chainage.tif")
facc <- raster("tuolumne_acc.tif")

river <- labels>0 
facc <- abs(facc*river)

nbnd <- 6 # number of boundary points
par(mfrow = c(3,2))

for (i in 1:nbnd) {
  
  riv <- facc*(labels == i)
  riv[riv==0] <- NA
  
  facc.vals <- riv[!is.na(riv)]
  felev.vals <- felev[!is.na(riv)]
  chainage.vals <- chainage[!is.na(riv)]
  df <- as.data.frame(cbind(facc.vals, felev.vals, chainage.vals))
  ordered.df <- df[order(df$facc.vals),]
  plot(ordered.df$chainage.vals,ordered.df$felev.vals, type="l", main = paste("Segment", i), 
       xlab = "Chainage (m)", ylab = "Elevation (m)", col = "blue", ylim=c(0,60))
  
}

# ----------------------------------------------------------------------------------------
# Repeat, but show bed elevation and smoothed elevation

bed <- raster("tuolumne_bed.tif")
selev <- raster("tuolumne_selev.tif")

i = 1

riv <- facc*(labels == i)
riv[riv==0] <- NA

facc.vals <- riv[!is.na(riv)]
bed.vals <- bed[!is.na(riv)]
selev.vals <- selev[!is.na(riv)]
chainage.vals <- chainage[!is.na(riv)]
df <- as.data.frame(cbind(facc.vals, selev.vals, bed.vals, chainage.vals))
ordered.df <- df[order(df$facc.vals),]
plot(ordered.df$chainage.vals,ordered.df$selev.vals, type="l", main = paste("Segment", i), 
     xlab = "Chainage (m)", ylab = "Elevation (m)", col = "blue", ylim=c(0,60))
lines(ordered.df$chainage.vals,ordered.df$bed.vals, col = "red")
legend("topright", c("Bank elev.","Bed elev."), fill = c("blue","red"))





