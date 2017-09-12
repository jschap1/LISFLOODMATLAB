# LISFLOOD pre-processing 1
#
# INPUTS
# HydroSHEDS filled DEM raster
# HydroSHEDS-based river width, depth, flow network shapefile
# HUC basin boundary shapefile
# 
# OUTPUTS
# DEM at modeling resolution, cropped to basin extent
# DEM at modeling resolution, clipped to basin boundary
# Flow network shapefile, clipped to basin boundary, using threshold upstream area
# Width and depth rasters
# Basin boundary (Mercator)

library(raster)
library(rgdal)
library(rgeos)

setwd("/Users/jschapMac/Desktop/Tuolumne/TuoSub/GIS")
dempath <- "/Users/jschapMac/Documents/Data/HydroSHEDs/na_dem_15s_grid/na_dem_15s/na_dem_15s"
bbpath <- "/Users/jschapMac/Desktop/Tuolumne/TuoSub/GIS"
rivdatpath <- "/Users/jschapMac/Documents/Data/HydroSHEDs/width_depth_data/namerica"

mres <- 1000 # modeling resolution (m)

# Import HydroSHEDS 15 arc-second filled DEM for North America
dem <- raster(file.path(dempath, "w001001.adf"))

# Project to Mercator coordinate system (epsg 3395)
merc <- CRS("+init=epsg:3395")
dem.merc <- projectRaster(dem, res = mres, crs = merc, method = "bilinear",
                           progress = "text", filename = "dem_merc.asc")

# Or just load the projected, filled DEM:
# dem.merc <- raster("/Users/jschapMac/Desktop/Tuolumne5/LF2/dem_merc.asc)

# Import basin shapefile
ogrListLayers(file.path(bbpath, "tuolumne_headwaters.shp"))
bb <- readOGR(dsn = bbpath, layer = "tuolumne_headwaters", verbose = F)

# Project to Mercator coordinate system (epsg 3395)
bb.merc <- spTransform(bb, merc)

# Crop DEM to extent of basin shapefile
dem.merc.cropped <- crop(dem.merc, extent(bb.merc))

# Clip DEM to basin boundary
dem.merc.clipped <- mask(dem.merc.cropped, bb.merc)

# Import river width and depth data
ogrListLayers(file.path(rivdatpath, "narivs.shp"))
narivs <- readOGR(dsn = rivdatpath, "narivs", verbose = F)
narivs.merc <- spTransform(narivs, merc)
narivs.merc.clipped <- intersect(x = narivs.merc, y = bb.merc)
# The above steps can take a few minutes

# Check the data contained in the river dataset
plot(bb.merc)
plot(narivs.merc.clipped, add = T)
names(narivs.merc.clipped@data)

# Change nodata value in the DEM
dem.merc.clipped.na <- dem.merc.clipped
dem.merc.clipped[is.na(dem.merc.clipped)] <- -9999

# ------------------------------------------------------------
# Save outputs

# cropped DEM (filled)
writeRaster(dem.merc.cropped, filename = "dem_merc_cropped.asc", 
            format = "ascii", overwrite = T)

# clipped DEM (filled)
writeRaster(dem.merc.clipped, filename = "dem_merc_clipped.asc", 
            format = "ascii", overwrite = T)

# clipped DEM (filled)
writeRaster(dem.merc.clipped, filename = "dem_merc_clipped.tif", 
            format = "GTiff", overwrite = T)

# full projected river network
writeOGR(narivs.merc, dsn = getwd(), layer = "narivs.merc", 
         driver = "ESRI Shapefile")

# clipped river network
writeOGR(narivs.merc.clipped, dsn = getwd(), layer = "narivs.merc.clipped", 
         driver = "ESRI Shapefile")

# basin boundary (Mercator)
writeOGR(bb.merc, dsn = getwd(), layer = "bb.merc", 
         driver = "ESRI Shapefile")

# ------------------------------------------------------------

# Do stuff in ArcMap, then finish running this script

ogrListLayers("TempGIS/streamlines_w_chars.shp")
streamsww <- readOGR(dsn = "TempGIS/streamlines_w_chars.shp", 
                     layer = "streamlines_w_chars")

# Rasterize the streamlines vector

# Make width and depth rasters
rwidth <- rasterize(streamsww, y = dem.merc.clipped, field = 
                      "Avg_WIDTH", progress = "text")

rdepth <- rasterize(streamsww, y = dem.merc.clipped, field = 
                      "Avg_DEPTH", progress = "text")

# Save width and depth rasters
writeRaster(rdepth, filename = "rdepth.tif", format = "GTiff", overwrite = T)
writeRaster(rwidth, filename = "rwidth.tif", format = "GTiff", overwrite = T)

# Convert flow direction, flow accumulation to GeoTiffs and save

fdir <- raster("TempGIS/mod_fdir/w001001.adf")
facc <- raster("TempGIS/mod_facc/w001001.adf")

# Change flow direction numbering convention in modfdir.R
fdir <- raster("modfdir")

writeRaster(fdir, filename = "fdir.tif", datatype = "INT2S", 
            format = "GTiff", overwrite = T)
writeRaster(facc, filename = "facc.tif", datatype = "INT2S",
            format = "GTiff", overwrite = T)

# ------------------------------------------------------------
# Load bci points and plot them to check if the boundary points are OK

bci <- read.table("TuoSub.bci")
names(bci) <- c("V1", "Easting", "Northing", "V4", 'V5')
bndpts <- cbind(bci$Easting[c(-1,-2,-3,-4)], bci$Northing[c(-1,-2,-3,-4)])
bnd.sp <- SpatialPoints(coords = bndpts, proj4string = merc)
FID <- 1:dim(bnd.sp@coords)[1]
bnd.sp <- SpatialPointsDataFrame(coords = bnd.sp, data = as.data.frame(FID))
writeOGR(bnd.sp, dsn = getwd(), 
         layer = "bnd.sp", driver = "ESRI Shapefile", overwrite_layer = TRUE)


plot(dem.merc.clipped.na)
plot(streamsww, add=T)
plot(bnd.sp, add=T)

######
#####
###
##
# Figure out what is going on with the bndpts:
fdir <- raster("fdir.tif")
widths <- raster("rwidth.tif")
river <- widths > 0

writeRaster(river,"river.asc", format = "ascii",datatype = "INT2S")

