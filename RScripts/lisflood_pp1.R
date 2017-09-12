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

setwd("/Users/jschapMac/Desktop/TuoSub/GIS")
fillpath <- "/Users/jschapMac/Desktop/Tuolumne5/LF2/HydroSHEDs/na_dem_15s_grid/na_dem_15s/na_dem_15s"
bbpath <- "/Users/jschapMac/Desktop/TuoSub/GIS"
rivpath <- "/Users/jschapMac/Documents/HydrologyData/HydroSHEDs/na_riv_15s"
rivdatpath <- "/Users/jschapMac/Documents/HydrologyData/HydroSHEDs/width_depth_data/namerica"

mres <- 1000 # modeling resolution (m)

# Import HydroSHEDS 15 arc-second filled DEM for North America
fill <- raster(file.path(fillpath, "w001001.adf"))

# Project to Mercator coordinate system (epsg 3395)
merc <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
fill.merc <- projectRaster(fill, res = mres, crs = merc, method = "bilinear",
                           progress = "text", filename = "fill_merc.asc")

# Or just load the projected, filled DEM:
# fill.merc <- raster("/Users/jschapMac/Desktop/Tuolumne5/LF2/fill_merc.asc)

# Import basin shapefile
ogrListLayers(file.path(bbpath, "tuolumne_headwaters.shp"))
bb <- readOGR(dsn = bbpath, layer = "tuolumne_headwaters", verbose = F)

# Project to Mercator coordinate system (epsg 3395)
bb.merc <- spTransform(bb, CRS(merc))

# Crop DEM to extent of basin shapefile
fill.merc.cropped <- crop(fill.merc, extent(bb.merc))

# Clip DEM to basin boundary
fill.merc.clipped <- mask(fill.merc.cropped, bb.merc)

# Import river width and depth data
ogrListLayers(file.path(rivdatpath, "narivs.shp"))
narivs <- readOGR(dsn = rivdatpath, "narivs", verbose = F)
narivs.merc <- spTransform(narivs, CRS(merc))
narivs.merc.clipped <- intersect(x = narivs.merc, y = bb.merc)
# The above steps can take a long time

# Check the data contained in the river dataset
plot(narivs.merc.clipped)
names(narivs.merc.clipped@data)

# ------------------------------------------------------------
# Select subset of river network with more than a threshold 
# number of contributing pixels. Note: this is only OK if 
# there are no gaps in the resulting river network. 
thres <- 100
ind <- as.numeric(narivs.merc.clipped@data$AREA) > thres
plot(narivs.merc.clipped[ind,])
plot(narivs.merc.clipped)
narivs.merc.thres <- narivs.merc.clipped[ind,]

# ------------------------------------------------------------

# Make width and depth rasters
rwidth <- rasterize(narivs.merc.thres, y = fill.merc.clipped, field = 
                      "WIDTH", progress = "text")

rdepth <- rasterize(narivs.merc.thres, y = fill.merc.clipped, field = 
                      "DEPTH", progress = "text")

# ------------------------------------------------------------
# Save outputs

# cropped DEM (filled)
writeRaster(fill.merc.cropped, filename = "fill_merc_cropped.asc", 
            format = "ascii", overwrite = T)

# clipped DEM (filled)
writeRaster(fill.merc.clipped, filename = "fill_merc_clipped.asc", 
            format = "ascii", overwrite = T)

# full projected river network
writeOGR(narivs.merc, dsn = getwd(), layer = "narivs.merc", 
         driver = "ESRI Shapefile")

# clipped river network
writeOGR(narivs.merc.clipped, dsn = getwd(), layer = "narivs.merc.clipped", 
         driver = "ESRI Shapefile")

# clipped river network (using threshold)
writeOGR(narivs.merc.thres, dsn = getwd(), layer = "narivs.merc.thres", 
         driver = "ESRI Shapefile")

# basin boundary (Mercator)
writeOGR(bb.merc, dsn = getwd(), layer = "bb.merc", 
         driver = "ESRI Shapefile")

# width and depth rasters
writeRaster(rdepth, filename = "rdepth.tif", format = "GTiff", overwrite = T)
writeRaster(rwidth, filename = "rwidth.tif", format = "GTiff", overwrite = T)

# ---------------------------------------------------------------------------
# Save flowdir and flowacc as GeoTIFFs
# (import them first)
writeRaster(flowdir, "flowdir_topotoolbox_1-8.tif", format = "GTiff")
writeRaster(flowacc, "flowacc_topotoolbox.tif", format = "GTiff")
# There may be an issue if flowacc is derived with TopoToolbox because the
# format is different from the format of an ArcMap-derived flowacc map.

# ---------------------------------------------------------------------------
# Attach width and depth attributes from narivs to the stream network 
# extracted from the flowacc file

# Do this in ArcMap

# Load the streamlines_w_width vector and rasterize it. Same for depth.
streamspath <- "/Users/jschapMac/Desktop/TuoSub/GIS/latest_set_of_hydrography_data/streamlines_w_width"
ogrListLayers(dsn = streamspath)
streamsww <- readOGR(dsn = streamspath, layer = "streamlines_w_width")

# Make width and depth rasters
rwidth <- rasterize(streamsww, y = fill.merc.clipped, field = 
                      "WIDTH", progress = "text")

rdepth <- rasterize(streamsww, y = fill.merc.clipped, field = 
                      "DEPTH", progress = "text")

# width and depth rasters
writeRaster(rdepth, filename = "rdepth.tif", format = "GTiff", overwrite = T)
writeRaster(rwidth, filename = "rwidth.tif", format = "GTiff", overwrite = T)
