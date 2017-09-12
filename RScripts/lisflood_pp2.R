# LISFLOOD pre-processing 2
#
# NOT WORKING. USE PYTHON SCRIPT, INSTEAD.
#
# INPUTS
# Rasters: flow direction, flow accumulation, widths
# VIC routing model output - discharge time series at upstream boundaries
#
# OUTPUTS
# Rasters: river, boundary_points, labels, chainage
# smoothed elevation raster (riverband)
# smoothed elevation raster with burned rivers (riverbed)
# initial depth (80% bank full, initial wd) raster

library(raster)
library(rgdal)
library(rgeos)

setwd("/Users/jschapMac/Desktop/TuoSub/GIS")
# dirpath <- "/Users/jschapMac/Documents/HydrologyData/HydroSHEDs/na_dir_15s_grid/na_dir_15s/na_dir_15s"
# accpath <- "/Users/jschapMac/Documents/HydrologyData/HydroSHEDs/na_acc_15s_grid/na_acc_15s/na_acc_15s"
dirpath <- getwd()
accpath <- getwd()
bbpath <- getwd()

# Calculate flowdir, correct flowdir, calculate flowacc. Do this using the
# upscaled, filled DEM produced in lisflood_pp1.R.

# Import flowdir, flowacc
# Use VIC convention for flowdir (1-8)
fdir <- raster(file.path(dirpath, "flowdir_topotoolbox_1-8.asc"))
facc <- raster(file.path(accpath, "flowacc_topotoolbox.asc"))

# Import widths raster
rwidth <- raster("rwidth.tif")

# Create a (0,1) raster called river consisting of all nonzero cells in widths
river <- rwidth > 0

# Using the flow direction and river rasters, create a raster called 
# boundary_points.
boundary_pts <- find_boundary_points(river, flowdir)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Some functions
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
find_boundary_points <- function(river, flowdir) {
  
  bnd1 <- vector()
  bnd2 <- vector()
  
  nr <- nrow(river); nc <- ncol(river)
  priver <- array(data = 0, dim = c(nr+2, nc+2))
  priver[2:dim(priver)[1]-2, 2:dim(priver)[2]-2] <- river[,,1]
  pflowdir <- array(data = 0, dim = c(nr+2, nc+2))
  pflowdir[2:dim(pflowdir)[1]-2, 2:dim(pflowdir)[2]-2] <- flowdir[,,1]
  
  for (i in 1:nc) {
    for (j in 1:nr) {
      if (river[i,j]) {
        if (!has_inflow(priver, pflowdir, i+1, j+1)) {
          bnd1 < append(bnd1, i)
          bnd2 < append(bnd1, j)
        }
      }
    }
  }
  
  bnd <- cbind(bnd1, bnd2)
  return(bnd)
}

# ---------------------------------------------------------------------------

has_inflow <- function(priver, pflowdir, i, j) {
  


  
  
  return(flows)
}

# ---------------------------------------------------------------------------

Create labels raster
Create chainage raster

Import elevation raster (at model resolution, cropped to basin extent)
Smooth bank heights using elevation raster, labels, and chainage as inputs to a Lowess location regression algorithm.
Import depths raster
Set negative depths to zero
Subtract depths from the smoothed elevation raster
Save the burned elevation raster 

Set a new raster equal to depths times 0.8
Set pixels with negative widths to zero
Save the resulting raster as initial water depth

Find upstream and lateral inflow points using river, chainage, labels, threshold, and flowacc as inputs
Use upstream and lateral inflow points to create bci file [inflows, upper left x, upper left y, xres, yres, nrows, ncols]. bci file may only include boundary points, check.

Save smoothed DEM as ASCII ArcInfo Grid at model resolution
Save widths raster as "
Save riverbed raster as "
Save initwd raster as "

Enter x y coordinates of stream gauge locations/location where you would like to produce river flow time series
Save as text file

Make bdy file using outputs from VIC routing model and locations in bci file


# ---------------------------------------------------------------------------

# Just a test to see how problematic it really is to upscale flowacc
mres <- 1000
merc <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
facc.merc <- projectRaster(facc, res = mres, crs = merc, method = "ngb",
progress = "text")
facc.cropped <- crop(facc, extent(bb.merc))