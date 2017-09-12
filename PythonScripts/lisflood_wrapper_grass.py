# LISFLOOD input file preparation.
# Created by Kostas Andreadis.
# Modified 8/31/2017 JRS
#
# Code is divided into sections run from the command prompt, GRASS GIS, and Python. In each case, keep the session open until the whole code has been run. Some GRASS commands seem to need to be run in the GRASS console and/or loaded into the map display in order to work. I am not sure why, but it seems kind of fussy. It is recommended to plot the map layers in GRASS as you go. The baseline code assumes you are using a 15 arc-sec HydroSHEDS DEM.
#
# Note: You cannot just run this script mindlessly. Requires careful reading through and babysitting/hand-holding as it is run :(
#
# INPUTS
# DEM
# Basin boundary
# River network dataset
#
# OUTPUTS
# VIC routing model input files
# LISFLOOD input files
#
# Useful EPSG codes:
# 4269 NAD83
# 4326 WGS84
# 3395 World Mercator

# --------------------------------------------------------------
# SYSTEM

cd /Users/jschapMac/Desktop/Tuolumne/Tuolumne7/GIS

# Transform DEM to GeoTIFF
gdal_translate -of GTiff "/Users/jschapMac/Documents/Data/HydroSHEDs/na_dem_15s_grid/na_dem_15s/na_dem_15s/w001001.adf" "na_dem_15s.tif"

# Reproject DEM to World Mercator
gdalwarp -dstnodata -9999 -t_srs epsg:3395 na_dem_15s.tif na_dem_15s_merc.tif

# Get extent of basin boundary in World Mercator
ogrinfo -al -so /Users/jschapMac/Desktop/Tuolumne/Tuolumne2/Shapefiles/upper_tuolumne.shp | grep Extent | awk -F'[(,]' '{printf("%f %f\n%f %f\n",$2,$5,$4,$3)}' | gdaltransform -s_srs epsg:4269 -t_srs epsg:3395 -output_xy
# Copy the output into the next command:

# Crop the DEM to the basin extent
gdal_translate -a_nodata -9999 -projwin -13489113.2481141 4585456.99715751 -13269202.5959274 4493210.47482815 na_dem_15s_merc.tif tuolumne_dem.tif

# Reproject the basin boundary to World Mercator
ogr2ogr -t_srs EPSG:3395 bb_merc.shp /Users/jschapMac/Desktop/Tuolumne/Tuolumne2/Shapefiles/upper_tuolumne.shp

# Reproject the river dataset to World Mercator and crop it
ogr2ogr -t_srs EPSG:3395 -clipdst -13489113.2481141 4585456.99715751 -13269202.5959274 4493210.47482815 narivs_merc.shp /Users/jschapMac/Documents/Data/HydroSHEDs/width_depth_data/namerica/narivs.shp

# --------------------------------------------------------------
# GRASS

# Create a LOCATION called Tuolumne.merc with the World Mercator projection

# Enter the following commands to do hydrologic analysis:
cd /Users/jschapMac/Desktop/Tuolumne/Tuolumne7/GIS

v.in.ogr in=bb_merc.shp out=basin

# This step seems to only work if the basin is displayed in GRASS
v.to.rast --overwrite in=basin out=basin use=val value=1

r.in.gdal --overwrite in=tuolumne_dem.tif out=hydrosheds.dem
g.region rast=hydrosheds.dem
g.region res=1000
r.resamp.stats --overwrite in=hydrosheds.dem out=hydrosheds.elev method=minimum
r.hydrodem -a --overwrite in=hydrosheds.elev out=hydrosheds.felev

# Choose an appropriate flow accumulation threshold to generate the stream network.

# --------------------------------------------------------------
# EXPERIMENTAL SECTION -> Skip if threshold is known
# Following Tarboton et al. (1991), we use the stream drop method to determine the upstream area where the constant drop law fails:

# Try a threshold, say 100, for the number of upstream pixels

# Generate flow direction, flow accumulation, and stream rasters using this threshold.
r.watershed -s elev=hydrosheds.felev accum=facc400 stream=river400 drain=fdir400 threshold=400

# Calculate stream order
r.stream.order stream_rast=river400 direction=fdir400 strahler=so400

# Export, and move to R
r.out.gdal in=hydrosheds.felev out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/fill.tif nodata=-9999
r.out.gdal in=facc400 out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/facc400.tif nodata=-9999
r.out.gdal in=river400 out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/river400.tif
r.out.gdal in=so400 out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/so400.tif

# (Do the following in R)
# Use a t-test of means to determine whether the stream drop stays constant across stream order at this area threshold.
# Choose the smallest area threshold where the constant stream drop property holds true.

# --------------------------------------------------------------

r.watershed --overwrite -s elev=hydrosheds.felev accum=hydrosheds.flowacc stream=hydrosheds.river drain=hydrosheds.flowdir threshold=400

r.out.gdal --overwrite in=hydrosheds.flowacc out=tuolumne_acc.tif nodata=-9999
r.out.gdal --overwrite in=hydrosheds.flowdir out=tuolumne_flowdir.tif nodata=-9999
r.out.gdal --overwrite type=Float64 in=hydrosheds.felev out=tuolumne_elev.tif nodata=-9999
r.out.gdal --overwrite type=Float64 in=hydrosheds.river out=tuolumne_river.tif nodata=-9999

# If there is an error about SetColorTable(), it can be ignored. It has do with GeoTIFF format not being able to reproduce colortables for floating-point data. See https://gis.stackexchange.com/questions/101926/grass-raster-map-types-uint16-and-float64.

# If you are going to do flow direction corrections, now is a good time -> but how to reincorporate the corrected flow direction into the workflow?
# Burning the river network into the DEM in flat areas might be useful, as an alternative.
# For now, just use the fdir results from GRASS directly.

# Import the width and depth database and attach the attributes to the nearest river segment in the derived network.

v.in.ogr --overwrite -r in=narivs_merc.shp out=nariv where='AREA>400'
# Use the same threshold for AREA as the threshold in r.watershed above

r.mask basin
r.to.vect --overwrite in=hydrosheds.river out=river type=line -s
r.mask -r
v.db.addcolumn map=river col='width real'
v.db.addcolumn map=river col='depth real'
v.distance from=river to=nariv upload=to_attr to_column=DEPTH column=depth
v.distance from=river to=nariv upload=to_attr to_column=WIDTH column=width
####### OUT: vector river width and depth

# Rasterize the river vector to get the width and depth rasters.
v.to.rast --overwrite in=river out=hydrosheds.width use=attr attribute_column=width
v.to.rast --overwrite in=river out=hydrosheds.depth use=attr attribute_column=depth
r.out.gdal --overwrite in=hydrosheds.width out=tuolumne_widths.tif nodata=-9999
r.out.gdal --overwrite in=hydrosheds.depth out=tuolumne_depths.tif nodata=-9999
r.out.gdal --overwrite in=basin out=basin.tif

# --------------------------------------------------------------
# PYTHON

####### The next steps involve preparing the actual LISFLOOD-FP input files. Given the uncertainty in the DEM, we will smooth the river bank heights to avoid any numerical instabilities during simulation. This is accomplished by first identifying the upstream and downstream boundary points in the domain, labeling the river channels and calculating each one's chainage (i.e. distance downstream).
import sys
sys.path.append("/Users/jschapMac/Documents/Codes/VICMATLAB/PythonScripts")
import os

import lisflood
flowdir = lisflood.read_raster("tuolumne_flowdir.tif")
flowacc = lisflood.read_raster("tuolumne_acc.tif")
widths = lisflood.read_raster("tuolumne_widths.tif")
river = widths > 0
bndpts = lisflood.find_boundary_points(river, flowdir)
labels, chainage = lisflood.calc_chainage(river, flowdir, abs(flowacc), bndpts, 1000)
####### river, bndpts, labels, chanage



####### The smoothing of the bank heights is done by using a LOWESS local regression, and the channel is burned in to the DEM by subtracting the depth raster.
elev = lisflood.read_raster("tuolumne_elev.tif")
selev = lisflood.smooth_bank_heights(elev, labels, chainage)
lisflood.write_raster(selev, "tuolumne_selev.tif", "tuolumne_elev.tif")
depths = lisflood.read_raster("tuolumne_depths.tif")
depths[depths < 0] = 0.0
belev = selev - depths
lisflood.write_raster(belev, "tuolumne_bed.tif", "tuolumne_elev.tif")
y0 = 0.8 * depths
y0[widths <= 0] = 0.0
lisflood.write_raster(y0, "tuolumne_initwd.tif", "tuolumne_elev.tif")
####### OUT: smoothed riverbank(arkansas_selev.tif), riverbed(arkansas_bad.tif), initial depth(80% bank full, arkansas_initwd.tif)



####### Next we identify the upstream and lateral inflow points, and generate the BCI file.
from osgeo import gdal
inflows = lisflood.identify_inflows(river, chainage, labels, abs(flowacc), 400)
f = gdal.Open("tuolumne_selev.tif")
xul, xres, _, yul, _, yres = f.GetGeoTransform()
f = None
nrows, ncols = elev.shape

# In lisflood.write_bci, there is a place where the prefix is specified for naming the boundary inflow points. This needs to be consistent with the prefixes used elsewhere so there is correspondence between the bci and bdy files. The default prefix is INF, as in INF001, INF002, etc.
lisflood.write_bci("/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Inputs/tuolumne.bci", inflows, xul, yul, xres, yres, nrows, ncols)
####### OUT: BCI file that include boundary points (arkansas.bci)

# ----------------------------------------------------------------------
# SYSTEM

####### Then we generate the DEM and sub-grid channel width files.
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 tuolumne_selev.tif /Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Inputs/tuolumne.dem
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 tuolumne_widths.tif /Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Inputs/tuolumne.width
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 tuolumne_bed.tif /Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Inputs/tuolumne.bed
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 -co DECIMAL_PRECISION=2 tuolumne_initwd.tif /Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Inputs/tuolumne.initwd

# ----------------------------------------------------------------------
# PYTHON

####### If we need LISFLOOD-FP to produce a time series of river flow at specific locations, we need to generate a virtual gauge file.

# Because of the way lisflood.write_gauges is set up, the gauges must be at points flowing in a non-diagonal direction (NEWS).

x = [-13479618, -13397374, -13400626]
y = [4494493, 4529860, 4529745]
lisflood.write_gauges(x, y, flowdir, "tuolumne_widths.tif", "/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Inputs/tuolumne.gauge")
####### OUT: gage file (arkansas.gauge)


####### The final step is generating the BDY file for LISFLOOD-FP, which contains the streamflow at the inflow points (in m2/s). For this simulation we will use the output of the VIC routing model forced by the NLDAS-2 VIC model output. We need to prepare the input files for the routing model, describing the flow direction, flow fraction and station locations. The station file is created from the inflows identified
from pyproj import Proj
wmerc = Proj("+init=EPSG:3395")
fout = open("tuolumne.stations", 'w')
for sta, xy in enumerate(inflows):
    # Enter xll, yll for the region, in lat/lon coordinates, as well as the VIC modeling resolution (degrees)
    x, y = wmerc(xul+xres*xy[1], yul+yres*xy[0], inverse=True)
    xi = int((x-(-121.1875))/0.0625) + 1
    yi = int((y-37.5625)/0.0625) + 1
    # The prefix AA{0:03d} is used for the bci/bdy files to talk to each other. They must correspond exactly.
    fout.write("1\tTUO{0:03d}\t{1}\t{2}\t-9999\nNONE\n".format(sta, xi, yi))
fout.close()

# --------------------------------------------------------------
# EXPERIMENTAL SECTION -> Use to determine threshold for VIC routing model. Skip if a good threshold is known

# Generate flow direction, flow accumulation, and stream rasters using an arbitrary threshold.
r.watershed -s elev=hydrosheds.felev accum=facc2 drain=fdir2 stream=river2 threshold=2

# Calculate stream order
r.stream.order stream_rast=river2 direction=fdir2 strahler=so2

# Export, and move to R to determine threshold
r.out.gdal in=hydrosheds.felev out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/vic/fill.tif nodata=-9999
r.out.gdal in=facc2 out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/vic/facc2.tif nodata=-9999
r.out.gdal in=river2 out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/vic/river2.tif
r.out.gdal in=so2 out=/Users/jschapMac/Desktop/Tuolumne/Tuolumne6/Threshold/vic/so2.tif

# --------------------------------------------------------------

# --------------------------------------------------------------
# GRASS

####### Within GRASS GIS, we create a new region with Lat/Long projection and work within that region to create the necessary files.

# Tuolumne.merc is the name of the LOCATION used for the previous GRASS commands. Now we create a new LOCATION called Tuolumne.latlon for creating the VIC routing model input files.

cd /Users/jschapMac/Desktop/Tuolumne/Tuolumne7/GIS

g.proj epsg=4326 location=Tuolumne.latlon
g.mapset mapset=PERMANENT location=Tuolumne.latlon

# Must run this from the GRASS command line, not the console.
g.region `r.proj location=Tuolumne.merc mapset=Tuolumne400 in=hydrosheds.elev -g`
# g.region n=39.5 s=32.25 e=-91.875 w=-106.75 res=0:00:30

r.proj mapset=Tuolumne400 location=Tuolumne.merc in=hydrosheds.felev method=bilinear
r.proj mapset=Tuolumne400 location=Tuolumne.merc in=basin
r.mapcalc exp='basin1=if(isnull(basin),0,1)'

# Enter the VIC modeling resolution here.
g.region res=0.0625 -a

# Choose a suitable threshold
# r.mapcalc exp="flowacc_abs=abs(hydrosheds.flowacc)"
# r.threshold flowacc_abs

r.watershed --overwrite -s elev=hydrosheds.felev accum=hydrosheds.flowacc drain=hydrosheds.flowdir stream=hydrosheds.river threshold=5

r.resamp.stats in=basin1 out=fract

# This line requires a "rules" file for converting between flow direction conventions
r.reclass in=hydrosheds.flowdir out=flowdir rules=/Users/jschapMac/Documents/Codes/VICMATLAB/PythonScripts/flowdir.rules
# !Potential source of error. Since I don't have the flowdir.rules files, I used modfdir.R to make the conversion, then reloaded the modified flow direction as flowdir with r.in.gdal.
# Also, there may be an error in the .rules file. Check this experimentally.

r.mapcalc --o exp='flowdir=flowdir'
r.null map=flowdir setnull=0
r.null map=fract null=0
r.out.gdal --overwrite in=fract out=tuolumne.fract format=AAIGrid nodata=0
r.out.gdal --overwrite in=flowdir out=tuolumne.flowdir format=AAIGrid nodata=0
r.out.gdal --overwrite in=hydrosheds.flowacc out=vic_flowacc.asc format=AAIGrid nodata=0
r.out.gdal --overwrite in=hydrosheds.river out=vic_river.asc format=AAIGrid nodata=0

# View the result with plotflowdir.m to check that the flow direction file looks good.

# Change lat/lon and resolution in the command below as necessary. See comment where the tuolumne.stations file generated, above.
awk 'BEGIN{OFS="|"}/^1/{print(-121.1875+0.0625*$3,37.5625+0.0625*$4,$2)}' tuolumne.stations | v.in.ascii in=- out=inflows col='x real,y real,name text'

v.out.ascii in=inflows out=inflows.asc

# Choose an appropriate snapping distance so the stations are placed on the "river channel"
r.stream.snap --o in=inflows out=stations accum=hydrosheds.flowacc stream=hydrosheds.river radius=3

v.out.ascii in=stations out=stations.asc

# Updates the station location file with the locations of the stations after they have been snapped to nearby river channel pixels
v.to.rast --o in=stations out=inflows use=cat
r.stats --o in=inflows -n -x | sort -k3 -n | awk '{printf("1\tTUO%03d\t%d\t%d\t-9999\nNONE\n",$3,$1,11-$2+1)}' > tuolumne.stations
# In the above line, be sure to change the number XX-$2+1 to nrows in your 1/16 degree flow direction file

v.db.addtable map=stations
v.db.addcolumn map=stations col='flowacc real'
v.what.rast map=stations raster=hydrosheds.flowacc column=flowacc

# Prints the flow accumulation values at the station locations (after snapping) to a file
# Scaling factor is situational: Kostas used 157.828125 for Arkansas simulation
v.db.select stations sep=" " | sort -k1 -n | awk 'function abs(x){return ((x < 0.0) ? -x : x)}NR>1{print(abs($2)*39.45703)}' > stations.flowacc

# -------------------------------------------------------------------
# VIC

# Run vicinputworkflow to prepare input files for VIC

# Use the bounding box of the tuolumne.flowdir file as the mask. We are running VIC for the rectangular extent, and using the fract file to choose which cells to put in the routing model.

# When loading the .asc extent file, it may not have a nodata value unless you specify one. See gdal_translate AAIGrid lines of code above.

# Run VIC 4.2 to (temporally) disaggregate the Livneh met. forcing data

# Run VIC 5 Classic (or image) to produce flux files with runoff at each 1/16 degree grid cell

# -------------------------------------------------------------------
# VIC ROUTING MODEL

# Copy flowdir, fract, and stations files to Rout_Inputs
# Shorten filenames, if necessary for Fortran (72 chars)
# Create routing model parameter file
# Remove extensions and headers from VIC fluxes.*.txt output files
# Run the routing model

#--------------------------------------------------------------------
# PYTHON

####### Finally, the BDY file is written using the VIC routing model's output.
import numpy as np

# Change the prefix AA{0:03d} to match the prefix used for boundary inflows in the bci file. The key needs to match.
stations = ["TUO{0:03d}".format(s+1) for s in range(len(inflows))]

a0 = np.array([abs(flowacc[i[0],i[1]]) for i in inflows])
a = np.loadtxt("stations.flowacc")

# This is a weighting factor to account for the station location snap AND?OR for the difference between flow acc calculated using the high resolution DEM vs. the low resolution DEM
area_mult = a0 / a

lisflood.write_bdy("/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/LF_Inputs/tuolumne.bdy", "/Users/jschapMac/Desktop/Tuolumne/Tuolumne7/Rout_Results", stations, area_mult)

# Create LISFLOOD parameter file
# Make sure the BCI file has the correct identifiers
# Run LISFLOOD using lisflood_double_all_channel_MacOS executable
