# Wrapper for creating LISFLOOD input files
#
# Run in the lisflood conda environment (contains gdal)
#
# Modified August 24, 2017 (JRS)

####### The next steps involve preparing the actual LISFLOOD-FP input files. Given the uncertainty in the DEM, we will smooth the river bank heights to avoid any numerical instabilities during simulation. This is accomplished by first identifying the upstream and downstream boundary points in the domain, labeling the river channels and calculating each one's chainage (i.e. distance downstream).
import sys
sys.path.append("/Users/jschapMac/Desktop/Tuolumne/TuoSub")
import lisflood
import os
os.chdir("/Users/jschapMac/Desktop/Tuolumne/TuoSub/GIS")

flowdir = lisflood.read_raster("fdir.tif")
flowacc = lisflood.read_raster("facc.tif")
widths = lisflood.read_raster("rwidth.tif")
river = widths > 0

##### Original code
bndpts = lisflood.find_boundary_points(river, flowdir)
#####

##### New code
# Find boundary points in Matlab and read the results here:
#import csv
#with open("bndpts.csv") as f:
#    bnd = csv.reader(f)
#    bnd2 = list(bnd)
#f.close()

#bnd3 = []]
#for i in range(0, len(bnd2)):
    bnd3.append(map(int, bnd2[i]))

#bndpts = bnd3
#####

labels, chainage = lisflood.calc_chainage(river, flowdir, abs(flowacc), bndpts, 1000)

lisflood.write_raster(labels, "labels.tif", "dem_merc_clipped.tif")
lisflood.write_raster(chainage, "chainage.tif", "dem_merc_clipped.tif")

####### river, bndpts, labels, chainage



####### The smoothing of the bank heights is done by using a LOWESS local regression, and the channel is burned in to the DEM by subtracting the depth raster.
elev = lisflood.read_raster("dem_merc_clipped.tif")
selev = lisflood.smooth_bank_heights(elev, labels, chainage)
lisflood.write_raster(selev, "TuoSub_selev.tif", "dem_merc_clipped.tif")
depths = lisflood.read_raster("rdepth.tif")
depths[depths < 0] = 0.0
belev = selev - depths
lisflood.write_raster(belev, "TuoSub_bed.tif", "dem_merc_clipped.tif")
y0 = 0.8 * depths
y0[widths <= 0] = 0.0
lisflood.write_raster(y0, "TuoSub_initwd.tif", "dem_merc_clipped.tif")
####### OUT: smoothed riverbank(arkansas_selev.tif), riverbed(arkansas_bad.tif), initial depth(80% bank full, arkansas_initwd.tif)



####### Next we identify the upstream and lateral inflow points, and generate the BCI file.
import gdal

##### Original code
flow_thres = 100 # minimum upstream contributing area
inflows = lisflood.identify_inflows(river, chainage, labels, abs(flowacc), flow_thres)
#####

##### New code
# Calculate inflows in Matlab and load results here
#####

f = gdal.Open("TuoSub_selev.tif")
xul, xres, _, yul, _, yres = f.GetGeoTransform()
f = None
nrows, ncols = elev.shape
lisflood.write_bci("TuoSub.bci", inflows, xul, yul, xres, yres, nrows, ncols)
####### OUT: BCI file that include boundary points (arkansas.bci)



####### Then we generate the DEM and sub-grid channel width files.
cd /Users/jschapMac/Desktop/TuoSub/GIS
# These commands require the GDAL version (2.1) in my lisflood conda package, or maybe one of the other utilities installed there.
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 TuoSub_selev.tif /Users/jschapMac/Desktop/TuoSub/LF/TuoSub.dem
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 rwidth.tif /Users/jschapMac/Desktop/TuoSub/LF/TuoSub.width
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 TuoSub_bed.tif /Users/jschapMac/Desktop/TuoSub/LF/TuoSub.bed
gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 -co DECIMAL_PRECISION=2 TuoSub_initwd.tif /Users/jschapMac/Desktop/TuoSub/LF/TuoSub.initwd
####### dem, width, bed, initwd

# I have managed to run the code up to here -JRS 8/25/2017

####### If we need LISFLOOD-FP to produce a time series of river flow at specific locations, we need to generate a virtual gauge file.
x = [-10646496, -10281401]
y = [4251469, 4111473]
lisflood.write_gauges(x, y, flowdir, "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_widths.tif", "/Users/kandread/Work/Swot/input/arkansas.gauge")
####### OUT: gage file (arkansas.gauge)



####### The final step is generating the BDY file for LISFLOOD-FP, which contains the streamflow at the inflow points (in m2/s). For this simulation we will use the output of the VIC routing model forced by the NLDAS-2 VIC model output. We need to prepare the input files for the routing model, describing the flow direction, flow fraction and station locations. The station file is created from the inflows identified
from pyproj import Proj
wmerc = Proj("+init=EPSG:3395")
fout = open("arkansas.stations", 'w')
for sta, xy in enumerate(inflows):
    x, y = wmerc(xul+xres*xy[1], yul+yres*xy[0], inverse=True)
    xi = int((x+106.625)/0.125) + 1
    yi = int((y-33.375)/0.125) + 1
    fout.write("1\tAR{0:03d}\t{1}\t{2}\t-9999\nNONE\n".format(sta, xi, yi))
fout.close()



####### Within GRASS GIS, we create a new region with Lat/Long projection and work within that region to create the necessary files.
g.proj epsg=4326 location=arkansas.latlon
g.mapset mapset=PERMANENT location=arkansas.latlon
g.region `r.proj location=arkansas mapset=arkansas in=hydrosheds.elev -g`
#g.region n=39.5 s=32.25 e=-91.875 w=-106.75 res=0:00:30
r.proj mapset=arkansas location=arkansas in=hydrosheds.felev method=bilinear
r.proj mapset=arkansas location=arkansas in=basin
r.mapcalc exp='basin1=if(isnull(basin),0,1)'
g.region res=0.125 -a
r.watershed -s elev=hydrosheds.felev accum=hydrosheds.flowacc drain=hydrosheds.flowdir stream=hydrosheds.river threshold=15
r.resamp.stats in=basin1 out=fract
r.reclass in=hydrosheds.flowdir out=flowdir rules=flowdir.rules
r.mapcalc --o exp='flowdir=flowdir'
r.null map=flowdir setnull=0
r.null map=fract null=0
r.out.gdal in=fract out=arkansas.fract format=AAIGrid nodata=0
r.out.gdal in=flowdir out=arkansas.flowdir format=AAIGrid nodata=0
awk 'BEGIN{OFS="|"}/^1/{print(-106.625+0.125*$3,33.375+0.125*$4,$2)}' arkansas.stations | v.in.ascii in=- out=inflows col='x real,y real,name text'
r.stream.snap in=inflows out=stations accum=hydrosheds.flowacc stream=hydrosheds.river radius=3
v.to.rast in=stations out=inflows use=cat
r.stats in=inflows -n -x | sort -k3 -n | awk '{printf("1\tAR%03d\t%d\t%d\t-9999\nNONE\n",$3,$1,49-$2+1)}' > arkansas.stations
v.db.addtable map=stations
v.db.addcolumn map=stations col='flowacc real'
v.what.rast map=stations raster=hydrosheds.flowacc column=flowacc
v.db.select stations sep=" " | sort -k1 -n | awk 'function abs(x){return ((x < 0.0) ? -x : x)}NR>1{print(abs($2)*157.828125)}' > stations.flowacc


####### Finally, the BDY file is written using the VIC routing model's output.
import numpy as np
stations = ["AR{0:03d}".format(s+1) for s in range(len(inflows))]
a0 = np.array([abs(flowacc[i[0],i[1]]) for i in inflows])
a = np.loadtxt("stations.flowacc")
area_mult = a0 / a
lisflood.write_bdy("/Users/kandread/Work/Swot/input/arkansas.bdy", "/Volumes/External2/nldas2", stations, area_mult)
