# Functions used by the lisflood wrapper
#
# Modified August 24, 2017 (JRS)

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import gdal

def has_inflow(river, flowdir, i, j):
    """Returns true if a river pixel flows into pixel (i, j)"""
    flows = False
    kl = [(k, l) for k in range(-1, 2) for l in range(-1, 2) if not (k == 0 and l == 0)]
    fdir = {kl[fi]: f for fi, f in enumerate([7, 6, 5, 8, 4, 1, 2, 3])}
    for k, l in kl:
        if river[i+k, j+l] and flowdir[i+k, j+l] == fdir[(k, l)]:
            flows = True
    return flows


def find_boundary_points(river, flowdir):
    """Find inflow and outflow points along a river network defined
    by the river raster file."""
    bnd = []
    nr, nc = river.shape
    priver = np.zeros((nr+2, nc+2), dtype='bool')
    priver[1:-1, 1:-1] = river
    pflowdir = np.zeros((nr+2, nc+2), dtype='int')
    pflowdir[1:-1, 1:-1] = flowdir
    for i in range(nr):
        for j in range(nc):
            if river[i, j]:
                if not has_inflow(priver, pflowdir, i+1, j+1):
                    bnd.append((i, j))
                # neighbors = [river[i+k, j+l] for k in range(-1, 2) for l in range(-1, 2) if i+k >= 0 and i+k < nr and j+l >= 0 and j+l < nc]
                # if len(np.where(neighbors)[0]) < 3:
                #     bnd.append((i, j))
    return bnd


def calc_chainage(river, flowdir, flowacc, bndpts, res):
    """Generate raster files of river channel chainage (i.e. downstream distance)
    and labels. Inputs include:
    river: boolean raster defining the river network
    flowdir: flow direction raster
    flowacc: flow accumulation raster
    bndpts: array of inflow/outflow coordinate tuples
    res: pixel resolution in meters
    """
    labels = np.zeros(river.shape)
    chainage = np.zeros(river.shape)
    nr, nc = river.shape
    riv = 0
    for i in np.argsort([flowacc[b[0], b[1]] for b in bndpts]):
        if labels[bndpts[i][0], bndpts[i][1]] == 0:
            riv += 1
            dist = 0.0
            ci, cj = bndpts[i]
            while labels[ci, cj] == 0 and river[ci, cj]:
                labels[ci, cj] = riv
                chainage[ci, cj] = dist
                dist += res
                downstream = [(ci-1, cj+1), (ci-1, cj), (ci-1, cj-1), (ci, cj-1), (ci+1, cj-1), (ci+1, cj), (ci+1, cj+1), (ci, cj+1)]
                if flowdir[ci, cj] > 0:
                    ci, cj = downstream[abs(flowdir[ci, cj])-1]
    return labels, chainage


def smooth_bank_heights(elev, labels, chainage):
    """Smooth river channel bank heights as a means of hydrologic conditioning"""
    selev = np.copy(elev)
    nrivers = int(np.amax(labels))
    for i in range(nrivers):
        selev[labels == i+1] = lowess(elev[labels == i+1], chainage[labels == i+1], return_sorted=False, frac=0.1)
    return selev


def burn_channel(elev, river, depth):
    """Burns river channel on DEM elevation using depth values"""
    belev = elev - depth
    return belev


def read_raster(filename):
    """Read data from raster file filename."""
    f = gdal.Open(filename)
    data = f.ReadAsArray()
    f = None
    return data


def write_raster(data, filename, basefile):
    """Writes data into Geotiff raster using extent and
    projection from basefile raster"""
    nr, nc = data.shape
    f = gdal.Open(basefile)
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(filename, nc, nr, 1, gdal.GDT_Float32)
    xul, xres, xrot, yul, yrot, yres = f.GetGeoTransform()
    dataset.SetGeoTransform((xul, xres, xrot, yul, yrot, yres))
    dataset.SetProjection(f.GetProjection())
    dataset.GetRasterBand(1).WriteArray(data)
    dataset.FlushCache()
    f = None
    dataset = None


def identify_inflows(river, chainage, labels, flowacc, threshold=10000):
    """Extract upstream and identify lateral inflow points"""
    # find upstream inflows as river points at zero chainage
    i, j = np.where(np.logical_and(river, chainage == 0))
    inflows = list(zip(i, j))
    # find potential lateral inflows
    nrivers = int(np.amax(labels))
    for riv in range(nrivers):
        c = chainage[labels == riv+1]
        a = abs(flowacc[labels == riv+1])
        # find river locations where flow accumulation increases above threshold
        breakpoints = np.where(np.diff(a[np.argsort(c)]) > threshold)
        if len(breakpoints) > 0:
            for b in breakpoints[0]:
                i = np.where(labels == riv+1)[0][np.argsort(c)][b+1]
                j = np.where(labels == riv+1)[1][np.argsort(c)][b+1]
                # find pixel with largest flow accumulation that doesn't belong to current river
                neighbors = [(i+k, j+l) for k in range(-1, 2) for l in range(-1, 2) if labels[i+k, j+l] != riv+1]
                if len(neighbors) > 0:
                    ni, nj = neighbors[np.argsort([abs(flowacc[n[0], n[1]]) for n in neighbors])[-1]]
                    # if identified point belongs to other river it is outlet and is not lateral inflow
                    if labels[ni, nj] == 0:
                        inflows.append((ni, nj))
    return inflows


def write_bci(bcifilename, inflows, xul, yul, xres, yres, nrows, ncols, uniform_slope=0.0001):
    """Write BCI file for LISFLOOD"""
    with open(bcifilename, 'w') as fout:
        fout.write("N {0:.1f} {1:.1f} FREE {2}\n".format(xul, xul+ncols*xres, uniform_slope))
        fout.write("S {0:.1f} {1:.1f} FREE {2}\n".format(xul, xul+ncols*xres, uniform_slope))
        fout.write("E {0:.1f} {1:.1f} FREE {2}\n".format(yul, yul+nrows*yres, uniform_slope))
        fout.write("W {0:.1f} {1:.1f} FREE {2}\n".format(yul, yul+nrows*yres, uniform_slope))
        for pi, pt in enumerate(inflows):
            x = xul + xres * pt[1]
            y = yul + yres * pt[0]
            fout.write("P {0:.1f} {1:.1f} QVAR TUO{2:03d}\n".format(x, y, pi+1))


def write_gauges(x, y, flowdir, widthfile, gaugefile):
    """Writes virtual gauge file"""
    f = gdal.Open(widthfile)
    xul, xres, _, yul, _, yres = f.GetGeoTransform()
    widths = f.ReadAsArray()
    nr, nc = widths.shape
    flowdirs = {8: "E", 2: "N", 6: "S", 4: "W"}
    with open(gaugefile, 'w') as fout:
        fout.write("{0}\n".format(len(x)))
        for c in range(len(x)):
            i = int((y[c] - yul) / yres)
            j = int((x[c] - xul) / xres)
            fout.write("{0:.2f} {1:.2f} {2} {3:.2f}\n".format(x[c], y[c], flowdirs[flowdir[i, j]], widths[i, j]))


def write_bdy(outfile, modelpath, stations, area_f=None, resolution=1000.0):
    """Reads flow files from VIC and write LISFLOOD BDY files"""
    with open(outfile, 'w') as fout:
        fout.write("QTBDY\n")
        for i in range(len(stations)):
            x = np.loadtxt("{0}/{1}.day".format(modelpath, stations[i]))
            fout.write("{0}\n{1}\tseconds\n".format(stations[i], len(x)))
            for t in range(len(x)):
                if area_f is None:
                    af = 1.0
                else:
                    af = area_f[i]
                # Note: conversion factor from cfs to m^3/s
                q = af * x[t, 3] * 0.0283168 / resolution
                fout.write("{0:.9f}\t{1:.0f}\n".format(q, t*86400))


def main():
    flowacc = read_raster("../data/hydrosheds/arkansas_acc.tif")
    widths = read_raster("../data/hydrosheds/arkansas_widths.tif")
    river = widths > 0
    bndpts = find_boundary_points(river)
    labels, chainage = calc_chainage(river, abs(flowacc), bndpts, 1000)


if __name__ == '__main__':
    main()
