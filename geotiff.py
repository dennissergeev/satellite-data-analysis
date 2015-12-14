# -*- coding: utf-8 -*-
from mpl_toolkits.basemap import Basemap
import numpy as np
from osgeo import gdal

def read_gtiff(fname):
    """
    Read GeoTiFF image and return coordinates and data array

    TODO: get projection and centre lon/lat automatically
    """
    # Read GeoTiff image using GDAL
    gtif = gdal.Open(fname)
    arr = gtif.ReadAsArray()
    trans = gtif.GetGeoTransform()
    
    # Calculate image dimensions 
    nx = gtif.RasterXSize*trans[1]
    ny = gtif.RasterYSize*trans[5]
    x = np.arange(0, gtif.RasterXSize*trans[1], trans[1])
    y = np.arange(0, gtif.RasterXSize*trans[5], trans[5])
    xx, yy = np.meshgrid(x, y)
    
    # Set a map using the Cassini-Soldner projection
    sat_m = Basemap(projection='cass',lon_0=10.,lat_0=74, width=nx, height=ny)
    #Convert image coordinates from meters to lon/lat
    lons, lats = sat_m(xx,yy,inverse=True)

    return lons, lats, arr
