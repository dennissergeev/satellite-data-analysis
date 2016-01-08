# -*- coding: utf-8 -*-
from mpl_toolkits.basemap import Basemap
import numpy as np
from osgeo import gdal
import re

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
    lon0, lat0 = get_lonlat_from_proj_str(gtif.GetProjection())
    sat_m = Basemap(projection='cass',lon_0=lon0,lat_0=lat0, width=nx, height=ny)
    #Convert image coordinates from meters to lon/lat
    lons, lats = sat_m(xx,yy,inverse=True)

    return lons, lats, arr


def get_lonlat_from_proj_str(proj_str, lon_key='central_meridian',
                                       lat_key='latitude_of_origin'):
    lon0 = _get_numvalue_by_str(proj_str, lon_key)
    lat0 = _get_numvalue_by_str(proj_str, lat_key)

    return lon0, lat0


def _get_numvalue_by_str(s, key):
    ibeg = s.find(key)
    if s[ibeg-2] != '[':
        warnings.warn('Uncertain results!')

    iend = s[ibeg:].find(']')
    str_till_next_bracket = s[ibeg:ibeg+iend]
    found = re.findall(r"[+-]?\d+(?:\.\d+)?", str_till_next_bracket)
    if len(found) > 0:
        return float(found[0])
    else:
        warnings.warn('Not found any numeric values')
