# -*- coding: utf-8 -*-
"""
@author: D. E. Sergeev
"""
from __future__ import division, print_function
import datetime
import dateutil.tz as tz
import h5py
import matplotlib as mpl
import numpy as np
import os
import pytz

def calipso_geodata(h5name, varnames=['Longitude', 'Latitude', 'Height', 'Profile_UTC_Time'], proftime2datetime=True, return_list=False):
    with h5py.File(h5name,'r') as f:
        var_dict = {}
        for ivar in varnames:
            if 'height' in ivar.lower():
                var_dict[ivar] = f['metadata'][0][29]
            else:
                var_dict[ivar] = f[ivar].value

        if 'Profile_UTC_Time' in var_dict:
            var_dict['Profile_UTC_Time'] = np.array([i[0] for i in var_dict['Profile_UTC_Time']])
            #datetime.datetime.strptime(f['metadata'][0][1], '%Y-%m-%dT%H:%M:%S.%fZ')
            if proftime2datetime:
                var_dict['Profile_UTC_Time'] = np.array([calipso_time2dt(t, tzinfo=None) for t in var_dict['Profile_UTC_Time']])

        if return_list:
            return [var_dict[i] for i in varnames]
        else:
            return var_dict


def calipso_read_data(h5name, data_field='Total_Attenuated_Backscatter_532', fillmask=False, return_units=True):
    with h5py.File(h5name,'r') as f:
        data_raw = f[data_field]
        data_dtype = ''.join(data_raw.attrs['format'].lower().split('_'))
        data = data_raw.value.astype(data_dtype)
        try:
            data = np.ma.masked_equal(data, data_raw.attrs['fillvalue'])
        except KeyError:
            pass
        try:
            valid_range = [np.float32(i) for i in data_raw.attrs['valid_range'].split('...')]
            data = np.ma.masked_outside(data, *valid_range)
        except KeyError:
            pass

        if fillmask:
            data = data.filled(np.nan)
        if return_units:
            return data, data_raw.attrs['units']
        else:
            return data

def cloudsat_geodata(h5name, varnames=['Longitude','Latitude', 'Height', 'Profile_time','DEM_elevation'], proftime2datetime=True, return_list=False):
    with h5py.File(h5name,'r') as f:
        for ikey in f.keys():
            if 'Data Fields' in f[ikey]:
                topkey = ikey
                break

        var_dict = {}
        for ivar in varnames:
            fld_dtype = f[topkey]['Geolocation Fields'][ivar].value[0][0].dtype
            var_dict[ivar] = f[topkey]['Geolocation Fields'][ivar].value.astype(fld_dtype)

        if 'Profile_time' in var_dict:
            start_time = f[topkey]['Swath Attributes']['start_time'][0][0].decode('UTF-8')
            start_time = datetime.datetime.strptime(start_time, '%Y%m%d%H%M%S')
            if proftime2datetime:
                var_dict['Profile_time'] = np.array([start_time + datetime.timedelta(seconds=float(i)) for i in var_dict['Profile_time']])

        if 'DEM_elevation' in var_dict:
            var_dict['DEM_elevation'][var_dict['DEM_elevation']<0] = 0.

        if return_list:
            return [var_dict[i] for i in varnames]
        else:
            return var_dict


def cloudsat_read_data(h5name, data_field='Radar_Reflectivity', limits=None, fillmask=False):
    with h5py.File(h5name, 'r') as f:
        for ikey in f.keys():
            if 'Data Fields' in f[ikey]:
                topkey = ikey
                break

        data = f[topkey]['Data Fields'][data_field].value.astype(np.float32)

        missval = f[topkey]['Swath Attributes'][data_field+'.missing'].value[0][0]
        fillval = f[topkey]['Swath Attributes']['_FV_'+data_field].value[0][0]

        for imask in [missval, fillval]:
            data = np.ma.masked_equal(data, imask)

        if limits is not None:
            lim_mask = np.logical_or(data < limits[0], data > limits[1])
            data = np.ma.masked_equal(data, lim_mask)

        if fillmask:
            data = data.filled(np.nan)

        data_factor = f[topkey]['Swath Attributes'][data_field+'.factor'][0][0]
        data_offset = f[topkey]['Swath Attributes'][data_field+'.offset'][0][0]
        data = (data - data_offset) / data_factor

    return data

def cloudsat_read_cldclass(h5name):
    with h5py.File(h5name, 'r') as f:
        cldclass = f['2B-CLDCLASS/Data Fields/cloud_scenario'].value
        # Convert cloud scenario codes to one of the 8 classes
        # source: http://disc.sci.gsfc.nasa.gov/measures/documentation/README.AIRS_CloudSat.pdf
        cldclass = np.right_shift(np.array(cldclass), 1) & int('1111', 2)
        cloudcodes = {
                      '0': 'clear',
                      '1': 'Ci',
                      '2': 'As',
                      '3': 'Ac',
                      '4': 'St',
                      '5': 'Sc',
                      '6': 'Cu',
                      '7': 'Ns',
                      '8': 'DC'
                     }

    return cldclass, cloudcodes

def cc_interp2d(data, X, Z, x1, x2, nx, z1, z2, nz, use_numba=True):
    if use_numba:
        try:
            import numba as nb
        except ImportError:
            print('Unsuccessful numba import, using pure Python')
            use_numba = False

    if use_numba:
        res = nb.jit()(_interp2d)(data, X, Z, x1, x2, nx, z1, z2, nz)
    else:
        res = _interp2d(data, X, Z, x1, x2, nx, z1, z2, nz)
    return res

def _interp2d(data, X, Z, x1, x2, nx, z1, z2, nz):
    """Interpolate 2D data with coordinates given by 1D and 2D arrays.

    data is a two-dimensional array of data to be interpolated.
    X and Z are one- and two-dimensional arrays, giving coordinates
    of data points along the first and second axis, resp.

    data, X and Z are expected to be C-contiguous float32 numpy arrays
    with no mask and no transformation (such as transposition) applied.
    """

    xs = (x2 - x1)/nx
    zs = (z2 - z1)/nz
    w = data.shape[0]
    h = data.shape[1]

    out = np.zeros((nx, nz), dtype=np.float32)
    q = np.zeros((nx, nz), dtype=np.int32)

    for i in range(w):
        n1 = ((X[i-1] + X[i])/2 - x1)/xs if i-1 >= 0 else -1
        n2 = ((X[i+1] + X[i])/2 - x1)/xs if i+1 < w else nx
        if n2 - n1 < 1: n1 = n2 = (X[i] - x1)/xs

        for j in range(h):
            m1 = ((Z[i,j-1] + Z[i,j])/2 - z1)/zs if j-1 >= 0 else -1
            m2 = ((Z[i,j+1] + Z[i,j])/2 - z1)/zs if j+1 < h else nz
            if m2 - m1 < 1: m1 = m2 = (Z[i,j] - z1)/zs

            for n in range(int(n1+0.5), int(n2+0.5+1)):
                for m in range(int(m1+0.5), int(m2+0.5+1)):
                    if n < 0 or n >= nx: continue
                    if m < 0 or m >= nz: continue
                    if np.isnan(data[i,j]): continue
                    out[n,m] += data[i,j]
                    q[n,m] += 1

    for n in range(nx):
        for m in range(nz):
            if q[n,m] == 0:
                out[n,m] = np.nan
            else:
                out[n,m] /= q[n,m]
    return out

def calipso_time2dt(time, tzinfo=pytz.utc):
    """Convert float in format yymmdd.ffffffff to datetime."""
    d = int(time % 100)
    m = int((time-d) % 10000)
    y = int(time-m-d)
    return datetime.datetime(2000 + y//10000, m//100, d, tzinfo=tzinfo) + datetime.timedelta(time % 1)


def get_cc_cmap(satname='cloudsat', cmap_dir=os.path.join(os.environ.get('HOME'),'anaconda3/share/ccplot/cmap/')):
    if 'cloudsat' in satname.lower():
        cmap_file = os.path.join(cmap_dir, 'cloudsat-reflectivity.cmap')
    elif 'calipso' in satname.lower():
        cmap_file = os.path.join(cmap_dir, 'calipso-backscatter.cmap')
    else:
        raise ValueError('Unrecognized satellite name')
    cmap = _cmap(cmap_file)
    cm = mpl.colors.ListedColormap(cmap['colors']/255.0)
    cm.set_under(cmap['under']/255.0)
    cm.set_over(cmap['over']/255.0)
    cm.set_bad(cmap['bad']/255.0)
    norm = mpl.colors.BoundaryNorm(cmap['bounds'], cm.N)
    return dict(cmap=cm, norm=norm)


def _cmap(filename):
    """Load colormap from file. The expected format of the file is:

        BOUNDS
        from1 to1 step1
        from2 to2 step2
        [...]

        TICKS
        from1 to1 step1
        from2 to2 step2
        [...]

        COLORS
        r1 g1 b1
        r2 g2 b2
        [...]

        UNDER_OVER_BAD_COLORS
        ro go bo
        ru gu bu
        rb gb bb

    where fromn, ton, stepn are floating point numbers as would be supplied
    to numpy.arange, and rn, gn, bn are the color components the n-th color
    stripe. Components are expected to be in base 10 format (0-255).
    UNDER_OVER_BAD_COLORS section specifies colors to be used for
    over, under and bad (masked) values, in that order.
    """
    bounds = []
    ticks = []
    colors = []
    special = []
    mode = "COLORS"
    white = [1, 1, 1, 1]

    try:
        with open(filename) as f:
            for n, s in enumerate(f.readlines()):
                s = s.strip()

                # Skip blank lines.
                if len(s) == 0:
                    continue

                if s in ("BOUNDS", "TICKS", "COLORS", "UNDER_OVER_BAD_COLORS"):
                    mode = s
                    continue

                a = s.split()
                if len(a) not in (3, 4):
                    raise ValueError("Invalid number of fields")

                if mode == "BOUNDS":
                    bounds += list(np.arange(float(a[0]), float(a[1]), float(a[2])))
                elif mode == "TICKS":
                    ticks += list(np.arange(float(a[0]), float(a[1]), float(a[2])))
                elif mode == "COLORS":
                    rgba = [int(c) for c in a]
                    if len(rgba) == 3:
                        rgba.append(255)
                    colors.append(rgba)
                elif mode == "UNDER_OVER_BAD_COLORS":
                    rgba = [int(c) for c in a]
                    if len(rgba) == 3:
                        rgba.append(255)
                    special.append(rgba)
    except IOError as e:
        raise e
    except ValueError as e:
        raise ValueError("Error reading `%s' on line %d: %s" %
                         (filename, n+1, e))

    return {
        'colors': np.array(colors),
        'bounds': np.array(bounds),
        'ticks': np.array(ticks),
        'under': np.array(special[0] if len(special) >= 1 else white),
        'over': np.array(special[1] if len(special) >= 2 else white),
        'bad': np.array(special[2] if len(special) >= 3 else white),
    }
