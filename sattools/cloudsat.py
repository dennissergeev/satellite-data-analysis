# -*- coding: utf-8 -*-
"""
Functions to read and visualise CloudSat data

Some parts include modified code from the following open-source projects:
    * ccplot (http://ccplot.org) Copyright (c) 2009-2015 Peter Kuma
    * https://github.com/a301-teaching/classcode
"""
from __future__ import division, print_function
import datetime
import h5py
import numpy as np

from .utils import cc_interp2d


def cloudsat_geodata(h5name, 
                     varnames=['Longitude','Latitude', 'Height', 'Profile_time','DEM_elevation'],
                     proftime2datetime=True, return_list=False):
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
    """    
    Convert cloud scenario codes to one of the 8 classes

    source: http://disc.sci.gsfc.nasa.gov/measures/documentation/README.AIRS_CloudSat.pdf
    """
    with h5py.File(h5name, 'r') as f:
        cldclass = f['2B-CLDCLASS/Data Fields/cloud_scenario'].value
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
