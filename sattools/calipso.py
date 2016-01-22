# -*- coding: utf-8 -*-
"""
Functions to read and visualise CALIPSO satellite data
"""
from __future__ import division, print_function
import h5py
import numpy as np

from .utils import calipso_time2dt


def geodata(h5name, varnames=['Longitude', 'Latitude', 'Height', 'Profile_UTC_Time'], proftime2datetime=True, return_list=False):
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


def read_data(h5name, data_field='Total_Attenuated_Backscatter_532', fillmask=False, return_units=True):
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
