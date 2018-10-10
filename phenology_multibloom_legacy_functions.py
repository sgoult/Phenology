#!/usr/bin/env python
import os
import numpy
import argparse
import netCDF4 as nc
import shutil
import scipy.ndimage as ndimage
import scipy.signal as signal
from scipy.interpolate import interp2d
from mpl_toolkits.basemap import Basemap

def regrid_sst_data(sst_file, chl_lat, chl_lon):
    """
    Takes an SST dataset, extracts the data and regrids it to a chlorophyll lat/lon grid

    direct version of extract_SST.jnl and regrid_sstc_4km.jnl
    """
    sst_ds = nc.Dataset(sst_file)
    lon_variable = [x for x in sst_ds.variables if "lon" in x.lower()][0]
    lat_variable = [x for x in sst_ds.variables if "lat" in x.lower()][0]
    time_variable = [x for x in sst_ds.variables if "time" in x.lower()][0]
    sst_lon = sst_ds.variables[lon_variable][:]
    sst_lat = sst_ds.variables[lat_variable][:]
    lon_regrid, lat_regrid = np.mgrid[0:len(sst_lat), 0:len(sst_lon)]
    sst_time = sst_ds.variables[time_variable][:]
    sst_variable = [x for x in sst_ds.variables if "sst" in x.lower()][0]
    sst_array = sst_ds.variables[sst_variable][:]
    sst_regridded = np.zeros((len(sst_time), len(chl_lat), len(chl_lon)))
    #this bit does not work, duped by interpolation rather than gridding
    for array in range(len(sst_time)):
        sub_sst_array = sst_array[array,:,:]
        X1 = lon_regrid[~sub_sst_array.mask]
        Y1 = lat_regrid[~sub_sst_array.mask]
        Z1 = sub_sst_array[~sub_sst_array.mask]
        f = interp2d(X1, Y1, Z1)
        sst_regridded[array,:,:] = f(chl_lon, chl_lat)
    #if it's over 100 assume its kelvin, since the sea would be boiling at this temp
    if numpy.all(sst_regridded > 100):
        sst_regridded = sst_regridded - 273.15
    return sst_regridded

def regrid_sst_on_time(sst_array, days=8):
    """
    if its a finer resolution than 8 day, regrid, if not then I guess regrid anyway
    """
    #TODO
    return True

def get_chl_array(chl_files):
    chl_arrays = []
    times = []
    for f in chl_files:
        chl_ds = nc.Dataset(f)
        lon_variable = [x for x in chl_ds.variables if "lon" in x.lower()][0]
        lat_variable = [x for x in chl_ds.variables if "lat" in x.lower()][0]
        time_variable = [x for x in chl_ds.variables if "time" in x.lower()][0]
        times.append(chl_ds.variables[time_variable])
        chl_variable = [x for x in chl_ds.variables if "chl" in x.lower()][0]
        chl_arrays.append(chl_ds.variables[chl_variable])
    if len(chl_arrays) == 1:
        chl_array = chl_arrays[0]
        times = times[0]
    elif not chl_arrays:
        raise "failed"
    else:
        chl_array = numpy.dstack(chl_arrays)
    
    return chl_array, times



def cumsum_der_medsea_OCCCI_8d():
    pyferret_runner(jnl_files["cumsum_der_MedSea3.jnl"])


def find_multi_peaks(array_like):
    peaks = signal.argrelmax(numpy.array(smooth_signal), numpy.greater)
    troughs = signal.argrmin(numpy.array(smooth_signal), numpy.less)


def create_phenology_two_blooms():
    """
    read the 
    """
    fix_value = -1.e+34

    ids = nc.Dataset("CCI_ALL-v3.0-8DAY-Sep97-Dec14_9kmba.nc")
    chl =  ids.variables["CHL"][:]
    #this doesn't actually implement the logic in bmed5_MedSea_OCCCI_8d.f, needs addressing
    med5 =  numpy.median(chl,axis = 0)*1.05
    #spit out into netcdf named bmed5_MedSea_OCCCI_8d.nc
    #is there a way this could be numpified
    cumsum_der_medsea_OCCCI_8d()

    


