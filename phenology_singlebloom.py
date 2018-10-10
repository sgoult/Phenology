#!/usr/bin/env python
import os
import numpy
import pyferret
import argparse
import netCDF4 as nc
import shutil

from datetime import date, datetime, timedelta

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

jnl_files = ['apply_sbx3weeks.jnl', 'fill_on_time_chl_v2.jnl', 'map_pheno_ANO_indices1bloom_upwelling_civ_v4.jnl', 'map_pheno_indices1bloom_upwelling_civ_v4.jnl', 'save_Sep97_Dec14_nofill_CHLOR.jnl']

jnl_files = {}
for f in os.listdir(os.path.join(__file__, "ferret")):
    jnl_files[f] = os.path.join(__file__, "ferret", f)    

def pyferret_runner(pyferret_file):
    """
    This is a wrapper around pyferret to run a jnl file one line at a time, should help with debugging. Will do very little with the actual file other than trying to run it line by line.

    failures are raised, normally if this happens we should go back and debug the jnl file.
    """
    if os.path.exists(pyferret_file):
        pyferret_file = list(open(pyferret_file))
    else:
        raise Exception("fatal, pyferret command file {} does not exist!".format(pyferret_file))
    
    pyferret.init(enterferret=False)
    for line in pyferret_file:
        err_int, err_msg = pyferret.run(command="line")
        if not err_int == pyferret.FERR_OK:
            raise Exception(err_msg)
    pyferret.stop()
    return True

def fill_chl_on_time(files):
    pyferret_runner(jnl_files["fill_on_time_chl_v2.jnl"])

def apply_sbx3weeks(files):
    pyferret_runner(jnl_files["apply_sbx3weeks.jnl"])

def save_Sep97_Dec14_nofill_CHLOR(files):
    """
    Applies a landmask
    """
    pyferret_runner(jnl_files["save_Sep97_Dec14_nofill_CHLOR.jnl"])

#everything after this definitely works

def get_start_index_and_duration(array_like):
    """
    takes a list of values and searches for the max, then the associated sign changes to indicate bloom start/end
    """
    array_like = list(array_like)
    max_val = max(array_like)
    max_idx = array_like.index(max(array_like[20:-1]))
    if max_idx == 0:
        #print "max is 0"
        return [-1000,-1000,-1000,-1000]
    forward = array_like[max_idx:-1]
    backward = array_like[0:max_idx]
    if not backward or not forward:
        print backward, forward
        return [-1000,-1000,-1000,-1000]
    try:
        end_flip = next(x for x, elem in enumerate(forward) if elem < 0)
    except:
        end_flip = next(x for x, elem in enumerate(forward) if elem <= min(forward))
    end_flip = max_idx + end_flip
    try:
        start_flip = next(x for x, elem in enumerate(list(reversed(backward))) if elem < 0)
    except:
        start_flip = next(x for x, elem in enumerate(list(reversed(backward))) if elem <= min(backward))
    start_flip = max_idx - (start_flip + 1)
    duration = end_flip - start_flip 
    return [start_flip  - 5, max_idx  - 5, end_flip  - 5, duration]

def create_phenology_one_bloom():
    """
    read the 
    """
    fix_value = -1.e+34
    #do sbxweek
    #to convert to work with a list of files replace the logic from here
    ids = nc.Dataset("CCI_ALL-v3.0-8DAY-Sep97-Dec14_9kmba.nc")
    chl =  ids.variables["CHL"][:]
    med5 =  numpy.ma.median(chl,axis = 0)*1.05
    ids.close()
    ids = nc.Dataset("OCCCI_v3-8DAY-97_14_9km_pheno_fav8_sbx3weeks.nc")
    chl = ids.variables["CHL"][:]
    #to here. Should be doable to fill the numpy arrays in python.
    data=[]
    create_chl_netcdf(med5, name="chl.nc")
    for year in range(0, chl.shape[0], 47):
        #create arrays from variable data again
        prev_year = year - 5 if year > 0 else 0
        next_year = year + 47
        print year
        array = chl[prev_year:next_year,:,:,:]
        max_dates = numpy.argmax(array, axis=0)
        min_dates = numpy.argmin(array, axis=0)
        max_vals = numpy.amax(array, axis=0)
        min_vals = numpy.amin(array, axis=0)
        med0_array = array - med5/1.05*1.2
        print med0_array.shape
        start_end_duration_array = numpy.apply_along_axis(get_start_index_and_duration, 0, med0_array)
        print start_end_duration_array.shape
        data.append(start_end_duration_array)
    create_phenology_netcdf(data)


def create_chl_netcdf(data, name="chl.nc"):
    ds = nc.Dataset(name,'w',format='NETCDF4_CLASSIC')
    ds.createDimension('LONGITUDE', data.shape[2])
    ds.createDimension('LATITUDE', data.shape[1])
    ds.createDimension('DEPTH', data.shape[0])
    ds.createDimension('TIME', 1)
    ds.createVariable('LATITUDE', 'float32', dimensions=['LATITUDE'])
    ds.variables['LATITUDE'].setncattr("units", "degrees north") 
    ds.createVariable('LONGITUDE', 'float32', dimensions=['LONGITUDE'])
    ds.variables['LONGITUDE'].setncattr("units", "degrees east") 
    ds.createVariable('DEPTH', 'float32', dimensions=['DEPTH'])
    ds.variables['DEPTH'].setncattr("units", "metres") 
    ds.createVariable('TIME', 'float32', dimensions=['TIME'])
    ds.variables['TIME'].setncattr("units", "weeks") 
    ds.createVariable('CHL', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['CHL'].setncattr("units", "mg m^3") 
    ds.variables['CHL'][:] = data[:]
    ds.close()

def create_phenology_netcdf(data, name="output2.nc"):
    ds = nc.Dataset(name,'w',format='NETCDF4_CLASSIC')
    ds.createDimension('LONGITUDE', data[0].shape[3])
    ds.createDimension('LATITUDE', data[0].shape[2])
    ds.createDimension('DEPTH', data[0].shape[1])
    ds.createDimension('TIME', len(data))
    ds.createVariable('LATITUDE', 'float32', dimensions=['LATITUDE'])
    ds.variables['LATITUDE'].setncattr("units", "degrees north") 
    ds.createVariable('LONGITUDE', 'float32', dimensions=['LONGITUDE'])
    ds.variables['LONGITUDE'].setncattr("units", "degrees east") 
    ds.createVariable('DEPTH', 'float32', dimensions=['DEPTH'])
    ds.variables['DEPTH'].setncattr("units", "metres") 
    ds.createVariable('TIME', 'float32', dimensions=['TIME'])
    ds.variables['TIME'].setncattr("units", "weeks") 
    ds.createVariable('date_start', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_start'].setncattr("units", "weeks") 
    ds.createVariable('date_max', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_start'].setncattr("units", "weeks") 
    ds.createVariable('date_end', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_start'].setncattr("units", "weeks") 
    ds.createVariable('duration', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['duration'].setncattr("units", "weeks") 
    ds.variables['date_start'][:] = numpy.stack([x[0] for x in data])
    ds.variables['date_max'][:] = numpy.stack([x[1] for x in data])
    ds.variables['date_end'][:] = numpy.stack([x[2] for x in data])
    ds.variables['duration'][:] = numpy.stack([x[3] for x in data])
    ds.close()

def save_eighteenyr_eightd(chlor_data_file):
    ids = nc.Dataset("CCI_ALL-v3.0-8DAY-Sep97-Dec15_9kmba_fillphenob.nc")
    ds = nc.Dataset("CCIv3-8DAY-Jan97-Dec14_9km_fillphenob_2.nc",'w',format='NETCDF4_CLASSIC')
    dims = ids.dimensions.keys()
    for dim in dims:
        if "LONGITUDE" in dim:
            dimname = "LONGITUDE"
        elif "LATITUDE" in dim:
            dimname = "LATITUDE"
        else:
            dimname = dim
        ds.createDimension(dimname, len(ids.variables[dim][:]))
    print ds.variables
    variables = ids.variables.keys()
    for var in variables:
        print var
        if "LONGITUDE" in var:
            varname = "LONGITUDE"
        elif "LATITUDE" in var:
            varname = "LATITUDE"
        else:
            varname = var
        dims = ids.variables[var].dimensions
        new_dims = []
        for dim in dims:
            if "LONGITUDE" in dim:
                new_dims.append("LONGITUDE")
            elif "LATITUDE" in dim:
                new_dims.append("LATITUDE")
            else:
                new_dims.append(dim)
        print "var creation"
        ds.createVariable(varname, ids.variables[var].datatype, dimensions=new_dims)
        print "test"
        for attr in ids.variables[var].ncattrs():
            print "attributes"
            ds.variables[varname].setncattr(attr, ids.variables[var].getncattr(attr))  
        print varname, var
        test = ids.variables[var][:]
        print ds.variables[varname]
        ds.variables[varname][:] = test
        print "preattributes"
    ds.close()
    ids.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser
    parser.add_argument("blooms", options=[1,2])
    parser.add_argument("filelist")
    parser.add_argument("file_directory")
    args = parser.parse_args()

    if args.blooms == 1:
        create_phenology_one_bloom()

    if args.blooms == 2:
        print "still wroking on it..."