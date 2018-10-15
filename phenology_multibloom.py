import os
import numpy
import argparse
import netCDF4 as nc
import shutil
import scipy.ndimage as ndimage
import scipy.signal as signal
import tempfile
import glob
from scipy.interpolate import interp2d
from scipy import ndimage
from scipy import signal

output_location = None

def boxcar_smoothing(x, window):
   return x*window

def centered_diff_derivative(array_like):
   """
   calculates the first derivative from centered difference
   """
   x = range(0, array_like.shape[0])
   dx = x[1] - x[0] # if we ever change to non uniform jumps in time use numpy.diff
   cf = numpy.convolve(array_like, [1,-1],'same') / dx
   return cf

def get_start_index_and_duration(array_like,depth=100):
    """
    takes a list of values and searches for the max, then the associated sign changes to indicate bloom start/end
    set depth to change limit of number of identifications
    """
    zero_crossings = numpy.where(numpy.diff(numpy.sign(array_like)))[0]
    #filter out crossings that are too close to each other:
    true_poss = [x for x in zero_crossings if not (x+1) in zero_crossings and not (x-1) in zero_crossings]
    #find out which way we're going
    starts = []
    ends = []
    for index in true_poss:
        forward_index = index + 1 if not (index + 1) > (len(array_like) + 1) else index
        backward_index =  index - 1 if not (index - 1) < 0 else index
        if array_like[forward_index] >= array_like[index] and array_like[backward_index] <= array_like[index]:
            starts.append(index)
        elif array_like[forward_index] <= array_like[index] and array_like[backward_index] >= array_like[index]:
            ends.append(index)
    #find an end for every start
    dates = []
    for start in starts:
        try:
           end = next(x for x in ends if x > start)
           max_idx = numpy.argmax(array_like[start:end]) + start
           dates.append([start,end,end-start,max_idx])
        except Exception as e:
           continue
    for pad in range(len(dates), depth):
        dates.append([None,None,None,None])
    return dates

def match_start_end_to_solar_cycle(array_like, start_end_list, chl_slice, reverse_search):
    """
    Corrects the start and end times based on solar activity
    """
    #remove durations below 2 weeks
    #possibly resort and create new durations based on remaining dates
    #look for sign changes in sst data, indicating high/low light periods
    array_like = numpy.squeeze(array_like)
    zero_crossings = numpy.where(numpy.diff(numpy.sign(array_like)))[0]
    #filter out crossings that are too close to each other:
    true_poss = [x for x in zero_crossings if not (x+1) in zero_crossings and not (x-1) in zero_crossings]
    #find out which way we're going
    highs = []
    lows = []
    for index in true_poss:
        forward_index = index + 1 if not (index + 1) > (len(array_like) + 1) else index
        backward_index =  index - 1 if not (index - 1) < 0 else index
        if array_like[forward_index] >= array_like[index] and array_like[backward_index] <= array_like[index]:
            highs.append(index)
        elif array_like[forward_index] <= array_like[index] and array_like[backward_index] >= array_like[index]:
            lows.append(index)
    highs = numpy.asarray(highs)
    lows = numpy.asarray(lows)
    maximum_sst = []
    for start in highs:
        try:
            low = next(x for x in lows if x > start)
            max_idx = numpy.argmax(array_like[start:low]) + start
            maximum_sst.append(max_idx)
        except StopIteration as e:
            continue
        except Exception as e:
            print repr(e)
            print e
            continue
    maximum_sst = numpy.asarray(maximum_sst)
    sst_bloom_timings = []
    true_durations = []
    start_end_list = numpy.squeeze(start_end_list)
    for bloom in start_end_list:
        if all(b is None for b in bloom):
            continue
        try:
            bloom_sst_max = (numpy.abs(maximum_sst - bloom[3])).argmin() + bloom[3]
            bloom_sst_start = (numpy.abs(highs - bloom[0])).argmin() + bloom[0]
            bloom_sst_end = (numpy.abs(lows - bloom[1])).argmin() + bloom[1]
            bloom_sst_duration = bloom_sst_end - bloom_sst_start
            #test the sst max is within the start (bloom[0])/end (bloom[1]) dates
            #I think this is wrong, shouldn't it be re testing the chlorophyll max between the two sst dates?
            if bloom_sst_max < bloom[1] and not bloom_sst_max < bloom[0]:
                #the chlorophyll max should be within the boundary of the sst
                true_max = bloom_sst_max
                true_max_val = chl_slice[bloom_sst_max]
            else:
                true_max = bloom[2]
                true_max_val =  chl_slice[bloom[2]]
            #check its not a bigger difference than 2 weeks, this hasn't actually been triggered in testing
            # and not bloom_sst_end > bloom[1] + 2
            if bloom_sst_end > bloom[1]:
                true_end = bloom_sst_end
            else:
                true_end = bloom[1]
            # and not bloom_sst_start < bloom[0] - 5
            if bloom_sst_start < bloom[0]:
                true_start = bloom_sst_start
            else:
                true_start = bloom[0]

            true_duration = true_end - true_start
            sst_bloom_timings.append([true_start-reverse_search, true_end-reverse_search,true_duration,true_max-reverse_search,true_max_val])
            true_durations.append(true_duration)
        except ValueError as e:
            if (maximum_sst.size > 1 and highs.size > 1 and lows.size > 1):
                print repr(e)
                print e
                print maximum_sst, highs, lows
            #possibly we should just put the chlorophyll dates in here?
    #this is some certified genuine python magic, *(sometimes called a splat) unpacks the zipped, sorted lists into a new zip that outputs to lists through a comprehension
    try:
        sst_sorted_durations, sst_sorted_durations_idx = (list(t) for t in zip(*sorted(zip(true_durations, range(0,len(true_durations))))))
    except ValueError:
        maximums = [[-1000,-1000,-1000,-1000,-1000], [-1000,-1000,-1000,-1000,-1000]]
    else:
     #sst_sorted_durations, sst_sorted_durations_idx = zip(*sorted(zip(sst_durations, range(0,len(sst_durations)))
        if len(sst_sorted_durations) >= 2:
            maximums = [sst_bloom_timings[sst_sorted_durations_idx[0]],sst_bloom_timings[sst_sorted_durations_idx[1]]]
        elif len(sst_sorted_durations) == 1:
            maximums = [sst_bloom_timings[sst_sorted_durations_idx[0]], [-1000,-1000,-1000,-1000,-1000]]
        else:
            maximums = [[-1000,-1000,-1000,-1000,-1000], [-1000,-1000,-1000,-1000,-1000]]
    return maximums

def prepare_sst_variables(sst_array, numpy_storage):
    #smoothed sst
    print "sst sbx"
    sst_boxcar = numpy.apply_along_axis(numpy.convolve, 0, sst_array, numpy.ones((8,))/8, mode='valid')
    sst_boxcar_map = numpy.memmap(os.path.join(numpy_storage, "sst_sbx"), mode="w+", shape=sst_boxcar.shape, dtype=sst_boxcar.dtype)
    sst_boxcar_map[:] = sst_boxcar[:]
    sst_boxcar = None
    sst_boxcar = sst_boxcar_map
    sst_boxcar_map = None
    sst_array = None

    #get sst derivative
    print "sst der"
    sst_der = numpy.apply_along_axis(centered_diff_derivative, 0, sst_boxcar)
    sst_der_map = numpy.memmap(os.path.join(numpy_storage, "sst_der"), mode="w+", shape=sst_der.shape, dtype=sst_der.dtype)
    sst_der_map[:] = sst_der[:]
    return sst_der.shape, sst_der.dtype

def prepare_chl_variables(chl_array, numpy_storage):
    #median * 1.05
    print "med5"
    med5 = numpy.ma.median(chl_array,axis = 0)*1.05
    #get anomalies
    print "anomaly"
    anomaly = chl_array - (med5/1.05*1.2)
    anomaly_map = numpy.memmap(os.path.join(numpy_storage, "chl_anomaly"), mode="w+", shape=anomaly.shape, dtype=anomaly.dtype)
    anomaly_map[:] = anomaly[:]
    anomaly = anomaly_map
    chl_array = None
    anomaly_map = None
    #need to ditch any empty entries here as they interefere with the cumsum

    #get cumsum of anomalies
    print "chl cumsum"
    chl_cumsum = numpy.ma.cumsum(anomaly,axis=1)
    chl_cumsum_map = numpy.memmap(os.path.join(numpy_storage, "chl_cumsum"), mode="w+", shape=chl_cumsum.shape, dtype=chl_cumsum.dtype)
    chl_cumsum_map[:] = chl_cumsum[:]
    chl_cumsum = chl_cumsum_map
    anomaly = None
    chl_cumsum_map = None


    #get centered derivative
    print "chl der"
    chl_der = numpy.apply_along_axis(centered_diff_derivative, 0, chl_cumsum)
    chl_der_map = numpy.memmap(os.path.join(numpy_storage, "chl_der"), mode="w+", shape=chl_der.shape, dtype=chl_der.dtype)
    chl_der_map[:] = chl_der[:]
    #chl_der = chl_der_map
    chl_cumsum = None
    chl_der_map = None

    #boxcar filter with width of 3 (sbx) should be something like this:
    print "chl sbx"
    chl_boxcar = numpy.apply_along_axis(numpy.convolve, 0, chl_der, numpy.ones((8,))/8, mode='valid')
    chl_boxcar_map = numpy.memmap(os.path.join(numpy_storage, "chl_sbx"), mode="w+", shape=chl_boxcar.shape, dtype=chl_boxcar.dtype)
    chl_boxcar_map[:] = chl_boxcar[:]
    chl_boxcar = None
    chl_boxcar = chl_boxcar_map
    chl_boxcar_map = None
    return chl_boxcar.shape, chl_boxcar.dtype

def create_phenology_netcdf(chl_lons, chl_lats, output_shape=None,name="phenology.nc"):
    global output_location
    output_location = name
    ds = nc.Dataset(name,'w',format='NETCDF4_CLASSIC')
    ds.createDimension('LONGITUDE', output_shape[3])
    ds.createDimension('LATITUDE', output_shape[2])
    ds.createDimension('DEPTH', output_shape[1])
    ds.createDimension('TIME', None)
    ds.createVariable('LATITUDE', 'float32', dimensions=['LATITUDE'])
    ds.variables['LATITUDE'].setncattr("units", "degrees north")
    ds.variables['LATITUDE'][:] = chl_lats
    ds.createVariable('LONGITUDE', 'float32', dimensions=['LONGITUDE'])
    ds.variables['LONGITUDE'].setncattr("units", "degrees east")
    ds.variables['LONGITUDE'][:] = chl_lons
    ds.createVariable('DEPTH', 'float32', dimensions=['DEPTH'])
    ds.variables['DEPTH'].setncattr("units", "metres")
    ds.createVariable('TIME', 'float32', dimensions=['TIME'])
    ds.variables['TIME'].setncattr("units", "years")
    ds.createVariable('date_start1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_start1'].setncattr("units", "weeks")
    ds.createVariable('date_max1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_max1'].setncattr("units", "weeks")
    ds.createVariable('date_end1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_end1'].setncattr("units", "weeks")
    ds.createVariable('duration1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['duration1'].setncattr("units", "weeks")
    ds.createVariable('date_start2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_start2'].setncattr("units", "weeks")
    ds.createVariable('date_max2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_max2'].setncattr("units", "weeks")
    ds.createVariable('date_end2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['date_end2'].setncattr("units", "weeks")
    ds.createVariable('duration2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['duration2'].setncattr("units", "weeks")
    ds.createVariable('max_val1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.createVariable('max_val2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=-1000)
    ds.variables['max_val2'].setncattr("units", "mgChl/m3")
    ds.variables['max_val1'].setncattr("units", "weeks")
    ds.close()
    print "created netcdf {}".format(name)

def write_to_output_netcdf(data):
    ds = nc.Dataset(output_location,'r+',format='NETCDF4_CLASSIC')
    print "pre-writing data shape: {}".format(data.shape)
    year = ds.variables['date_start1'][:].shape[0]
    """
    if not ds.variables['date_start1'][:].shape[0]:
        netcdf_data_shape = (+1, 1, data.shape[2], data.shape[3])
        print "start1 shapes: {} {}".format(ds.variables['date_start1'][:].shape, data[0][0].shape)
        ds.variables['date_start1'][:] = numpy.append(ds.variables['date_start1'][:], data[0][0]).reshape(netcdf_data_shape)
    else:
    """
    ds.variables['date_start1'][year] = data[:,:,0,0]
    print "max1 shapes: {} {}".format(ds.variables['date_max1'][:].shape, data[1].shape)
    ds.variables['date_max1'][year] = data[:,:,0,3]
    ds.variables['date_end1'][year] = data[:,:,0,1]
    ds.variables['duration1'][year] = data[:,:,0,2]
    ds.variables['date_start2'][year] = data[:,:,1,0]
    ds.variables['date_max2'][year] = data[:,:,1,3]
    ds.variables['date_end2'][year] = data[:,:,1,1]
    ds.variables['duration2'][year] = data[:,:,1,2]
    ds.variables['max_val1'][year] = data[:,:,0,4]
    ds.variables['max_val2'][year] = data[:,:,1,4]
    print ds.variables['TIME'][:]
    ds.variables['TIME'][year] = year
    """
    ds.variables['date_start1'][:] = numpy.append(ds.variables['date_start1'][:], data[0][0])
    print "max1 shapes: {} {}".format(ds.variables['date_max1'][:].shape, data[0][1].shape)
    ds.variables['date_max1'][:] = numpy.append(ds.variables['date_max1'][:], data[0][1])
    ds.variables['date_end1'][:] = numpy.append( ds.variables['date_end1'][:], data[0][2])
    ds.variables['duration1'][:] = numpy.append(ds.variables['duration1'][:], data[0][3])
    ds.variables['date_start2'][:] = numpy.append(ds.variables['date_start2'][:], data[1][0])
    ds.variables['date_max2'][:] = numpy.append(ds.variables['date_max2'][:], data[1][1])
    ds.variables['date_end2'][:] = numpy.append(ds.variables['date_end2'][:], data[1][2])
    ds.variables['duration2'][:] = numpy.append(ds.variables['duration2'][:], data[1][3])
    """
    ds.close()



def get_multi_year_two_blooms_output(numpy_storage, chl_shape, chl_dtype, sst_shape, sst_dtype, date_seperation_per_yr=47, start_date=0, reverse_search=20):
    #this all works on the assumption the axis 0 is time
    print "reading variables"
    chl_boxcar = numpy.memmap(os.path.join(numpy_storage, "chl_sbx"), mode="r", dtype=chl_dtype, shape=chl_shape)
    sst_der = numpy.memmap(os.path.join(numpy_storage, "sst_der"), mode="r", dtype=sst_dtype, shape=sst_shape)
    print "shapes after reading sst: {} chl: {}".format(chl_boxcar.shape, sst_der.shape)
    print "reshaping to sst: {} chl: {}".format(sst_shape, chl_shape)
    sst_der.shape = sst_shape
    chl_boxcar.shape = chl_shape

    print "doing yearly evaluation, this may take up a lot of memory, if so double check all memmaps have been flushed"
    for year in range(start_date, chl_boxcar.shape[0], date_seperation_per_yr):
        print "doing year {}".format(year / date_seperation_per_yr)
        prev_year = year - reverse_search if year > (reverse_search - 1) else start_date
        next_year = year + date_seperation_per_yr
        #get start and ends, this works
        year_chl_boxcar = chl_boxcar[prev_year:next_year,:,:,:]
        year_chl_boxcar = numpy.ma.masked_where((year_chl_boxcar == -9.999999999999998e+33), year_chl_boxcar)
        year_sst_der = sst_der[prev_year:next_year,:,:,:]
        year_sst_der = numpy.ma.masked_where((year_sst_der == -9.999999999999998e+33), year_sst_der)
        print "doing chlorophyll initiations"
        start_end_duration_array = numpy.apply_along_axis(get_start_index_and_duration, 0, year_chl_boxcar)
        print "done chlorophyll initiations"
        year_true_start_end_array = numpy.ndarray((start_end_duration_array.shape[2],start_end_duration_array.shape[3], 2,5))
        year_true_start_end_array.fill(-1000)
        print "doing sst initiations and correction"
        for ix in numpy.ndindex(start_end_duration_array.shape[2]):
            for iy in numpy.ndindex(start_end_duration_array.shape[3]):
                year_true_start_end_array[ix,iy] = match_start_end_to_solar_cycle(year_sst_der[:,:,ix,iy], start_end_duration_array[:,:,ix,iy], year_chl_boxcar[:,:,ix,iy], reverse_search)
        print "done sst initiations and correction"
        print "writing to netcdf"
        write_to_output_netcdf(year_true_start_end_array)
        #year_storage = numpy.memmap(os.path.join(numpy_storage, "year_{}".format(year)), mode="w+", shape=year_true_start_end_array.shape, dtype=year_true_start_end_array.dtype)

    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("chl_location", help="A chlorophyll file or glob for files (e.g. chl*.nc) that can be stacked into an array of shape (time,depth,lat,lon). These must have a variable called chl, lat and lon (or contain those strings) which will be used to define the array shape")
    parser.add_argument("--sst_location", help="An sst file, or glob of files (e.g. sst*.nc) that matches the chlorophyll observations. if it does not match some interpolation can be attempted with --sst_date_interp", required=False)
    parser.add_argument("--output", help="output location, default is phenology.nc in current folder", default="phenology.nc", required=False)
    parser.add_argument("--date_seperation_per_year", help="how many temporal observations we have in a year, if not specified will be guessed", default=47, required=False)
    parser.add_argument("--first_date_index", help="specify if the first date you want to include is not the first date present in the date stack", default=0, required=False)
    parser.add_argument("--intermediate_file_store", help="change where intermediate numpy files are placed, if not specified then /tmp is assumed - you should specify somewhere else if your tmp cannot handle the array sizes needed (and currently this program will fill it until it cannot).", required=False)
    args = parser.parse_args()
    if not args.intermediate_file_store:
        numpy_storage = tempfile.mkdtemp(prefix="phenology_")
    else:
        numpy_storage = args.intermediate_file_store


    if args.sst_location:
        print "sst file provided, reading array"
        sst_files = glob.glob(args.sst_location)
        if len(sst_files) == 1:
            print "only one file found, assuming full stack of observations"
            sst_ds = nc.Dataset(sst_files[0])
            sst_variable = [x for x in sst_ds.variables if "sst" in x.lower()][0]
            sst_array = sst_ds.variables[sst_variable][:]
        else:
            raise NotImplementedError
        sst_shape, sst_dtype = prepare_sst_variables(sst_array, numpy_storage)
        print "sst_shape: {}".format(sst_shape)
        print "sst_dtype: {}".format(sst_dtype)
        sst_array = None

    chl_files = glob.glob(args.chl_location)
    if len(chl_files) == 1:
        print "only one chl file found, assuming full stack of observations"
        chl_ds = nc.Dataset(chl_files[0])
        chl_variable = [x for x in chl_ds.variables if "chl" in x.lower()][0]
        chl_array = chl_ds.variables[chl_variable][:]
    else:
        raise NotImplementedError
    chl_shape, chl_dtype = prepare_chl_variables(chl_array, numpy_storage)
    print "chl_shape: {}".format(chl_shape)
    print "chl_dtype: {}".format(chl_dtype)

    chl_ds = nc.Dataset(args.chl_location)
    chl_variable = [x for x in chl_ds.variables if "chl" in x.lower()][0]
    chl_lon_var = [x for x in chl_ds.variables if "lon" in x.lower()][0]
    chl_lat_var = [x for x in chl_ds.variables if "lat" in x.lower()][0]
    chl_lons = chl_ds.variables[chl_lon_var][:]
    chl_lats = chl_ds.variables[chl_lat_var][:]
    chl_shape = (821, 1, 396, 1297)
    sst_shape = (821, 1, 396, 1297)
    sst_dtype = 'float64'
    chl_dtype = 'float64'
    print "creating output netcdf {}".format(args.output)
    create_phenology_netcdf(chl_lons, chl_lats, chl_shape, args.output)

    get_multi_year_two_blooms_output(numpy_storage, chl_shape, chl_dtype, sst_shape, sst_dtype, date_seperation_per_yr=args.date_seperation_per_yr, start_date=args.first_date_index, reverse_search=20)

    
    

