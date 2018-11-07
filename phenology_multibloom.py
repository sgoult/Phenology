from __future__ import print_function, division
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

FILL_VAL = -9.999999999999998e+33

MEDIAN_THRESHOLD_DEFAULT = 5

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

def get_start_index_and_duration(array_like,chl_values,date_offset,depth=5, pad_values=False, verbose=False):
    """
    takes a list of values and searches for the max, then the associated sign changes to indicate bloom start/end
    set depth to change limit of number of identifications
    
!comment[A] my approach was to save the  max Chl value found between each of the start and end time, and also the duration, i.e., estimated as number of steps between start and end times
! with the first derivative light availability or SST is increasing when PAR or SST first derivative is positive, and vice versa
 
    """
    #array_like = numpy.squeeze(array_like)
    #if it's all gone horribly wrong then this will quit out of it straight away
    if len(array_like):
        zero_crossings = list(numpy.where(numpy.diff(numpy.sign(array_like)))[0])
    else:
        zero_crossings = []
    true_poss = zero_crossings


    if verbose:
        print("zero crossings")
        print(zero_crossings)
        print("chl sbx")
        print(array_like)
    #find out which way we're going
    starts = []
    ends = []
    #works much the same qas the SST one below
    for index in true_poss:
        forward_index = index + 1 if not (index + 1) >= (len(array_like)) else index
        backward_index =  index - 1 if not (index - 1) < 0 else index
        if array_like[forward_index] >= array_like[index] and array_like[backward_index] <= array_like[index]:
            starts.append(index)
        elif array_like[forward_index] <= array_like[index] and array_like[backward_index] >= array_like[index]:
            ends.append(index)
    #we know the last entry will be an end
    ends.append(len(array_like))
    if verbose:
        print("starts and ends")
        print(starts)
        print(ends)
    #find an end for every start
    dates = []
    for start in starts:
        try:
           end = next(x for x in ends if x > start)
           if verbose:
               print("chl values")
               print(chl_values[start:end])
           max_idx = numpy.nanargmax(chl_values[start:end]) + start
           dates.append([start + date_offset,end + date_offset,end-start,max_idx,chl_values[max_idx]])
        except StopIteration:
            continue
        except Exception as e:
           print(e)
           print(repr(e))
           continue
    if pad_values:
        for pad in range(len(dates), depth):
            dates.append([None,None,None,None,None])
    if verbose:
        print("end dates")
        print(dates)
    return dates

def phen_records_to_one_val_on_max(records, date_correction=False, index=4):
    """
    Reduces a list to one value, and corrects dates based on variable. Selects max chlorophyll as default, change index to whatever to sort by something else
    """
    if len(records):
        if len(records) == 1:
            output_record = records[0]
        else:
            maxes = [x[index] for x in records]
            maximum = maxes.index(max(maxes))
            output_record = records[maximum]
        if date_correction:
            output_record[0] = output_record[0] - date_correction
            output_record[1] = output_record[1] - date_correction
            #TODO double check this is actually correcting the max idx value
            output_record[3] = output_record[3] - date_correction
        return output_record
    else:
        return [None,None,None,None,None]

def match_start_end_to_solar_cycle(array_like, chl_sbx_slice, chl_slice, date_seperation_per_year, reverse_search, start_date=0, verbose=False):
    """
    Attributes the start and end times in relation to the SST or solar cycle, takes an sst array (array_like), smoothed chlorophyll derivative slice (chl_sbx_slice) and the original chlorophyll data.

    Slices up the data based on high/low periods of SST (or otherwise), then feeds each period into get_start_index_and_duration, once finished it will output an array of shape (x, y, time, 2, 5)
    verbose will spame the terminal with information about whats going on, best to establish a few pixels you want to inspect rather than having this on all the time.
    """
    
    #possibly resort and create new durations based on remaining dates
    #look for sign changes in sst or PAR data, indicating high/low SST or light periods
    array_like = numpy.squeeze(array_like)
    chl_slice = numpy.squeeze(chl_slice)
    chl_sbx_slice = numpy.squeeze(chl_sbx_slice)
    zero_crossings = numpy.where(numpy.diff(numpy.sign(array_like)))[0]

    #find out which way we're going
    highs = []
    lows = []
    for index in zero_crossings:
        #checks if the values are increasing or decreasing
        forward_index = index + 1 if not (index + 1) > (len(array_like) + 1) else index
        backward_index =  index - 1 if not (index - 1) < 0 else index
        #add to high period or low periods
        if array_like[forward_index] >= array_like[index] and array_like[backward_index] <= array_like[index]:
            highs.append(index)
        elif array_like[forward_index] <= array_like[index] and array_like[backward_index] >= array_like[index]:
            lows.append(index)

    #print(everything thus far
    if verbose:
        print("***********************")
        print("highs and lows")
        print(highs)
        print(lows)
    maximum_sst = []
    max_yr_idx = len(array_like)
    activity_period = None
    try:
        """
        if we have identified some high and low periods check which one we start with, then add the opposite one to pad back to the beginning of the year (otherwise we can miss a lot of data)

        This generally won't get used if you are working with more than one year of data
        """
        if len(highs) and len(lows):
            if not highs[0] == 0 and not lows[0] == 0:
                if highs[0] > lows[0]:
                    highs = [0] + highs
                else:
                    lows = [0] + lows
            if highs[-1] > lows[-1]:
                highs.append(len(array_like))
            else:
                lows.append(len(array_like))
        elif len(highs) and not len(lows):
            lows = [0, len(array_like)]
        elif len(lows) and not len(highs):
            highs = [0, len(array_like)]
        else:
            return [[None,None,None,None,None], [None,None,None,None,None]]
    except Exception as e:
        #triggered once, but was helpful to know what the contents were
        print(highs)
        print(lows)
        print(e)

    #chuck out the results
    if verbose:
        print("updated highs and lows")
        print(highs)
        print(lows)

    high_records = []
    low_records = []
    for index in zero_crossings:
        #we've classified them, this just lets us select out the high and low period
        if index in highs:
            activity_period = 1
            try:
                end_date = next(x for x in lows if x > index)
            except StopIteration as e:
                continue
        else:
            activity_period = 0
            try:
                end_date = next(x for x in highs if x > index)
            except StopIteration as e:
                continue
        first = False
        
        chl_sbx_period_slice = chl_sbx_slice[index:end_date].flatten()
        chl_period_slice = chl_slice[index:end_date].flatten()
        #get the phenology for this period, depth pads extra data if needed for numpy (we don't use this for SST model)
        period_chl_phenology = get_start_index_and_duration(chl_sbx_period_slice,chl_period_slice,index,depth=5,verbose=verbose)
        #if we found anything
        if len(period_chl_phenology):
            #loop through them and add them to the high/low mega lists
            for record in period_chl_phenology:
                if activity_period:
                    high_records.append(record)
                else:
                    low_records.append(record)
    
    #to remind ourselves what the phenology records look like
    #[start,end,end-start,max_idx,chl_values[max_idx]]
    blooms = []
    for year in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year):
        #find blooms that start after the year - our reverse search, end before the end of the year, and end during the current year
        possible_high_blooms = [x for x in high_records if x[0] > (year - reverse_search) and x[1] < (year + date_seperation_per_year) and x[1] > year]
        possible_low_blooms = [x for x in low_records if x[0] > (year - reverse_search) and x[1] < (year + date_seperation_per_year) and x[1] > year]
        if verbose:
            print(year)
            print("found ", len(possible_high_blooms), " highs")
            print("found ", len(possible_low_blooms), " lows")

        #reduce them to one record
        low = phen_records_to_one_val_on_max(possible_low_blooms, year)
        high = phen_records_to_one_val_on_max(possible_high_blooms, year)
        if not low or not high:
            print("***************")
            print(high_records)
            print(low_records)
            print(year - reverse_search)
            print(possible_high_blooms)
            print(possible_low_blooms)
            print(low, high)
        #spit out the low period and high period for this year
        blooms.append([low,high])

    return blooms

def prepare_sst_variables(sst_array, numpy_storage):
    """
    Creates smoothed sst, currently has a large portion commented out as source file is already the centered derivative diff data.
    """

    #smoothed sst
    print("sst sbx")
    """
    sst_boxcar = numpy.apply_along_axis(numpy.convolve, 0, sst_array, numpy.ones((8,))/8, mode='valid')
    sst_boxcar_map = numpy.memmap(os.path.join(numpy_storage, "sst_sbx"), mode="w+", shape=sst_boxcar.shape, dtype=sst_boxcar.dtype)
    sst_boxcar_map[:] = sst_boxcar[:]
    sst_boxcar = None
    sst_boxcar = sst_boxcar_map
    sst_boxcar_map = None
    sst_array = None
    #get sst derivative
    print("sst der")
    sst_der = numpy.apply_along_axis(centered_diff_derivative, 0, sst_boxcar)
    """
    sst_der = sst_array
    sst_der_map = numpy.memmap(os.path.join(numpy_storage, "sst_der"), mode="w+", shape=sst_der.shape, dtype=sst_der.dtype)
    sst_der_map[:] = sst_der[:]
    return sst_der.shape, sst_der.dtype

def prepare_chl_variables(chl_array, numpy_storage, median_threshold=None):
    """
    Creates the smoothed anomaly chlorophyll data, saves a file to the temporary directory that is read as a mem map later to conserve resources.
    """
    #median * 1.05
    
    #! here i would prefer we output the media value (without adding 5% to it)
    
    print("med5")
    med5 = numpy.ma.median(chl_array,axis = 0)
    print("median value: {}".format(med5))
    print("median threshold: {}".format(median_threshold))


    #! after estimating the median, i apply a filling to the chl time-series, e.g., window of 8 time steps, but it would be good to be able to choose other width of the window, e.g., 5 time-steps or 3 time-steps...

    #! here it would be good to give the option to select the value of the median threshold, e.g., median plus 5%, 10%, 15%...
    
    #get anomalies
    print("anomaly")
    anomaly = chl_array - (med5*1.2)
    anomaly_map = numpy.memmap(os.path.join(numpy_storage, "chl_anomaly"), mode="w+", shape=anomaly.shape, dtype=anomaly.dtype)
    anomaly_map[:] = anomaly[:]
    anomaly = anomaly_map
    chl_array = None
    anomaly_map = None
    #need to ditch any empty entries here as they interefere with the cumsum

    #get cumsum of anomalies
    print("chl cumsum")
    chl_cumsum = numpy.ma.cumsum(anomaly,axis=1)
    chl_cumsum_map = numpy.memmap(os.path.join(numpy_storage, "chl_cumsum"), mode="w+", shape=chl_cumsum.shape, dtype=chl_cumsum.dtype)
    chl_cumsum_map[:] = chl_cumsum[:]
    chl_cumsum = chl_cumsum_map
    anomaly = None
    chl_cumsum_map = None


    #get centered derivative        
    print("chl der")
    chl_der = numpy.apply_along_axis(centered_diff_derivative, 0, chl_cumsum)
    chl_der_map = numpy.memmap(os.path.join(numpy_storage, "chl_der"), mode="w+", shape=chl_der.shape, dtype=chl_der.dtype)
    chl_der_map[:] = chl_der[:]
    #chl_der = chl_der_map
    chl_cumsum = None
    chl_der_map = None

    #boxcar filter with width of 3 (sbx) should be something like this:
    print("chl sbx")
    chl_boxcar = numpy.apply_along_axis(numpy.convolve, 0, chl_der, numpy.ones((8,))/8, mode='valid')
    chl_boxcar_map = numpy.memmap(os.path.join(numpy_storage, "chl_sbx"), mode="w+", shape=chl_boxcar.shape, dtype=chl_boxcar.dtype)
    chl_boxcar_map[:] = chl_boxcar[:]
    chl_boxcar = None
    chl_boxcar = chl_boxcar_map
    chl_boxcar_map = None
    return chl_boxcar.shape, chl_boxcar.dtype

def create_phenology_netcdf(chl_lons, chl_lats, output_shape=None,name="phenology.nc"):
    """
    Creates the skeleton of the netcdf file to be used by write_to_output_netcdf, all of this is metadata.
    """

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
    ds.createVariable('date_start1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['date_start1'].setncattr("units", "weeks")
    ds.createVariable('date_max1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['date_max1'].setncattr("units", "weeks")
    ds.createVariable('date_end1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['date_end1'].setncattr("units", "weeks")
    ds.createVariable('duration1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['duration1'].setncattr("units", "weeks")
    ds.createVariable('date_start2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['date_start2'].setncattr("units", "weeks")
    ds.createVariable('date_max2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['date_max2'].setncattr("units", "weeks")
    ds.createVariable('date_end2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['date_end2'].setncattr("units", "weeks")
    ds.createVariable('duration2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['duration2'].setncattr("units", "weeks")
    ds.createVariable('max_val1', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.createVariable('max_val2', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ds.variables['max_val2'].setncattr("units", "mgChl/m3")
    ds.variables['max_val1'].setncattr("units", "weeks")
    ds.close()
    print("created netcdf {}".format(name))

def write_to_output_netcdf(data):
    """
    Loops through each year in the numpy array and writes the data to the netcdf file, this should work faster if we get rid of the loop but I can't seem to grock the logic to fix it right now.
    """
    ds = nc.Dataset(output_location,'r+',format='NETCDF4_CLASSIC')
    data = data.astype(np.float32)
    print(output_location)
    print("pre-writing data shape: {}".format(data.shape))
    year = ds.variables['date_start1'][:].shape[0]
    print(data[:,:,:,0,0].shape)
    ds.variables['TIME'][:] = range(0, data.shape[2])
    for year in range(0, data.shape[2] -1):
        ds.variables['date_start1'][year] = data[:,:,year,0,0]
        ds.variables['date_max1'][year] = data[:,:,year,0,3]
        ds.variables['date_end1'][year] = data[:,:,year,0,1]
        ds.variables['duration1'][year] = data[:,:,year,0,2]
        ds.variables['date_start2'][year] = data[:,:,year,1,0]
        ds.variables['date_max2'][year] = data[:,:,year,1,3]
        ds.variables['date_end2'][year] = data[:,:,year,1,1]
        ds.variables['duration2'][year] = data[:,:,year,1,2]
        ds.variables['max_val1'][year] = data[:,:,year,0,4]
        ds.variables['max_val2'][year] = data[:,:,year,1,4]
    print(ds.variables['TIME'][:])
    ds.close()



def get_multi_year_two_blooms_output(numpy_storage, chl_shape, chl_dtype, chl_data, sst_shape, sst_dtype, date_seperation_per_year=47, start_date=0, reverse_search=20):
    #this all works on the assumption the axis 0 is time
    print("reading variables")
    chl_boxcar = numpy.memmap(os.path.join(numpy_storage, "chl_sbx"), mode="r", dtype=chl_dtype, shape=chl_shape)
    sst_der = numpy.memmap(os.path.join(numpy_storage, "sst_der"), mode="r", dtype=sst_dtype, shape=sst_shape)
    print("shapes after reading sst: {} chl: {}".format(chl_boxcar.shape, sst_der.shape))
    print("reshaping to sst: {} chl: {}".format(sst_shape, chl_shape))
    sst_der.shape = sst_shape
    chl_boxcar.shape = chl_shape

    print("doing yearly evaluation, this may take up a lot of memory, if so double check all memmaps have been flushed")
    #get start and ends, this works
    chl_boxcar = numpy.ma.masked_where((chl_boxcar == FILL_VAL), chl_boxcar)
    sst_der = numpy.ma.masked_where((sst_der == FILL_VAL), sst_der)
    #print("doing chlorophyll initiations")
    #start_end_duration_array = numpy.apply_along_axis(get_start_index_and_duration, 0, year_chl_boxcar)
    year_true_start_end_array = numpy.ndarray((chl_data.shape[2],chl_data.shape[3], int(abs(chl_data.shape[0] / date_seperation_per_year)), 2,5))
    year_true_start_end_array.fill(FILL_VAL)
    completion_points = range(0, chl_data.shape[2], chl_data.shape[2] // 10)
    print("doing sst initiations and correction")
    for ix in numpy.ndindex(chl_data.shape[2]):
        for iy in numpy.ndindex(chl_data.shape[3]):
            try:
                verbose=False
                if ix[0] > 85 and ix[0] < 95 and iy[0] > 180 and iy[0] < 190:
                    verbose = True
                year_true_start_end_array[ix,iy] = match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=False, start_date=start_date)
                if verbose:
                    print("end duration array")
                    print(year_true_start_end_array[ix,iy])
            except Exception as e:
                print(e)
                print(repr(e))
                print(match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=False, start_date=start_date))
                print()
        if ix[0] in completion_points:
            print(completion_points.index(ix[0]) * 10, "% complete")

    print("done sst initiations and correction")
    print("writing to netcdf")
    #needs to be extended to be able to output 3 files: sorted by calendar year, one with primary and secondary chloropjhyll maximum
    write_to_output_netcdf(year_true_start_end_array)

    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("chl_location", help="A chlorophyll file or glob for files (e.g. chl*.nc) that can be stacked into an array of shape (time,depth,lat,lon). These must have a variable called chl, lat and lon (or contain those strings) which will be used to define the array shape")
    parser.add_argument("--sst_location", help="An sst file, or glob of files (e.g. sst*.nc) that matches the chlorophyll observations. if it does not match some interpolation can be attempted with --sst_date_interp", required=False)
    parser.add_argument("--output", help="output location, default is phenology.nc in current folder", default="phenology.nc", required=False)
    parser.add_argument("--date_seperation_per_year", help="how many temporal observations we have in a year, if not specified will be guessed", default=47, required=False)
    parser.add_argument("--first_date_index", help="specify if the first date you want to include is not the first date present in the date stack", default=0, required=False)
    parser.add_argument("--intermediate_file_store", help="change where intermediate numpy files are placed, if not specified then /tmp is assumed - you should specify somewhere else if your tmp cannot handle the array sizes needed (and currently this program will fill it until it cannot).", required=False)
    parser.add_argument("--median_threshold", default=MEDIAN_THRESHOLD_DEFAULT, help="change median threshold", required=False)
    args = parser.parse_args()
    if not args.intermediate_file_store:
        numpy_storage = tempfile.mkdtemp(prefix="phenology_")
    else:
        numpy_storage = args.intermediate_file_store
    #remember to change median threshold to percentage!

    if args.sst_location:
        print("sst file provided, reading array")
        sst_files = glob.glob(args.sst_location)
        if len(sst_files) == 1:
            print("only one file found, assuming full stack of observations")
            sst_ds = nc.Dataset(sst_files[0])
            sst_variable = [x for x in sst_ds.variables if "sst" in x.lower()][0]
            sst_array = sst_ds.variables[sst_variable][:]
        else:
            raise NotImplementedError
        sst_shape, sst_dtype = prepare_sst_variables(sst_array, numpy_storage)
        print("sst_shape: {}".format(sst_shape))
        print("sst_dtype: {}".format(sst_dtype))
        sst_array = None

    chl_files = glob.glob(args.chl_location)
    if len(chl_files) == 1:
        print("only one chl file found, assuming full stack of observations")
        chl_ds = nc.Dataset(chl_files[0])
        chl_variable = [x for x in chl_ds.variables if "chl" in x.lower()][0]
        chl_array = chl_ds.variables[chl_variable][:]
    else:
        raise NotImplementedError
    chl_shape, chl_dtype = prepare_chl_variables(chl_array, numpy_storage)
    print("chl_shape: {}".format(chl_shape))
    print("chl_dtype: {}".format(chl_dtype))

    chl_ds = nc.Dataset(args.chl_location)
    chl_variable = [x for x in chl_ds.variables if "chl" in x.lower()][0]
    chl_lon_var = [x for x in chl_ds.variables if "lon" in x.lower()][0]
    chl_lat_var = [x for x in chl_ds.variables if "lat" in x.lower()][0]
    chl_lons = chl_ds.variables[chl_lon_var][:]
    chl_lats = chl_ds.variables[chl_lat_var][:]
    sst_dtype = 'float64'
    chl_dtype = 'float64'
    chl_data = chl_ds[chl_variable]
    print("creating output netcdf {}".format(args.output))

    print(chl_shape)
    create_phenology_netcdf(chl_lons, chl_lats, chl_shape, args.output)

    get_multi_year_two_blooms_output(numpy_storage, chl_shape, chl_dtype, chl_data, sst_shape, sst_dtype, date_seperation_per_year=int(args.date_seperation_per_year), start_date=args.first_date_index, reverse_search=20)



    """
    old code from match_start_end
    for start in highs:
        #this looks for the end (ie the start of a low period) to every start of a high period
        try:
            #get the end date
            low = next(x for x in lows if x > start)
            #select the data between the start and end of the high period
            max_idx = numpy.argmax(array_like[start:low]) + start
            #append the sst maximum
            
            #! not sure if we need to save the SST maximum, but need to save the Chl maximum and the duration
            
            maximum_sst.append(max_idx)
        except StopIteration as e:
            #we've run out of values to sequence through
            continue
        except Exception as e:
            #we don't recognise this error - print(it!
            print(repr(e)
            print(e
            continue
    maximum_sst = numpy.asarray(maximum_sst)
    sst_bloom_timings = []
    true_durations = []
    start_end_list = numpy.squeeze(start_end_list)
    for bloom in start_end_list:
        #bloom in this case is a 0 indexed list representing [start, end, duration, max_idx]
        if all(b is None for b in bloom):
            continue
        try:
            bloom_sst_max = (numpy.abs(maximum_sst - bloom[3])).argmin() + bloom[3]
            bloom_sst_start = (numpy.abs(highs - bloom[0])).argmin() + bloom[0]
            bloom_sst_end = (numpy.abs(lows - bloom[1])).argmin() + bloom[1]
            bloom_sst_duration = bloom_sst_end - bloom_sst_start
            #test the sst max is within the start (bloom[0])/end (bloom[1]) dates
            #I think this is wrong, shouldn't it be re testing the chlorophyll max between the two sst dates?

            # !yes, please see !comment[A] at the beginning
            
            if bloom_sst_max < bloom[1] and not bloom_sst_max < bloom[0]:
                #the chlorophyll max should be within the boundary of the sst

                #!the Chl max should be within the boundary of the start and end times
                
                true_max = bloom_sst_max
                true_max_val = chl_slice[bloom_sst_max]
            else:
                true_max = bloom[2]
                true_max_val =  chl_slice[bloom[2]]
            #check its not a bigger difference than 2 weeks, this hasn't actually been triggered in testing
            # and not bloom_sst_end > bloom[1] + 2

            # !please see my comment above - i don't think i had put a condition on the bloom duration length at this point. The condition was to keep the start and end dates for the longest bloom duration within one period of high SST or PAR and the same for a period of low SST and PAR 
            
            if bloom_sst_end > bloom[1]:
                true_end = bloom_sst_end
            else:
                true_end = bloom[1]
            # and not bloom_sst_start < bloom[0] - 5

            # !not sure i understand the -5 here
            
            if bloom_sst_start < bloom[0]:
                true_start = bloom_sst_start
            else:
                true_start = bloom[0]

            true_duration = true_end - true_start
            sst_bloom_timings.append([true_start-reverse_search, true_end-reverse_search,true_duration,true_max-reverse_search,true_max_val])
            true_durations.append(true_duration)
        except ValueError as e:
            if (maximum_sst.size > 1 and highs.size > 1 and lows.size > 1):
                print(repr(e)
                print(e
                print(maximum_sst, highs, lows
            #possibly we should just put the chlorophyll dates in here?

            # ! sorry not clear to me what the code is doing here, we can discuss when we meet
            
    #this is some certified genuine python magic, * (sometimes called a splat) unpacks the zipped, sorted lists into a new zip that outputs to lists through a comprehension
    try:
        #check the selection criteria for the first and second selection
        #think this should actually be finding the longest bloom in the high and then the longest duration in the low periods

        # ! yes correct :)
        #add some logic here to allow user to select either durtation or maximum chlorophyll (or something else maybe)
        sst_sorted_durations, sst_sorted_durations_idx = (list(t) for t in zip(*sorted(zip(true_durations, range(0,len(true_durations))))))
    except ValueError:
        maximums = [[FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL], [FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL]]
    else:
     #sst_sorted_durations, sst_sorted_durations_idx = zip(*sorted(zip(sst_durations, range(0,len(sst_durations)))
        if len(sst_sorted_durations) >= 2:
            maximums = [sst_bloom_timings[sst_sorted_durations_idx[0]],sst_bloom_timings[sst_sorted_durations_idx[1]]]
        elif len(sst_sorted_durations) == 1:
            maximums = [sst_bloom_timings[sst_sorted_durations_idx[0]], [FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL]]
        else:
            maximums = [[FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL], [FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL,FILL_VAL]]
    return maximums

        # ! do you attribute the timing according to the calendar year? I think i do it first for the SST and PAR cycles, and then the timing of chl start and end (based on the SST or PAR cycles) 
    """
