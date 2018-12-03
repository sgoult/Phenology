from __future__ import print_function, division
import os
import numpy
import argparse
import netCDF4 as nc
import shutil
import tempfile
import glob
import sys
import math

#TODO Set dynamically from the input netcdf or user specified from command line
FILL_VAL = -9.999999999999998e+33

#TODO set this from the user inputs
MEDIAN_THRESHOLD_DEFAULT = 20
LAT_IDX = 2
LON_IDX = 3

output_location = None

"""
Solar zenith stuff taken from jad
"""
def computeSunrise(jday, lat):
    theta = (2 * math.pi * jday) / 365.
    delta = 0.006918 - 0.399912 * math.cos(theta) + 0.070257 * math.sin(theta) - 0.006758 * math.cos(2. * theta) + 0.000907 * math.sin(
        2. * theta) - 0.002697 * math.cos(3. * theta) + 0.001480 * math.sin(3. * theta)
    phi = numpy.deg2rad(lat)
    phidel = -math.tan(phi) * math.tan(delta)
    if phidel < -1:
        phidel = -1  # 24 hour sunlight
    elif phidel > 1:
        phidel = 1  # 24 hour darkness
    return 12. - numpy.rad2deg(math.acos(phidel)) / 15., delta, phi

# use time for each element of the time array and the delta/phi are returned from computeSunrise()
def computeZenith(local_time, delta, phi):
    th = (local_time - 12) * (math.pi / 12)
    zen = math.sin(delta) * math.sin(phi) + math.cos(delta) * math.cos(phi) * math.cos(th)
    if zen < -1:
        zen = -1.
    elif zen > 1:
        zen = 1.0
    zen = (math.pi / 2.) - math.asin(zen)
    return zen # radians

def zenithreturn(jday, lat):
    sunrise, delta, phi = computeSunrise(jday, lat)
    #change to 11.5 or 12
    zeniths = [math.degrees(computeZenith(time, delta, phi)) for time in [11, 11.5, 12]]
    return sum(zeniths) / len(zeniths)
"""
End of solar zentih stuff
"""



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
 
    in a run using global data that took 30 minutes, this function made up 513 seconds of the processing time
    """
    #array_like = numpy.squeeze(array_like)
    #if it's all gone horribly wrong then this will quit out of it straight away
    if len(array_like):
        #this is 38.9% of time spent in this function
        zero_crossings = numpy.where(numpy.diff(numpy.sign(array_like)))[0]
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
        #20% of time in function
        if array_like[forward_index] >= array_like[index] and array_like[backward_index] <= array_like[index]:
            starts.append(index)
        #20% of time in function
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
           max_idx = numpy.nanargmax(chl_values[start:end])
           dates.append([start + date_offset,end + date_offset,end-start,max_idx + date_offset + start,chl_values[max_idx + start]])
           max_idx = None
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
            output_record[3] = output_record[3] - date_correction
        return output_record
    else:
        return [None,None,None,None,None]

def match_start_end_to_solar_cycle(array_like, chl_sbx_slice, chl_slice, date_seperation_per_year, reverse_search, start_date=0, verbose=False):
    """
    Attributes the start and end times in relation to the SST or solar cycle, takes an sst array (array_like), smoothed chlorophyll derivative slice (chl_sbx_slice) and the original chlorophyll data.

    Slices up the data based on high/low periods of SST (or otherwise), then feeds each period into get_start_index_and_duration, once finished it will output an array of shape (x, y, time, 2, 5)
    verbose will spame the terminal with information about whats going on, best to establish a few pixels you want to inspect rather than having this on all the time.

    in a run using global data that too 30 minutes, this function made up 703 seconds of the processing time
    I would guess that 500 of those seconds can be attributed to get_start_index_and_duration
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
            return [[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year)]
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
        #10% of time in function
        #each low/high period is treated as a seperate entity, so we can look forward and backward in time without worrying too much about the overlap
        pad = int(round((end_date - index) * 0.20))
        end_date = int(round(end_date + pad))
        start = index - pad if index >= pad else index
        chl_sbx_period_slice = chl_sbx_slice[start:end_date].flatten()
        chl_period_slice = chl_slice[start:end_date].flatten()
        #get the phenology for this period, depth pads extra data if needed for numpy (we don't use this for SST model)
        #73% of time in function
        period_chl_phenology = get_start_index_and_duration(chl_sbx_period_slice,chl_period_slice,start,depth=5,verbose=verbose)
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
    ngds = []

    for year in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year):
        #find blooms that start after the year - our reverse search, end before the end of the year, and end during the current year
        possible_high_blooms = [x for x in high_records if x[0] > (year - reverse_search) and x[1] < (year + date_seperation_per_year) and x[1] > year]
        possible_low_blooms = [x for x in low_records if x[0] > (year - reverse_search) and x[1] < (year + date_seperation_per_year) and x[1] > year]

        high_removals = []
        low_removals = []
        #filters out the blooms that might overlap
        for hindex, high_bloom in enumerate(possible_high_blooms):
            append_high = True
            for lindex, low_bloom in enumerate(possible_low_blooms):
                if any([low_bloom[x] == low_bloom[x] for x in [0, 1]]):
                    #whichever max is higher we should select
                    if low_bloom[4] > high_bloom[4]:
                        high_removals.append(hindex)
                    else:
                        low_removals.append(lindex)
        
        for idx in sorted(set(high_removals), reverse = True):
            del possible_high_blooms[idx]
        
        for idx in sorted(set(low_removals), reverse = True):
            del possible_low_blooms[idx]

                    

        if verbose:
            print(year)
            print("found ", len(possible_high_blooms), " highs")
            print("found ", len(possible_low_blooms), " lows")

        #reduce them to one record
        #additionally look at using max vs duration for the bloom selection
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
        if low[4] > high[4]:
            blooms.append([low,high])
        else:
            blooms.append([high, low])
        #establish if the date is within the year - does this need to be tested for?
        #alternative is (date_seperation_per_year // 0.630136986) to get in first 230 days but this seems wrong
        #we have all of the blooms in a year, could establish how many total bloom peaks over a year vs 2 blooms - is this necessarily much higher than
        #TODO reimplement this to reflect 1/2 bloom data - we can establish this over one year at a time

        if low[3] or high[3]:
            if low[3]:
                n1 = 1 if abs(low[3]) <= (date_seperation_per_year) else 0
            else:
                n1 = 0
            if high[3]:
                n2 = 1 if abs(high[3]) <= (date_seperation_per_year) else 0
            else:
                n2 = 0

            ngd = n1 + n2
        else:
            ngd = 0
        ngds.append(ngd)
    return blooms, ngds

def prepare_sst_variables(sst_array, numpy_storage, skip=False):
    """
    Creates smoothed sst, currently has a large portion commented out as source file is already the centered derivative diff data.
    """
    #smoothed sst
    if not skip:
        sst_boxcar = numpy.apply_along_axis(numpy.convolve, 0, sst_array, numpy.ones((8,))/8, mode='valid')
        fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])), mask=numpy.ones((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])))
        sst_boxcar = numpy.concatenate([fill_arr, fill_arr, fill_arr, sst_boxcar, fill_arr, fill_arr, fill_arr, fill_arr])
        sst_boxcar_map = numpy.memmap(os.path.join(numpy_storage, "sst_sbx"), mode="w+", shape=sst_boxcar.shape, dtype=sst_boxcar.dtype)
        sst_boxcar_map[:] = sst_boxcar[:]
        sst_boxcar = None
        sst_boxcar = sst_boxcar_map
        sst_boxcar_map = None
        sst_array = None
        #get sst derivative
        sst_der = numpy.apply_along_axis(centered_diff_derivative, 0, sst_boxcar)
    else:
        print("Skipped sst data preparation")
        sst_der = sst_array
    sst_der_map = numpy.memmap(os.path.join(numpy_storage, "sst_der"), mode="w+", shape=sst_der.shape, dtype=sst_der.dtype)
    sst_der_map[:] = sst_der[:]
    return sst_der.shape, sst_der.dtype

def prepare_chl_variables(chl_array, numpy_storage, date_seperation, chl_lats, chl_lons, do_model_med=False, median_threshold=1.2):
    """
    Creates the smoothed anomaly chlorophyll data, saves a file to the temporary directory that is read as a mem map later to conserve resources.
    """
    #median * 1.05
    
    #! here i would prefer we output the media value (without adding 5% to it)
    
    #TODO insert the solar zenith angle establishment, only if handed a modelled input (so add flag for --modelled-input)
    if do_model_med:
        days_per_ob = round(365 / date_seperation)
        half_entry = days_per_ob / 2
        ob_dates = [(d * days_per_ob) - half_entry for d in range(1,date_seperation+1)]
        #use 70 degrees as cut off
        date_masks = []
        true_zens = []
        for d in ob_dates:
            date_zeniths = []
            true_zen = []
            for index, lat in enumerate(chl_lats):
                zen = zenithreturn(d, lat)
                if zen < 70:
                    date_zeniths.append(0)
                else:
                    date_zeniths.append(1)
                true_zen.append(zen)
            true_zens.append(true_zen)
            date_masks.append(date_zeniths)
        
        temp_chl_array = chl_array.copy()
        ods = nc.Dataset("lat_zen_angle.nc", "w")
        ods.createDimension('LATITUDE', chl_lats.shape[0])
        ods.createDimension('LONGITUDE', chl_lons.shape[0])
        ods.createDimension('TIME', date_seperation)
        ods.createVariable('LATITUDE', 'float64', dimensions=['LATITUDE'])
        ods.variables['LATITUDE'].setncattr("units", "degrees_north")
        ods.variables['LATITUDE'][:] = chl_lats

        ods.createVariable('LONGITUDE', 'float64', dimensions=['LONGITUDE'])
        ods.variables['LONGITUDE'].setncattr("units", "degrees_north")
        ods.variables['LONGITUDE'][:] = chl_lons

        ods.createVariable('TIME', 'float32', dimensions=['TIME'])
        ods.variables['TIME'].setncattr("units", "years")
        ods.createVariable('zen', 'float32', dimensions=['TIME', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
        ods.variables['zen'].setncattr("units", "degrees")
        for year in range(0, date_seperation, chl_array.shape[0] + 1):
            for index, date_mask in enumerate(date_masks):
                for row, row_mask in enumerate(date_mask):
                    if not row_mask:
                        if LAT_IDX == 2:
                            temp_chl_array.mask[year + index,0,row,:] = True
                        if LAT_IDX == 3:
                            temp_chl_array.mask[year + index,0,:,row] = True
                ods.variables['zen'][index,row,:] = true_zens[index][row]
        ods.close()

        #TODO add gap filling, --gap-filling with a few choices for interpolation options, if not specified then don't do it at all

        print("med5")
        med5 = numpy.ma.median(temp_chl_array,axis = 0)
        print("median value: {}".format(med5))
        print("median threshold: {}".format(median_threshold))
    else:
        print("med5")
        med5 = numpy.ma.median(chl_array,axis = 0)
        print("median value: {}".format(med5))
        print("median threshold: {}".format(median_threshold))

    #! after estimating the median, i apply a filling to the chl time-series, e.g., window of 8 time steps, but it would be good to be able to choose other width of the window, e.g., 5 time-steps or 3 time-steps...

    #! here it would be good to give the option to select the value of the median threshold, e.g., median plus 5%, 10%, 15%...
    
    #get anomalies
    print("anomaly")
    anomaly = chl_array - (med5*median_threshold)
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
    fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])), mask=numpy.ones((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])))
    chl_boxcar = numpy.concatenate([fill_arr, fill_arr, fill_arr, chl_boxcar, fill_arr, fill_arr, fill_arr, fill_arr])
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

    ds.createDimension('LONGITUDE', chl_lons.shape[0])
    ds.createDimension('LATITUDE', chl_lats.shape[0])
    ds.createDimension('DEPTH', output_shape[1])
    ds.createDimension('TIME', None)
    ds.createVariable('LATITUDE', 'float64', dimensions=['LATITUDE'])
    ds.variables['LATITUDE'].setncattr("units", "degrees_north")
    ds.variables['LATITUDE'][:] = chl_lats
    ds.createVariable('LONGITUDE', 'float64', dimensions=['LONGITUDE'])
    ds.variables['LONGITUDE'].setncattr("units", "degrees_east")
    ds.variables['LONGITUDE'][:] = chl_lons
    ds.createVariable('DEPTH', 'float32', dimensions=['DEPTH'])
    ds.variables['DEPTH'].setncattr("units", "meters")
    ds.variables['DEPTH'].setncattr("positive", "down")
    ds.variables['DEPTH'][:] = [0.1]
    ds.createVariable('TIME', 'float32', dimensions=['TIME'])
    ds.variables['TIME'].setncattr("units", "years")
    #switch back to LATITIUDE and LONGITUDE, establish why the flipping of the axis makes everything go screwey
    ds.createVariable('date_start1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_start1'].setncattr("units", "weeks")
    ds.createVariable('date_max1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_max1'].setncattr("units", "weeks")
    ds.createVariable('date_end1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_end1'].setncattr("units", "weeks")
    ds.createVariable('duration1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['duration1'].setncattr("units", "weeks")
    ds.createVariable('date_start2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_start2'].setncattr("units", "weeks")
    ds.createVariable('date_max2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_max2'].setncattr("units", "weeks")
    ds.createVariable('date_end2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_end2'].setncattr("units", "weeks")
    ds.createVariable('duration2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['duration2'].setncattr("units", "weeks")
    ds.createVariable('max_val1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.createVariable('max_val2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['max_val2'].setncattr("units", "mgChl/m3")
    ds.variables['max_val1'].setncattr("units", "mgChl/m3")
    ds.createVariable('total_blooms', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['total_blooms'].setncattr("units", "observations")
    ds.createVariable('probability', 'float32', dimensions=['DEPTH', 'LONGITUDE', 'LATITUDE'],fill_value=FILL_VAL)
    ds.variables['total_blooms'].setncattr("units", "likelihood")
    ds.close()
    print("created netcdf {}".format(name))

def write_to_output_netcdf(data, total_blooms=None, probability=None):
    """
    Loops through each year in the numpy array and writes the data to the netcdf file, this should work faster if we get rid of the loop but I can't seem to grock the logic to fix it right now.
    """
    ds = nc.Dataset(output_location,'r+',format='NETCDF4_CLASSIC')
    data = data.astype(numpy.float32)
    data = numpy.ma.fix_invalid(data)
    print(output_location)
    print("pre-writing data shape: {}".format(data.shape))
    print(data[:,:,:,0,0].shape)
    ds.variables['TIME'][:] = range(0, data.shape[2])
    for year in range(0, data.shape[2]):
        ds.variables['date_start1'][year] = data[:,:,year,0,0]
        ds.variables['date_max1'][year] = data[:,:,year,0,3]
        ds.variables['date_end1'][year] = data[:,:,year,0,1]
        ds.variables['duration1'][year] = data[:,:,year,0,2]
        print(data[:,:,year,1,0])
        ds.variables['date_start2'][year] = data[:,:,year,1,0]
        ds.variables['date_max2'][year] = data[:,:,year,1,3]
        ds.variables['date_end2'][year] = data[:,:,year,1,1]
        ds.variables['duration2'][year] = data[:,:,year,1,2]
        ds.variables['max_val1'][year] = data[:,:,year,0,4]
        ds.variables['max_val2'][year] = data[:,:,year,1,4]
        if total_blooms is not None:
            ds.variables['total_blooms'][year] = total_blooms[:,:,year]
        if probability is not None:
            ds.variables['probability'][:] = probability
    print(ds.variables['TIME'][:])
    ds.close()


def total_blooms_to_probability(array_like):
    return numpy.count_nonzero(array_like == 2) / array_like.size

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
    year_true_start_end_array = numpy.ndarray((chl_data.shape[2],chl_data.shape[3], int((chl_data.shape[0] -start_date) // date_seperation_per_year), 2,5))
    print(int((chl_data.shape[0] -start_date) // date_seperation_per_year))
    total_blooms = numpy.ndarray((chl_data.shape[2],chl_data.shape[3], int((chl_data.shape[0] -start_date) // date_seperation_per_year)))
    year_true_start_end_array.fill(FILL_VAL)
    total_blooms.fill(FILL_VAL)
    completion_points = range(0, chl_data.shape[2], chl_data.shape[2] // 10)
    print("doing sst initiations and correction")
    for ix, iy in numpy.ndindex(chl_data.shape[2], chl_data.shape[3]):
        try:
            verbose=False
            results = match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=False, start_date=start_date)
            year_true_start_end_array[ix,iy] = results[0]
            total_blooms[ix,iy] = results[1]
            if verbose:
                print("end duration array")
                print(year_true_start_end_array[ix,iy])
        except Exception as e:
            print(e)
            print(repr(e))
            print(match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=False, start_date=start_date))
            print()
        if ix in completion_points and iy == 0:
            print(completion_points.index(ix) * 10, "% complete")
    
    probability_array = numpy.apply_along_axis(total_blooms_to_probability, 2, total_blooms)

    print("done sst initiations and correction")
    print("writing to netcdf")
    #needs to be extended to be able to output 3 files: sorted by calendar year, one with primary and secondary chloropjhyll maximum
    write_to_output_netcdf(year_true_start_end_array, total_blooms=total_blooms, probability=probability_array)

    
def extend_array(array_like, entries_per_year, start_date=0):
    """
    Takes the start and end year entries ( assuming it ends on the final index) and adds those to the array given.

    A bit fragile since it makes a few assumptions about the time axis, but we should be fixing that elsewhere in the program
    """
    #take first year
    first_entry = array_like[start_date:start_date + entries_per_year]
    last_entry = array_like[-entries_per_year:]
    #stick em together
    output = numpy.concatenate([first_entry, array_like, last_entry], axis = 0)
    return output, start_date + entries_per_year



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("chl_location", nargs='+', help="A chlorophyll file or list of chlorophyll files. These must have a variable called chl, depth, lat and lon (or contain those strings) which will be used to define the array shape")
    parser.add_argument("--chl_var", help="specify chl variable, otherwise is guessed based on variables that contain'chl' in the name", default="chl", required=False)
    parser.add_argument("--sst_location", help="An sst file, or glob of files (e.g. sst*.nc) that matches the chlorophyll observations. if it does not match some interpolation can be attempted with --sst_date_interp", required=False)
    parser.add_argument("--sst_var", help="specify the sst variable name, otherwise is guessed based on variables containing 'sst'", default="sst", required=False)
    parser.add_argument("--output", help="output filename, if not specified defaults to (chlorophyll_filename)_phenology.nc in the current folder", default=None, required=False)
    parser.add_argument("--date_seperation_per_year", help="how many temporal observations we have in a year, if not specified will be guessed", default=47, required=False)
    parser.add_argument("--first_date_index", help="specify if the first date you want to include is not the first date present in the date stack", default=0, required=False)
    parser.add_argument("--chl_climatology", action="store_true", help="extend the input chlorophyll array by creatingidentical copied for the year previous and year next", default=0, required=False)
    parser.add_argument("--intermediate_file_store", help="change where intermediate numpy files are placed, if not specified then /tmp is assumed - you should specify somewhere else if your tmp cannot handle the array sizes needed (and currently this program will fill it until it cannot).", required=False)
    parser.add_argument("--reverse_search", default=False, help="specify the number of observations to search in the previous year, if not specified will be calculated as a representation of 100 days (date_seperation_per_year / 0.27).", required=False)
    parser.add_argument("--skip_sst_prep", action="store_true", default=False, help="skip sst preparation, instead the program will assume that the sst input is already in a state where low/high period fluctuation can be identified", required=False)
    parser.add_argument("--median_threshold", default=MEDIAN_THRESHOLD_DEFAULT, help="change median threshold, set as percentage value e.g. 20 = 20%", required=False)
    parser.add_argument("--modelled_median", action='store_true',default=False, help="test for solar zenith of inputs", required=False)
    #probably needs a better description!
    #give options for both since we might end up in a situation where sst is 3 years and chl is one year (or vice versa)
    parser.add_argument("--extend_chl_data", default=False, action="store_true", help="extends chlorophyll by copying the (central) chl array to the year previous and year next")
    parser.add_argument("--extend_sst_data", default=False, action="store_true", help="extends sea surfaace temperature by copying the (central) chl array to the year previous and year next")
    parser.add_argument("--reshape", default=False, action="store_true", help="reshape to be t, 1, x, y")
    args = parser.parse_args()
    med_thresh = 1+ (float(args.median_threshold) / 100)
    if not args.intermediate_file_store:
        numpy_storage = tempfile.mkdtemp(prefix="phenology_")
    else:
        numpy_storage = args.intermediate_file_store
    #remember to change median threshold to percentage!
    #reverse search should be ratio of 100 days so 5 days * 20 = 100 or 8 days * 12.5 (so 13) = 100 days
    #as 100 days is 0.27 of 365 we can just do that
    date_seperation_per_year = int(args.date_seperation_per_year)
    if not args.reverse_search:
        reverse_search = int(round(int(args.date_seperation_per_year) * 0.28))
    print("Reverse search:{}".format(reverse_search))

    #TODO list of files or file specified (mid november)
    for chl_location in args.chl_location:
        if args.sst_location:
            print("sst file provided, reading array")
            print("only one file found, assuming full stack of observations")
            sst_ds = nc.Dataset(args.sst_location)
            sst_variable = [x for x in sst_ds.variables if args.sst_var in x.lower()][0]
            sst_lon_var = [x for x in sst_ds.variables if "lon" in x.lower()][0]
            sst_lat_var = [x for x in sst_ds.variables if "lat" in x.lower()][0]
            sst_lons = sst_ds.variables[sst_lon_var][:]
            sst_lats = sst_ds.variables[sst_lat_var][:]
            sst_array = sst_ds.variables[sst_variable][:]
            if args.extend_sst_data:
                sst_array, _ = extend_array(sst_array, date_seperation_per_year, args.first_date_index)
            if args.reshape:
                sst_array.shape = (sst_array.shape[0], 1, sst_array.shape[1], sst_array.shape[2])
            sst_shape, sst_dtype = prepare_sst_variables(sst_array, numpy_storage, skip=args.skip_sst_prep)
            print("sst_shape: {}".format(sst_shape))
            print("sst_dtype: {}".format(sst_dtype))
            sst_array = None

        chl_files = chl_location
        print("calculating on {}".format(chl_files))
        chl_filename = chl_location
        chl_ds = nc.Dataset(chl_location)
        #test for more than 1 variable, if so quit out and complain that it doesn't know which to use
        chl_variable = [x for x in chl_ds.variables if args.chl_var in x.lower()][0]
        chl_lon_var = [x for x in chl_ds.variables if "lon" in x.lower()][0]
        chl_lat_var = [x for x in chl_ds.variables if "lat" in x.lower()][0]
        chl_lons = chl_ds.variables[chl_lon_var][:]
        chl_lats = chl_ds.variables[chl_lat_var][:]
        print("lats shape",chl_lats.shape)
        print("lons shape",chl_lons.shape)
        LAT_IDX = chl_ds.variables[chl_variable].dimensions.index(chl_lat_var)
        LON_IDX = chl_ds.variables[chl_variable].dimensions.index(chl_lon_var)
        if LAT_IDX > LON_IDX:
            DIM_ORDER = ['TIME', 'DEPTH', 'LONGITUDE', 'LATITUDE']
        else:
            DIM_ORDER = ['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE']
        chl_array = chl_ds.variables[chl_variable][:]
        if args.reshape:
            chl_array.shape = (chl_array.shape[0], 1, chl_array.shape[1], chl_array.shape[2])

        #check if there are any nans in the data
        chl_array = numpy.ma.masked_invalid(chl_array)

        """
        if not (chl_array.shape[2] == chl_lats.shape[0] and chl_array.shape[3] == chl_lons.shape[0]):
            print("adjusting to flip lat and lon")
            chl_array.shape = (chl_array.shape[0], chl_array.shape[1], chl_array.shape[3], chl_array.shape[2])
        """
        start_date = None
        if args.extend_chl_data:
            chl_array, start_date = extend_array(chl_array, date_seperation_per_year, args.first_date_index)

        start_date = args.first_date_index if not start_date else start_date

        chl_shape, chl_dtype = prepare_chl_variables(chl_array, numpy_storage, date_seperation_per_year, chl_lats, chl_lons, do_model_med=args.modelled_median, median_threshold=med_thresh)
        print("chl_shape: {}".format(chl_shape))
        print("chl_dtype: {}".format(chl_dtype))

        if sst_shape[2:] != chl_shape[2:]:
            print("sst and chlorophyll x,y array shapes do not match got:")
            print("chlorophyll:",chl_shape[2:])
            print("sst:",sst_shape[2:])
            print("quitting!")
            sys.exit()

        #TODO set dynamically
        sst_dtype = 'float64'
        chl_dtype = 'float64'

        if args.output:
            output = args.output
        else:
            output = chl_filename.replace(".nc", "_phenology.nc")
        
        print("creating output netcdf {}".format(output))
        #simple regridding
        print(chl_shape)
        create_phenology_netcdf(chl_lons, chl_lats, chl_shape, output)

        get_multi_year_two_blooms_output(numpy_storage,
                                        chl_shape,
                                        chl_dtype, 
                                        chl_array, 
                                        sst_shape, 
                                        sst_dtype, 
                                        date_seperation_per_year=date_seperation_per_year, 
                                        start_date=start_date, 
                                        reverse_search=reverse_search)
