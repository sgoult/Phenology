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
import tqdm
import calendar
import datetime
import scipy.misc
import functools
from fast_polarity.polarity import polarity_edge_finder_optimised as polaritiser

#TODO Set dynamically from the input netcdf or user specified from command line
FILL_VAL = -9.999999999999998e+33

#TODO set this from the user inputs
MEDIAN_THRESHOLD_DEFAULT = 20
LAT_IDX = 2
LON_IDX = 3
REF_MONTH = 'January'
END_MONTH = 'December'
USE_DTYPE = "float64"

output_location = None
output_location_date = None

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
    squeezed = array_like.squeeze()
    cf = numpy.convolve(squeezed, [1,0,-1],'same') / 1
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
        zero_crossings = numpy.where(polaritiser(array_like.astype(numpy.float)))[0]
    else:
        zero_crossings = []
    true_poss = zero_crossings

    if verbose:
        print("chl sbx for current period")
        print(array_like)
        print("zero crossings in chl sbx")
        print(zero_crossings)

    
    """
    if not first_real_val in zero_crossings:
        print(first_real_val)
        print(array_like[first_real_val])
        zero_crossings = [first_real_val] + zero_crossings
        zero_crossings.sort()
    """

    #find out which way we're going
    starts = []
    ends = []
    #works much the same qas the SST one below
    #TODO work out why we have dropped the first index (1) from starts - theres no reason it should be doing this!!
    poss_starts = [x for x in zero_crossings if array_like[x] < 0]
    poss_ends = [x for x in zero_crossings if array_like[x] > 0]

    starts = []
    #flip this around!
    already_tested = False
    for idx, date in enumerate(poss_starts):
        if not idx == len(poss_starts) - 1:
            #TODO change this to a command line option
            if already_tested:
                already_tested = False
                continue
            if (poss_starts[idx+1] - date) <= 5:
                starts.append(date)
                already_tested = True
            else:
                starts.append(date)
        else:
            if (date - poss_starts[idx-1]) > 5:
                starts.append(date)
    ends = []
    for start in starts:
        try:
            ends.append(next(x for x in poss_ends if x > start))
        except:
            continue

    #we know the last entry will be an end
    if starts and ends:
        if starts[-1] > ends[-1]:
            ends.append(len(array_like))

    durations = []
    for start in starts:
        end = next(x for x in ends if x > start)
        dura = end - start +1 
        durations.append(dura)
        
            


    if verbose:
        print("starts and ends")
        print("durations:", durations)
        print("starts: ", starts)
        print("ends: ", ends)

    #find an end for every start
    dates = []
    for idx, start in enumerate(starts):
        if durations[idx] <= 2:
            continue
        try:
           end = next(x for x in ends if x > start)
           max_idx = numpy.nanargmax(chl_values[start:end + 1])
           dates.append([start + date_offset,end + date_offset,end-start,max_idx + date_offset + start,chl_values[max_idx + date_offset + start]])
           max_idx = None
        except ValueError:
            if chl_values[start:end + 1]:
                print("val error encountered at max estimation, sequence was:")
                print(chl_values[start:end + 1])
                print("start end")
                print(start, end)
            continue
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

        print("maxes")
        print([x[4] for x in dates])
    
        print("tmaxes")
        print([x[3] for x in dates])
    return dates

def phen_records_to_one_val_on_max(records, date_correction=False, index=4,verbose=False):
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

def polarity_edge_finder(input_array):
    return (numpy.diff(numpy.sign(input_array)) != 0)*1

def match_start_end_to_solar_cycle(array_like, chl_sbx_slice, chl_slice, date_seperation_per_year, reverse_search, start_date=0, reference_date=0, verbose=False):
    """
    Attributes the start and end times in relation to the SST or solar cycle, takes an sst array (array_like), smoothed chlorophyll derivative slice (chl_sbx_slice) and the original chlorophyll data.

    Slices up the data based on high/low periods of SST (or otherwise), then feeds each period into get_start_index_and_duration, once finished it will output an array of shape (x, y, time, 2, 5)
    verbose will spame the terminal with information about whats going on, best to establish a few pixels you want to inspect rather than having this on all the time.

    in a run using global data that too 30 minutes, this function made up 703 seconds of the processing time
    I would guess that 500 of those seconds can be attributed to get_start_index_and_duration
    """ 
    if array_like.mask.all():
        if verbose:
            print("array all -32616.30273438")
            print(array_like)
        #we can stop here, there's not point continuing with an empty array
        return [[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]]

    #possibly resort and create new durations based on remaining dates
    #look for sign changes in sst or PAR data, indicating high/low SST or light periods
    array_like = numpy.squeeze(array_like)
    chl_slice = numpy.squeeze(chl_slice)
    chl_sbx_slice = numpy.squeeze(chl_sbx_slice)
    zero_crossings = numpy.where(polaritiser(array_like.astype(numpy.float)))[0]

    if verbose:
        print("sst values for pixel")
        print(array_like)
        print("zero crossings of sst for pixel")
        print(zero_crossings)

    #find out which way we're going
    highs = []
    lows = []
    highs = [x for x in zero_crossings if array_like[x] < 0]
    lows = [x for x in zero_crossings if array_like[x] > 0]
    
    #print(everything thus far
    if verbose:
        print("***********************")
        print("starts of high and low periods")
        print("highs: ", highs)
        print("lows: ", lows)
    maximum_sst = []
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
            if highs[-1] < lows[-1]:
                highs.append(len(array_like))
            else:
                lows.append(len(array_like))
        elif len(highs) and not len(lows):
            lows = [0, len(array_like)]
        elif len(lows) and not len(highs):
            highs = [0, len(array_like)]
        else:
            return [[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]]
    except Exception as e:
        #triggered once, but was helpful to know what the contents were
        print(highs)
        print(lows)
        print(e)

    #chuck out the results
    if verbose:
        print("updated highs and lows after filtering and correcting for missing dates")
        print("highs: ", highs)
        print("lows: ", lows)

    period_chl_phenology = get_start_index_and_duration(chl_sbx_slice,chl_slice,0,depth=5,verbose=verbose)

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
        #pad = 0
        if verbose:
            print("running bloom detection for period: ", index, end_date)
            print("this is a ", "high" if activity_period else "low", "period")

        if activity_period:
            high_records.extend([x for x in period_chl_phenology if x[3] > index and x[3] < end_date])
        else:
            low_records.extend([x for x in period_chl_phenology if x[3] > index and x[3] < end_date])


        """
        old stuff
        end_date = int(round(end_date + pad))
        start = index - pad if index >= pad else index
        if verbose:
            print("padding period with ", pad, "entries")
            print("resultant dates (start, end): ", start, end_date)
        chl_sbx_period_slice = chl_sbx_slice[start:end_date]
        chl_period_slice = chl_slice[start:end_date]
        #get the phenology for this period, depth pads extra data if needed for numpy (we don't use this for SST model)
        #73% of time in function
        #run this before then sort into high/low periods
        period_chl_phenology = get_start_index_and_duration(chl_sbx_period_slice,chl_period_slice,start,depth=5,verbose=verbose)
        #if we found anything
        if len(period_chl_phenology):
            #loop through them and add them to the high/low mega lists
            for record in period_chl_phenology:
                if activity_period:
                    high_records.append(record)
                else:
                    low_records.append(record)
        """
    
    #to remind ourselves what the phenology records look like
    #[start,end,end-start,max_idx,chl_values[max_idx]]
    blooms = []
    ngds = []
    ngds_date = []
    total_blooms = []
    blooms_by_date = []
    for year in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year):
        if year + date_seperation_per_year > chl_sbx_slice.shape[0]:
            continue
        #find blooms that start after the year - our reverse search, end before the end of the year, and end during the current year
        #this is where I have concerns that there are problems with the selection of blooms, if we're doing a year with reference month of June
        #then this currently selects june - june, rather than january - december, is this correct? The central month 
        possible_high_blooms = [x for x in high_records if x[3] > (year) and x[3] < (year + date_seperation_per_year) and not x[0] < (year - date_seperation_per_year)]
        possible_low_blooms = [x for x in low_records if x[3] > (year) and x[3] < (year + date_seperation_per_year) and not x[0] < (year - date_seperation_per_year)]

        #filters out the blooms that might overlap
        high_removals = []
        low_removals = []
         #filters out the blooms that might overlap
        for hindex, high_bloom in enumerate(possible_high_blooms):
            append_high = True
            for lindex, low_bloom in enumerate(possible_low_blooms):
                if any([low_bloom[x] == high_bloom[x] for x in [0, 1]]):
                    #whichever max is higher we should select
                    if low_bloom[4] > high_bloom[4]:
                        high_removals.append(hindex)
                    else:
                        low_removals.append(lindex)
        
        for idx in sorted(set(high_removals), reverse = True):
            del possible_high_blooms[idx]
        
        for idx in sorted(set(low_removals), reverse = True):
            del possible_low_blooms[idx]

        high_blooms = possible_high_blooms    
        low_blooms = possible_low_blooms

        if verbose:
            print("working on year: ", year)
            print("found ", len(possible_high_blooms), " high blooms")
            print("found ", len(possible_low_blooms), " low blooms")
        
        if verbose:
            print("highs")
            print(high_blooms)
            print("lows")
            print(low_blooms)
        #reduce them to one record
        #additionally look at using max vs duration for the bloom selection
        low = phen_records_to_one_val_on_max(low_blooms, year + reference_date)
        high = phen_records_to_one_val_on_max(high_blooms, year + reference_date)
        if not low or not high:
            print("***************")
            print(high_records)
            print(low_records)
            print(year - reverse_search)
            print(possible_high_blooms)
            print(possible_low_blooms)
            print(low, high)
        #spit out the low period and high period for this year
        #
        if verbose:
            print("high bloom")
            print(high)
            print("low bloom")
            print(low)
        high_val = high[4] if high[4] else -1000
        low_val = low[4] if low[4] else -1000
        if low_val > high_val:
            blooms.append([low,high])
        else:
            blooms.append([high,low])

        
        #Alternative date sorting, based purely on high/low max date
        high_date = high[3] if high[3] else -1000
        low_date = low[3] if low[3] else -1000
        if low_date > high_date:
            blooms_by_date.append([high,low])
        else:
            blooms_by_date.append([low,high])
        

        """
        date_sorted = [None, None]
        if low[3]:
            if low[3] > (date_seperation_per_year // 2):
                date_sorted[1] = low
            else:
                date_sorted[0] = low
        
        #change here to change date sorting


        if high[3]:
            if high[3] > (date_seperation_per_year // 2):
                if high and (not date_sorted[1] or high[4] > low[4]):
                    date_sorted[1] = high
            else:
                if high and (not date_sorted[0] or high[4] > low[4]):
                    date_sorted[0] = high
        
        if not date_sorted[0]:
            date_sorted[0] = [None,None,None,None,None]
        if not date_sorted[1]:
            date_sorted[1] = [None,None,None,None,None]
        

        blooms_by_date.append(date_sorted)
        """
        ngd_date = int(not low[0] is None) + int(not high[0] is None)
        ngds_date.append(ngd_date)
        
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
        total_blooms.append(len(possible_high_blooms) + len(possible_low_blooms))
    return blooms, ngds, blooms_by_date, ngds_date, total_blooms

def prepare_sst_variables(sst_array, numpy_storage, skip=False):
    """
    Creates smoothed sst, currently has a large portion commented out as source file is already the centered derivative diff data.
    """
    #smoothed sst
    if not skip:
        print("smoothing sst")
        #print(sst_array[:20,:,100,281])
        print(sst_array[:5,:,100,281])
        sst_boxcar = numpy.apply_along_axis(lambda m: numpy.convolve(m, numpy.ones(8)/8, mode='valid'), axis=0, arr=sst_array)
        print("shape after sst sbx")
        print(sst_boxcar.shape)
        fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])), mask=numpy.ones((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])))
        sst_boxcar_map = numpy.memmap(os.path.join(numpy_storage, "sst_sbx"), mode="w+", shape=sst_boxcar.shape, dtype=USE_DTYPE)
        sst_boxcar_map[:] = sst_boxcar[:]
        sst_boxcar = None
        sst_boxcar = sst_boxcar_map
        sst_boxcar_map = None
        sst_array = None
        #get sst derivative
        print("doing sst_derivative")
        sst_der = numpy.apply_along_axis(centered_diff_derivative, 0, sst_boxcar[5:,:,:,:])
        sst_der = numpy.concatenate([fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, sst_der, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr,])
        print("shape after sst der")
        print(sst_der.shape)
    else:
        print("Skipped sst data preparation")
        sst_der = sst_array
        numpy.ma.set_fill_value(sst_array, fill_value=FILL_VAL)
        print(sst_array.shape)
    sst_der_map = numpy.memmap(os.path.join(numpy_storage, "sst_der"), mode="w+", shape=sst_der.shape, dtype=USE_DTYPE)
    print(sst_der)
    sst_der_map[:] = sst_der[:]
    return sst_der.shape, 'float64'

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
        #use 75 degrees as cut off
        date_masks = []
        true_zens = []
        for d in ob_dates:
            date_zeniths = []
            true_zen = []
            for index, lat in enumerate(chl_lats):
                zen = zenithreturn(d, lat)
                if zen <= 75:
                    date_zeniths.append(0)
                else:
                    date_zeniths.append(1)
                true_zen.append(zen)
            true_zens.append(true_zen)
            date_masks.append(date_zeniths)
        
        temp_chl_array = chl_array.copy()
        ods = nc.Dataset("lat_zen_angle_{}.nc".format(datetime.datetime.now().strftime("%H%M")), "w")
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
        temp_chl_array
        med5 = numpy.ma.median(temp_chl_array,axis=0)
        print("median value: {}".format(med5))
        print("median threshold: {}".format(median_threshold))
    else:
        print("med5")
        med5 = numpy.ma.median(chl_array,axis = 0)
        print("median value: {}".format(med5))
        print("median threshold: {}".format(median_threshold))

    ods = nc.Dataset(os.path.join(numpy_storage,"median_output_{}.nc".format(datetime.datetime.now().strftime("%H%M"))), "w")
    ods.createDimension('LATITUDE', chl_lats.shape[0])
    ods.createDimension('LONGITUDE', chl_lons.shape[0])
    ods.createDimension('TIME', 1)
    ods.createVariable('LATITUDE', 'float64', dimensions=['LATITUDE'])
    ods.variables['LATITUDE'].setncattr("units", "degrees_north")
    ods.variables['LATITUDE'][:] = chl_lats

    ods.createVariable('LONGITUDE', 'float64', dimensions=['LONGITUDE'])
    ods.variables['LONGITUDE'].setncattr("units", "degrees_north")
    ods.variables['LONGITUDE'][:] = chl_lons

    ods.createVariable('TIME', 'float32', dimensions=['TIME'])
    ods.variables['TIME'].setncattr("units", "years")
    ods.createVariable('median', 'float32', dimensions=['TIME', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL)
    ods.variables['median'].setncattr("units", "mg chl m^3")
    ods.variables['median'][:] = med5

    #! after estimating the median, i apply a filling to the chl time-series, e.g., window of 8 time steps, but it would be good to be able to choose other width of the window, e.g., 5 time-steps or 3 time-steps...

    #! here it would be good to give the option to select the value of the median threshold, e.g., median plus 5%, 10%, 15%...
    
    #get anomalies
    print("anomaly")
    anomaly = chl_array - (med5*median_threshold)
    print(anomaly[:5,:,100,281])
    anomaly_map = numpy.memmap(os.path.join(numpy_storage, "chl_anomaly"), mode="w+", shape=anomaly.shape, dtype=USE_DTYPE)
    anomaly_map[:] = anomaly[:]
    #anomaly = anomaly_map
    chl_array = None
    anomaly_map = None
    #need to ditch any empty entries here as they interefere with the cumsum

    #get cumsum of anomalies
    print("chl cumsum")
    chl_cumsum = numpy.ma.cumsum(anomaly,axis=0)
    print(chl_cumsum[:5,:,100,281])
    chl_cumsum_map = numpy.memmap(os.path.join(numpy_storage, "chl_cumsum"), mode="w+", shape=chl_cumsum.shape, dtype=USE_DTYPE)
    chl_cumsum_map[:] = chl_cumsum[:]
    #chl_cumsum = chl_cumsum_map
    anomaly = None
    chl_cumsum_map = None

    #get centered derivative        
    print("chl der")
    chl_der = numpy.apply_along_axis(centered_diff_derivative, 0, chl_cumsum)
    print(chl_der[:5,:,100,281])
    chl_der_map = numpy.memmap(os.path.join(numpy_storage, "chl_der"), mode="w+", shape=chl_der.shape, dtype=USE_DTYPE)
    chl_der_map[:] = chl_der[:]
    #chl_der = chl_der_map
    chl_cumsum = None
    chl_der_map = None
    

    #boxcar filter with width of 3 (sbx) should be something like this:
    print("chl sbx")
    chl_boxcar = numpy.apply_along_axis(lambda m: numpy.convolve(m, numpy.ones(3)/3, mode='full'), axis=0, arr=chl_der)
    print(chl_boxcar[:5,:,100,281])
    fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])), mask=numpy.ones((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])))
    print("chl shape after boxcar")
    print(chl_boxcar.shape)
    chl_boxcar = numpy.concatenate([chl_boxcar, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr,])
    print("chl shape after boxcar padding")
    print(chl_boxcar.shape)
    chl_boxcar_map = numpy.memmap(os.path.join(numpy_storage, "chl_sbx"), mode="w+", shape=chl_boxcar.shape, dtype=USE_DTYPE)
    chl_boxcar_map[:] = chl_boxcar[:]
    chl_boxcar = None
    chl_boxcar = chl_boxcar_map
    chl_boxcar_map = None

    return chl_boxcar.shape, USE_DTYPE

def create_phenology_netcdf(chl_lons, chl_lats, output_shape=None,name="phenology_{}.nc".format(datetime.datetime.now().strftime("%H%M")), date=False):
    """
    Creates the skeleton of the netcdf file to be used by write_to_output_netcdf, all of this is metadata.
    """
    if date:
        global output_location_date
        output_location_date = name
    else:
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
    week_descriptor = 'weeks from {}'.format(REF_MONTH)
    #description = ' data between {} and {} for years {} to {}'.format(REF_MONTH, END_MONTH, START_YEAR, END_YEAR)
    ds.variables['date_start1'].setncattr("units", week_descriptor)
    ds.createVariable('date_max1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_max1'].setncattr("units", week_descriptor)
    ds.createVariable('date_end1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_end1'].setncattr("units", week_descriptor)
    ds.createVariable('duration1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['duration1'].setncattr("units", week_descriptor)
    ds.createVariable('date_start2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_start2'].setncattr("units", week_descriptor)
    ds.createVariable('date_max2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_max2'].setncattr("units", week_descriptor)
    ds.createVariable('date_end2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['date_end2'].setncattr("units", week_descriptor)
    ds.createVariable('duration2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['duration2'].setncattr("units", week_descriptor)
    ds.createVariable('max_val1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.createVariable('max_val2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['max_val2'].setncattr("units", "mgChl/m3")
    ds.variables['max_val1'].setncattr("units", "mgChl/m3")
    ds.createVariable('total_blooms', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL)
    ds.variables['total_blooms'].setncattr("units", "observations")
    ds.createVariable('probability', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL)
    ds.variables['total_blooms'].setncattr("units", "likelihood")
    ds.close()
    print("created netcdf {}".format(name))

def write_to_output_netcdf(data, total_blooms=None, probability=None, date=False):
    """
    Loops through each year in the numpy array and writes the data to the netcdf file, this should work faster if we get rid of the loop but I can't seem to grock the logic to fix it right now.
    """
    if date:
        ds = nc.Dataset(output_location_date,'r+',format='NETCDF4_CLASSIC')
        print("writing to:", output_location_date)
    else:
        ds = nc.Dataset(output_location,'r+',format='NETCDF4_CLASSIC')
        print("writing to:", output_location)
    data = data.astype(numpy.float32)
    data = numpy.ma.fix_invalid(data)
    print("pre-writing data shape: {}".format(data.shape))
    ds.variables['TIME'][:] = range(0, data.shape[2])
    for year in range(0, data.shape[2]):
        ds.variables['date_start1'][year] = data[:,:,year,0,0]
        ds.variables['date_max1'][year] = data[:,:,year,0,3]
        ds.variables['date_end1'][year] = data[:,:,year,0,1]
        ds.variables['duration1'][year] = data[:,:,year,0,2]
        ds.variables['max_val1'][year] = data[:,:,year,0,4]
        ds.variables['date_start2'][year] = data[:,:,year,1,0]
        ds.variables['date_max2'][year] = data[:,:,year,1,3]
        ds.variables['date_end2'][year] = data[:,:,year,1,1]
        ds.variables['duration2'][year] = data[:,:,year,1,2]
        ds.variables['max_val2'][year] = data[:,:,year,1,4]
        if total_blooms is not None:
            ds.variables['total_blooms'][year] = total_blooms[:,:,year]
        if probability is not None:
            ds.variables['probability'][:] = probability
    print("wrote ", len(ds.variables['TIME'][:]), "years of data")
    ds.close()


def total_blooms_to_probability(array_like):
    return numpy.count_nonzero(array_like == 2) / array_like.size

def get_multi_year_two_blooms_output(numpy_storage, chl_shape, chl_dtype, chl_data, sst_shape, sst_dtype, date_seperation_per_year=47, start_date=0, reverse_search=20, reference_index=0):
    #this all works on the assumption the axis 0 is time
    print("reading variables")
    chl_boxcar = numpy.memmap(os.path.join(numpy_storage, "chl_sbx"), mode="r", dtype=USE_DTYPE, shape=chl_shape)
    chl_boxcar = chl_boxcar.copy()
    print(sst_dtype)
    sst_der = numpy.memmap(os.path.join(numpy_storage, "sst_der"), mode="r", dtype=USE_DTYPE, shape=sst_shape)
    sst_der = sst_der.copy()
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
    blooms_by_date = year_true_start_end_array.copy()
    total_blooms.fill(FILL_VAL)
    total_blooms_date = total_blooms.copy()
    print("doing sst initiations and correction")
    for ix, iy in tqdm.tqdm(numpy.ndindex(chl_data.shape[2], chl_data.shape[3]), total=(chl_data.shape[2] * chl_data.shape[3])):
        try:
            verbose=False
            if iy == args.debug_pixel[0] and ix == args.debug_pixel[1]:
                print(ix,iy)
                verbose = True
            else:
                pass
            results = match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=verbose, start_date=start_date, reference_date=reference_index)
            year_true_start_end_array[ix,iy] = results[0]
            total_blooms[ix,iy] = results[1]
            blooms_by_date[ix, iy] = results[2]
            total_blooms_date[ix,iy] = results[3]
            if verbose:
                print("end duration array")
                print(year_true_start_end_array[ix,iy])
        except Exception as e:
            print(e)
            print(repr(e))
            print(match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=False, start_date=start_date, reference_date=reference_index))
            print()
        """
        if ix in completion_points and iy == 0:
            print(completion_points.index(ix) * 10, "% complete")
        """
    
    print(total_blooms.shape)
    probability_array = numpy.apply_along_axis(total_blooms_to_probability, 2, total_blooms)
    probability_array_date = numpy.apply_along_axis(total_blooms_to_probability, 2, total_blooms_date)

    print("done sst initiations and correction")
    print("writing to netcdf")
    #needs to be extended to be able to output 3 files: sorted by calendar year, one with primary and secondary chloropjhyll maximum
    write_to_output_netcdf(year_true_start_end_array, total_blooms=total_blooms, probability=probability_array)
    write_to_output_netcdf(blooms_by_date, total_blooms=total_blooms_date, probability=probability_array_date, date=True)

    
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
    parser.add_argument("--first_date_index", help="Specify if the first date you want to include is not the first date present in the date stack.", default=1, required=False)
    parser.add_argument("--reference_index", help="Date index in relation to first_date_index to use as the reference from which weeks are measured in the phenology output. If specified is used as first_date_index + reference_index, if not set will measure from first_date_index (default 1)", default=1, required=False)
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
    parser.add_argument("--start_year", default=0, help="What year to use as the start point for metadata output, if not specified will use 0, affects no processing.")
    parser.add_argument("--debug_pixel",  nargs='+', default=[143,112], type=int, help="pixle in x, y (lat, lon), these entries are 0 indexed.")
    parser.add_argument("--reshape", default=False, action="store_true", help="reshape to be t, 1, x, y")
    parser.add_argument("--reshape_sst", default=False, action="store_true", help="reshape to be t, 1, x, y")
    args = parser.parse_args()
    med_thresh = 1+ (float(args.median_threshold) / 100)
    if not args.intermediate_file_store:
        numpy_storage = tempfile.mkdtemp(prefix="phenology_")
    else:
        numpy_storage = args.intermediate_file_store
        os.mkdir(numpy_storage)
    #remember to change median threshold to percentage!
    #reverse search should be ratio of 100 days so 5 days * 20 = 100 or 8 days * 12.5 (so 13) = 100 days
    #as 100 days is 0.27 of 365 we can just do that
    date_seperation_per_year = int(args.date_seperation_per_year)
    if not args.reverse_search:
        reverse_search = int(round(int(args.date_seperation_per_year) * 0.28))
    print("Reverse search:{}".format(reverse_search))

    start_date = None
    start_date = int(args.first_date_index) - 1 if not start_date else start_date
    ref_index = int(args.reference_index) - 1

    if (ref_index + start_date) > (date_seperation_per_year + start_date):
        print("reference index is too great and would result in mostly or all negative values.")
        sys.exit()

    mon_index = int((start_date + ref_index) // (date_seperation_per_year / 12))
    print(mon_index)
    REF_MONTH = calendar.month_name[mon_index]
    print("reference month (central):", REF_MONTH)
    START_YEAR = args.start_year

    #TODO list of files or file specified (mid november)
    for chl_location in args.chl_location:
        if args.sst_location:
            print("sst file provided, reading array")
            print("only one file found, assuming full stack of observations")
            sst_ds = nc.Dataset(args.sst_location)
            sst_variable = [x for x in sst_ds.variables if args.sst_var in x.lower()][0]
            try:
                sst_lon_var = [x for x in sst_ds.variables if "lon" in x.lower()][0]
                sst_lat_var = [x for x in sst_ds.variables if "lat" in x.lower()][0]
            except:
                print("trying to get lat/lon from standard name")
                for va in sst_ds.variables:
                    if "standard_name" in sst_ds.variables[va].__dict__.keys():
                        if sst_ds.variables[va].__dict__["standard_name"] == "latitude":
                            sst_lat_var= va
                        elif sst_ds.variables[va].__dict__["standard_name"] == "longitude":
                            sst_lon_var = va
            sst_lons = sst_ds.variables[sst_lon_var][:]
            sst_lats = sst_ds.variables[sst_lat_var][:]
            sst_array = sst_ds.variables[sst_variable][:]
            if args.extend_sst_data:
                sst_array, _ = extend_array(sst_array, date_seperation_per_year, start_date)
            if args.reshape_sst:
                sst_array.shape = (sst_array.shape[0], 1, sst_lats.shape[0], sst_lons.shape[0])
            sst_shape, sst_dtype = prepare_sst_variables(sst_array, numpy_storage, skip=args.skip_sst_prep)
            print("sst_shape: {}".format(sst_shape))
            print("sst_dtype: {}".format(USE_DTYPE))
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

        print("using debug pixel:", args.debug_pixel, " which equates to:", chl_lats[args.debug_pixel[1]], "N", chl_lons[args.debug_pixel[0]], "E", "(zero indexed so you may need to add 1 to get reference in other software)")

        """
        if not (chl_array.shape[2] == chl_lats.shape[0] and chl_array.shape[3] == chl_lons.shape[0]):
            print("adjusting to flip lat and lon")
            chl_array.shape = (chl_array.shape[0], chl_array.shape[1], chl_array.shape[3], chl_array.shape[2])
        """
        if args.extend_chl_data:
            chl_array, start_date = extend_array(chl_array, date_seperation_per_year, start_date)

        chl_shape, chl_dtype = prepare_chl_variables(chl_array, numpy_storage, date_seperation_per_year, chl_lats, chl_lons, do_model_med=args.modelled_median, median_threshold=med_thresh)
        print("chl_shape: {}".format(chl_shape))
        print("chl_dtype: {}".format(USE_DTYPE))

        if sst_shape[2:] != chl_shape[2:]:
            print("sst and chlorophyll x,y array shapes do not match got:")
            print("chlorophyll:",chl_shape[2:])
            print("sst:",sst_shape[2:])
            print("quitting!")
            sys.exit()

        #TODO set dynamically

        if args.output:
            output = args.output
        else:
            output = chl_filename.replace(".nc", "_phenology_{}.nc".format(datetime.datetime.now().strftime("%H%M")))
        
        print("creating output netcdf {}".format(output))
        #simple regridding
        create_phenology_netcdf(chl_lons, chl_lats, chl_shape, output.replace(".nc", "_by_maxval.nc"))
        create_phenology_netcdf(chl_lons, chl_lats, chl_shape, output.replace(".nc", "_by_date.nc"), date=True)
        print("using start date",start_date)
        get_multi_year_two_blooms_output(numpy_storage,
                                        chl_shape,
                                        chl_dtype, 
                                        chl_array, 
                                        sst_shape, 
                                        sst_dtype, 
                                        date_seperation_per_year=date_seperation_per_year, 
                                        start_date=start_date, 
                                        reverse_search=reverse_search,
                                        reference_index=ref_index)
