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
import multiprocessing
import logging
import time
import traceback
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
median_output_name = None

logging.basicConfig()
logger = logging.getLogger("phenology_logger")
logger.setLevel(logging.DEBUG)



numpy_storage=None

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

    logger.debug("chl sbx for current period")
    logger.debug(array_like)
    logger.debug("zero crossings in chl sbx")
    logger.debug(zero_crossings)

    
    """
    if not first_real_val in zero_crossings:
        logger.debug(first_real_val)
        logger.debug(array_like[first_real_val])
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
    try:
        for start in starts:
            end = next(x for x in ends if x > start)
            dura = end - start +1 
            durations.append(dura)
    except StopIteration:
        pass
        
            



    logger.debug("starts and ends")
    logger.debug("durations: {}".format(durations))
    logger.debug("starts: {}".format(starts))
    logger.debug("ends: {}".format(ends))

    #find an end for every start
    dates = []
    for idx, start in enumerate(starts):
        try:
            if durations[idx] <= 2:
                continue
        except IndexError:
            continue
        try:
           end = next(x for x in ends if x > start)
           max_idx = numpy.nanargmax(chl_values[start:end + 1])
           chl_val = chl_values[max_idx + date_offset + start] if not chl_values.mask[max_idx + date_offset + start] else numpy.nan
           dates.append([start + date_offset,end + date_offset,end-start,max_idx + date_offset + start,chl_val])
           max_idx = None
        except ValueError:
            if chl_values[start:end + 1]:
                logger.debug("val error encountered at max estimation, sequence was:")
                logger.debug(chl_values[start:end + 1])
                logger.debug("start end")
                logger.debug(start, end)
            continue
        except StopIteration:
            continue
        except Exception as e:
           logger.error(e)
           logger.error(repr(e))
           continue
    if pad_values:
        for pad in range(len(dates), depth):
            dates.append([None,None,None,None,None])

    logger.debug("end dates")
    logger.debug(dates)

    logger.debug("maxes")
    logger.debug([x[4] for x in dates])

    logger.debug("tmaxes")
    logger.debug([x[3] for x in dates])
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
        logger.info("array all -32616.30273438")
        logger.info(array_like)
        #we can stop here, there's not point continuing with an empty array
        return [[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]],[[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]]

    #possibly resort and create new durations based on remaining dates
    #look for sign changes in sst or PAR data, indicating high/low SST or light periods
    array_like = numpy.squeeze(array_like)
    chl_slice = numpy.squeeze(chl_slice)
    chl_sbx_slice = numpy.squeeze(chl_sbx_slice)
    zero_crossings = numpy.where(polaritiser(array_like.astype(numpy.float)))[0]

    logger.debug("sst values for pixel")
    logger.debug(array_like)
    logger.debug("zero crossings of sst for pixel")
    logger.debug(zero_crossings)

    #find out which way we're going
    highs = []
    lows = []
    highs = [x for x in zero_crossings if array_like[x] < 0]
    lows = [x for x in zero_crossings if array_like[x] > 0]
    
    #logger.debug(everything thus far
    logger.debug("***********************")
    logger.debug("starts of high and low periods")
    if highs:
        logger.debug("highs: {}".format(*highs))
    else:
        logger.debug("no highs found!")
    if lows:
        logger.debug("lows: {}".format(*lows))
    else:
        logger.debug("no lows found!")
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
            return [[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]],[[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]]
    except Exception as e:
        #triggered once, but was helpful to know what the contents were
        logger.error(highs)
        logger.error(lows)
        logger.error(e)

    #chuck out the results
    logger.debug("updated highs and lows after filtering and correcting for missing dates")
    logger.debug("highs: {}".format(*highs))
    logger.debug("lows: {}".format(*lows))

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
        logger.debug("running bloom detection for period: {} {}".format(index, end_date))
        logger.debug("this is a {}".format("high" if activity_period else "low", "period"))

        if activity_period:
            high_records.extend([x for x in period_chl_phenology if x[3] > index and x[3] < end_date])
        else:
            low_records.extend([x for x in period_chl_phenology if x[3] > index and x[3] < end_date])
        
    logger.debug(high_records)
    logger.debug(low_records)
    
    #to remind ourselves what the phenology records look like
    #[start,end,end-start,max_idx,chl_values[max_idx]]
    blooms = []
    ngds = []
    ngds_date = []
    total_blooms = []
    blooms_by_date = []
    logger.debug("performing year filtering from {}".format(start_date))
    for year in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year):
        logger.debug("doing year {}".format(year))
        if year + date_seperation_per_year > chl_sbx_slice.shape[0]:
            continue
        #find blooms that start after the year - our reverse search, end before the end of the year, and end during the current year
        #this is where I have concerns that there are problems with the selection of blooms, if we're doing a year with reference month of June
        #then this currently selects june - june, rather than january - december, is this correct? The central month 
        possible_high_blooms = [x for x in high_records if x[3] > (year) and x[3] < (year + date_seperation_per_year) and not x[0] < (year - date_seperation_per_year)]
        possible_low_blooms = [x for x in low_records if x[3] > (year) and x[3] < (year + date_seperation_per_year) and not x[0] < (year - date_seperation_per_year)]

        logger.debug("possible_high_blooms")
        logger.debug(possible_high_blooms)
        logger.debug("possible_low_blooms")
        logger.debug(possible_low_blooms)

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

        logger.debug("after filtering")
        logger.debug("possible_high_blooms")
        logger.debug(possible_high_blooms)
        logger.debug("possible_low_blooms")
        logger.debug(possible_low_blooms)

        logger.debug("working on year: {}".format(year))
        logger.debug("found {} high blooms".format(len(possible_high_blooms)))
        logger.debug("found {} low blooms".format(len(possible_low_blooms)))
        
        logger.debug("highs")
        logger.debug(high_blooms)
        logger.debug("lows")
        logger.debug(low_blooms)
        #reduce them to one record
        #additionally look at using max vs duration for the bloom selection
        low = phen_records_to_one_val_on_max(low_blooms, year + reference_date)
        high = phen_records_to_one_val_on_max(high_blooms, year + reference_date)
        if not low or not high:
            logger.debug("***************")
            logger.debug(high_records)
            logger.debug(low_records)
            logger.debug(year - reverse_search)
            logger.debug(possible_high_blooms)
            logger.debug(possible_low_blooms)
            logger.debug(low, high)
        #spit out the low period and high period for this year
        #
        logger.debug("high bloom")
        logger.debug(high)
        logger.debug("low bloom")
        logger.debug(low)
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
    logger.debug((blooms, ngds, blooms_by_date, ngds_date, total_blooms))

    [[None,None,None,None,None], [None,None,None,None,None]],
    [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]],
    [[None,None,None,None,None], [None,None,None,None,None]], 
    [None for x in range(start_date, chl_sbx_slice.shape[0], date_seperation_per_year) if not x + date_seperation_per_year > chl_sbx_slice.shape[0]]
    return blooms, ngds, blooms_by_date, ngds_date, total_blooms

def prepare_sst_variables(sst_array, numpy_storage, chunk, skip=False):
    """
    Creates smoothed sst, currently has a large portion commented out as source file is already the centered derivative diff data.
    """
    global debug_pixel
    #smoothed sst
    if not skip:
        logger.info("smoothing sst")
        #logger.info(sst_array[:20,:,debug_pixel[0],debug_pixel[1]])
        logger.info(sst_array[:5,:,debug_pixel[0],debug_pixel[1]])
        sst_boxcar = numpy.apply_along_axis(lambda m: numpy.convolve(m, numpy.ones(8)/8, mode='valid'), axis=0, arr=sst_array)
        logger.info("shape after sst sbx")
        logger.info(sst_boxcar.shape)
        fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])), mask=numpy.ones((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])))
        sst_boxcar_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "sst_sbx"), mode="w+", shape=sst_boxcar.shape, dtype=USE_DTYPE)
        sst_boxcar_map[:] = sst_boxcar[:]
        sst_boxcar = None
        sst_boxcar = sst_boxcar_map
        sst_boxcar_map = None
        sst_array = None
        #get sst derivative
        logger.info("doing sst_derivative")
        sst_der = numpy.apply_along_axis(centered_diff_derivative, 0, sst_boxcar[5:,:,:,:])
        sst_der = numpy.concatenate([fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, sst_der, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr,])
        logger.info("shape after sst der")
        logger.info(sst_der.shape)
    else:
        logger.info("Skipped sst data preparation")
        sst_der = sst_array
        numpy.ma.set_fill_value(sst_array, fill_value=FILL_VAL)
        logger.info(sst_array.shape)
    sst_der_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "sst_der"), mode="w+", shape=sst_der.shape, dtype=USE_DTYPE)
    logger.info(sst_der)
    sst_der_map[:] = sst_der[:]
    return sst_der.shape, 'float64'

def prepare_chl_variables(chl_array, numpy_storage, chunk, date_seperation, chl_lats, chl_lons, do_model_med=False, median_threshold=1.2, median_filename="median_output.nc"):
    """
    Creates the smoothed anomaly chlorophyll data, saves a file to the temporary directory that is read as a mem map later to conserve resources.
    """
    #median * 1.05
    global debug_pixel
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

        logger.info("med5")
        temp_chl_array
        med5 = numpy.ma.median(temp_chl_array,axis=0)
        logger.info("median value: {}".format(med5))
        logger.info("median threshold: {}".format(median_threshold))
    else:
        logger.info("med5")
        med5 = numpy.ma.median(chl_array,axis = 0)
        logger.info("median value: {}".format(med5))
        logger.info("median threshold: {}".format(median_threshold))
    
    std_dev = numpy.std(chl_array, axis=0)

    """
    median_output_name = median_filename

    if not os.path.exists(median_output_name):
        ods = nc.Dataset(median_output_name, "w")
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
    else:
        ods = nc.Dataset(median_output_name, "w")
        ods.variables['median'][:] = med5
        ods.close()
    """

    #! after estimating the median, i apply a filling to the chl time-series, e.g., window of 8 time steps, but it would be good to be able to choose other width of the window, e.g., 5 time-steps or 3 time-steps...

    #! here it would be good to give the option to select the value of the median threshold, e.g., median plus 5%, 10%, 15%...
    
    #get anomalies
    logger.info("anomaly")
    anomaly = chl_array - (med5*median_threshold)
    logger.info(anomaly[:5,:,debug_pixel[0],debug_pixel[1]])
    anomaly_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_anomaly"), mode="w+", shape=anomaly.shape, dtype=USE_DTYPE)
    anomaly_map[:] = anomaly[:]
    #anomaly = anomaly_map
    chl_array = None
    anomaly_map = None
    #need to ditch any empty entries here as they interefere with the cumsum

    #get cumsum of anomalies
    logger.info("chl cumsum")
    chl_cumsum = numpy.ma.cumsum(anomaly,axis=0)
    logger.info(chl_cumsum[:5,:,debug_pixel[0],debug_pixel[1]])
    chl_cumsum_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_cumsum"), mode="w+", shape=chl_cumsum.shape, dtype=USE_DTYPE)
    chl_cumsum_map[:] = chl_cumsum[:]
    #chl_cumsum = chl_cumsum_map
    anomaly = None
    chl_cumsum_map = None

    #get centered derivative        
    logger.info("chl der")
    chl_der = numpy.apply_along_axis(centered_diff_derivative, 0, chl_cumsum)
    logger.info(chl_der[:5,:,debug_pixel[0],debug_pixel[1]])
    chl_der_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_der"), mode="w+", shape=chl_der.shape, dtype=USE_DTYPE)
    chl_der_map[:] = chl_der[:]
    #chl_der = chl_der_map
    chl_cumsum = None
    chl_der_map = None
    

    #boxcar filter with width of 3 (sbx) should be something like this:
    logger.info("chl sbx")
    chl_boxcar = numpy.apply_along_axis(lambda m: numpy.convolve(m, numpy.ones(3)/3, mode='full'), axis=0, arr=chl_der)
    logger.info(chl_boxcar[:5,:,debug_pixel[0],debug_pixel[1]])
    fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])), mask=numpy.ones((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])))
    logger.info("chl shape after boxcar")
    logger.info(chl_boxcar.shape)
    chl_boxcar = numpy.concatenate([chl_boxcar, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr, fill_arr,])
    logger.info("chl shape after boxcar padding")
    logger.info(chl_boxcar.shape)
    chl_boxcar_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_sbx"), mode="w+", shape=chl_boxcar.shape, dtype=USE_DTYPE)
    chl_boxcar_map[:] = chl_boxcar[:]
    chl_boxcar = None
    chl_boxcar = chl_boxcar_map
    chl_boxcar_map = None

    return chl_boxcar.shape, USE_DTYPE, med5, std_dev

def create_phenology_netcdf(chl_lons, chl_lats, output_shape=None,name="phenology_{}.nc".format(datetime.datetime.now().strftime("%H%M")), date=False, median=None, std=None):
    """
    Creates the skeleton of the netcdf file to be used by write_to_output_netcdf, all of this is metadata.
    """
    if date:
        output_location_date = name
    else:
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
    ds.variables['probability'].setncattr("units", "likelihood")
    ds.createVariable('median_chlorophyll', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL)
    ds.createVariable('chlorophyll_std_dev', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL)
    ds.variables['median_chlorophyll'].setncattr("units", "mg chl m^3")
    ds.variables['chlorophyll_std_dev'].setncattr("units", "mg chl m^3")
    ds.setncattr("generation command", str(" ".join(sys.argv)))
    ds.setncattr("run location", str(os.getcwd()))
    if isinstance(median, numpy.ndarray):
        ds.variables['median_chlorophyll'][:] = median
    if isinstance(std, numpy.ndarray):
        ds.variables['chlorophyll_std_dev'][:] = std
    ds.close()
    logger.info("created netcdf {}".format(name))

def write_to_output_netcdf(data, total_blooms=None, probability=None, date=False, output_location_date=output_location_date, output_location=output_location):
    """
    Loops through each year in the numpy array and writes the data to the netcdf file, this should work faster if we get rid of the loop but I can't seem to grock the logic to fix it right now.
    """
    if date:
        ds = nc.Dataset(output_location_date,'r+',format='NETCDF4_CLASSIC')
        logger.info("writing to: {}".format(output_location_date))
    else:
        ds = nc.Dataset(output_location,'r+',format='NETCDF4_CLASSIC')
        logger.info("writing to: {}".format(output_location))
    data = data.astype(numpy.float32)
    data = numpy.ma.fix_invalid(data)
    logger.info("pre-writing data shape: {}".format(data.shape))
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
    logger.info("wrote {} years of data".format(len(ds.variables['TIME'][:])))
    ds.close()


def total_blooms_to_probability(array_like):
    return numpy.count_nonzero(array_like == 2) / array_like.size

def get_multi_year_two_blooms_output(numpy_storage, chunk, chl_shape, chl_dtype, chl_data, sst_shape, sst_dtype, date_seperation_per_year=47, start_date=0, reverse_search=20, reference_index=0, out_netcdf=output_location, out_date_netcdf=output_location_date):
    #this all works on the assumption the axis 0 is time
    logger.info("reading variables")
    chl_boxcar = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_sbx"), mode="r", dtype=USE_DTYPE, shape=chl_shape)
    chl_boxcar = chl_boxcar.copy()
    logger.info(sst_dtype)
    sst_der = numpy.memmap(os.path.join(numpy_storage, str(chunk), "sst_der"), mode="r", dtype=USE_DTYPE, shape=sst_shape)
    sst_der = sst_der.copy()
    logger.info("shapes after reading sst: {} chl: {}".format(chl_boxcar.shape, sst_der.shape))
    logger.info("reshaping to sst: {} chl: {}".format(sst_shape, chl_shape))
    sst_der.shape = sst_shape
    chl_boxcar.shape = chl_shape

    logger.info("doing yearly evaluation, this may take up a lot of memory, if so double check all memmaps have been flushed")
    #get start and ends, this works
    chl_boxcar = numpy.ma.masked_where((chl_boxcar == FILL_VAL), chl_boxcar)
    sst_der = numpy.ma.masked_where((sst_der == FILL_VAL), sst_der)
    #logger.info("doing chlorophyll initiations")
    #start_end_duration_array = numpy.apply_along_axis(get_start_index_and_duration, 0, year_chl_boxcar)
    logger.info((chl_data.shape[2],chl_data.shape[3], int((chl_data.shape[0] -start_date) // date_seperation_per_year), 2,5))
    year_true_start_end_array = numpy.ndarray((chl_data.shape[2],chl_data.shape[3], int((chl_data.shape[0] -start_date) // date_seperation_per_year), 2,5))
    logger.info(int((chl_data.shape[0] -start_date) // date_seperation_per_year))
    total_blooms = numpy.ndarray((chl_data.shape[2],chl_data.shape[3], int((chl_data.shape[0] -start_date) // date_seperation_per_year)))
    year_true_start_end_array.fill(FILL_VAL)
    blooms_by_date = year_true_start_end_array.copy()
    total_blooms.fill(FILL_VAL)
    total_blooms_date = total_blooms.copy()
    logger.info("doing sst initiations and correction")
    logger.debug("start date : {}".format(start_date))
    for ix, iy in tqdm.tqdm(numpy.ndindex(chl_data.shape[2], chl_data.shape[3]), total=(chl_data.shape[2] * chl_data.shape[3]), disable=True):
        try:
            verbose=False
            #if iy == debug_pixel[0] and ix == debug_pixel[1]:
            if iy == 11 and ix == 70:
                if logger.level <= logging.INFO:
                    logger.setLevel(logging.DEBUG)
                    logger.info("debug pixel {} {} encountered".format(ix,iy))
                    verbose = True
            else:
                logger.setLevel(default_logging)
                pass
            results = match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=verbose, start_date=start_date, reference_date=reference_index)
            year_true_start_end_array[ix,iy] = results[0]
            total_blooms[ix,iy] = results[1]
            blooms_by_date[ix, iy] = results[2]
            total_blooms_date[ix,iy] = results[3]

            logger.debug("end duration array")
            logger.debug(year_true_start_end_array[ix,iy])
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            logger.error(e)
            logger.error(repr(e))
            logger.error(results)
            logger.error((ix, iy))
            logger.error(exc_tb.tb_lineno)
            match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],chl_boxcar[:,:,ix,iy], chl_data[:,:,ix,iy], date_seperation_per_year, reverse_search, verbose=False, start_date=start_date, reference_date=reference_index)
        """
        if ix in completion_points and iy == 0:
            logger.info(completion_points.index(ix) * 10, "% complete")
        """
    
    logger.info(total_blooms.shape)
    probability_array = numpy.apply_along_axis(total_blooms_to_probability, 2, total_blooms)
    probability_array_date = numpy.apply_along_axis(total_blooms_to_probability, 2, total_blooms_date)

    logger.info("done sst initiations and correction")
    logger.info("writing to netcdf")
    #needs to be extended to be able to output 3 files: sorted by calendar year, one with primary and secondary chloropjhyll maximum
    write_to_output_netcdf(year_true_start_end_array, total_blooms=total_blooms, probability=probability_array, output_location=out_netcdf)
    write_to_output_netcdf(blooms_by_date, total_blooms=total_blooms_date, probability=probability_array_date, date=True, output_location_date=out_date_netcdf)

    
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
    output = numpy.ma.masked_invalid(output)
    return output, start_date + entries_per_year

def multi_proc_init(c, s, l):
    global sst_lock
    global chl_lock
    global log_lock
    sst_lock = s
    chl_lock = c
    log_lock = l

def test_process_for_exception(result):
    print("results")
    print("result ", result)

def chunk_by_chunk_handler(chunk_idx, chunk):
    global default_logging
    try:
        log = log_lock.acquire(block=False)
        if log:
            handler.setLevel(logging.DEBUG)
            logger.debug("set logging to debug on filehandler")
        logger.debug("processing chunk {} with bounds {}".format(chunk_idx, str(chunk)))
        chunk_by_chunk(chunk_idx, chunk)
        if log:
            log_lock.release()
            handler.setLevel(default_logging)

    except Exception as e:
        traceback.print_exc()
        log_lock.release()
        handler.setLevel(default_logging)
        logger.error(e)
    raise e

def chunk_by_chunk(chunk_idx, chunk):
    global chl_ds
    global args
    global start_date
    global debug_chunk
    global numpy_storage

    chunk_start_date = start_date

    if chunk_idx == debug_chunk:
        pass

    slc = [slice(None)] * len(chl_ds.variables[chl_variable].shape)
    x,y = chunk
    y = y if y[0] != 1 else (0, y[1])
    x = x if x[0] != 1 else (0, x[1])
    slc[LON_IDX] = slice(x[0], x[1])
    slc[LAT_IDX] = slice(y[0], y[1])

    output = chl_filename.replace(".nc", "_phenology_{}_chunk{}.nc".format(time_of_run, chunk_idx))

    chl_lock.acquire()
    chl_array = chl_ds.variables[chl_variable][slc]
    chl_lons = chl_ds.variables[chl_lon_var][x[0]:x[1]]
    chl_lats = chl_ds.variables[chl_lat_var][y[0]:y[1]]
    #mask if it isn't, else will raise errors
    if not numpy.ma.is_masked(chl_array):
        chl_array = numpy.ma.masked_array(chl_array, numpy.isnan(chl_array))
    if chl_array.mask.all():
        logger.info(numpy.isnan(chl_array).all())
        logger.info("skipping as empty")
        chl_lock.release()
        #output empty netcdf
        empty_med = numpy.empty((chl_lats.shape[0], chl_lons.shape[0]), dtype=chl_array.dtype) if LAT_IDX >LON_IDX else numpy.empty(shape=(chl_lons.shape[0], chl_lats.shape[0]), dtype=chl_array.dtype)
        create_phenology_netcdf(chl_lons, chl_lats, [1,1], output.replace(".nc", "_by_maxval.nc"), median=empty_med, std=empty_med)
        create_phenology_netcdf(chl_lons, chl_lats, [1,1], output.replace(".nc", "_by_date.nc"), date=True, median=empty_med, std=empty_med)
        return True
    chl_lock.release()
    if args.reshape:
        chl_array.shape = (chl_array.shape[0], 1, chl_array.shape[1], chl_array.shape[2])

    #check if there are any nans in the data
    chl_array = numpy.ma.masked_invalid(chl_array)
    logger.info("making temp storage")
    os.mkdir(os.path.join(numpy_storage,str(chunk_idx)))


    if args.sst_location:
        sst_lock.acquire()
        logger.info("sst file provided, reading array")
        logger.info("only one file found, assuming full stack of observations")
        sst_ds = nc.Dataset(args.sst_location)
        sst_variable = [x for x in sst_ds.variables if args.sst_var in x.lower()][0]
        try:
            sst_lon_var = [x for x in sst_ds.variables if "lon" in x.lower()][0]
            sst_lat_var = [x for x in sst_ds.variables if "lat" in x.lower()][0]
            sst_time_var = [x for x in sst_ds.variables if "time" in x.lower()][0]
        except:
            logger.info("trying to get lat/lon from standard name")
            for va in sst_ds.variables:
                if "standard_name" in sst_ds.variables[va].__dict__.keys():
                    if sst_ds.variables[va].__dict__["standard_name"] == "latitude":
                        sst_lat_var= va
                    elif sst_ds.variables[va].__dict__["standard_name"] == "longitude":
                        sst_lon_var = va
                    elif sst_ds.variables[va].__dict__["standard_name"] == "time":
                        sst_time_var = va
        sst_lons = sst_ds.variables[sst_lon_var][:]
        sst_lats = sst_ds.variables[sst_lat_var][:]
        SST_LAT_IDX = sst_ds.variables[sst_variable].dimensions.index(sst_lat_var)
        SST_LON_IDX = sst_ds.variables[sst_variable].dimensions.index(sst_lon_var)
        SST_TIME_IDX = sst_ds.variables[sst_variable].dimensions.index(sst_time_var)
        sst_slc = [slice(None)] * len(sst_ds.variables[sst_variable].shape)
        sst_slc[SST_LON_IDX] = slice(x[0], x[1])
        sst_slc[SST_LAT_IDX] = slice(y[0], y[1])
        sst_slc[SST_TIME_IDX] = slice(args.sst_start_index, args.sst_end_index)
        sst_array = sst_ds.variables[sst_variable][sst_slc]
        if not numpy.ma.is_masked(sst_array):
            sst_array = numpy.ma.masked_array(sst_array, numpy.isnan(sst_array))

        sst_lock.release()

        if args.extend_sst_data:
            sst_array, _ = extend_array(sst_array, date_seperation_per_year, chunk_start_date)

        if args.reshape_sst:
            sst_array.shape = (sst_array.shape[0], 1, sst_array.shape[1], sst_array.shape[2])

        sst_shape, sst_dtype = prepare_sst_variables(sst_array, numpy_storage, chunk_idx, skip=args.skip_sst_prep)
        logger.info("sst_shape: {}".format(sst_shape))
        logger.info("sst_dtype: {}".format(USE_DTYPE))
        sst_array = None

    """
    if not (chl_array.shape[2] == chl_lats.shape[0] and chl_array.shape[3] == chl_lons.shape[0]):
        logger.info("adjusting to flip lat and lon")
        chl_array.shape = (chl_array.shape[0], chl_array.shape[1], chl_array.shape[3], chl_array.shape[2])
    """
    if args.extend_chl_data:
        chl_array, chunk_start_date = extend_array(chl_array, date_seperation_per_year, chunk_start_date)
        logger.info("start date after extension: {}".format(chunk_start_date))

    chl_shape, chl_dtype, chl_median, chl_std_dev = prepare_chl_variables(chl_array, numpy_storage, chunk_idx, date_seperation_per_year, chl_lats, chl_lons, do_model_med=args.modelled_median, median_threshold=med_thresh, median_filename=chl_filename.replace(".nc", "_median_{}_chunk{}.nc".format(time_of_run, chunk_idx)))
    logger.info("chl_shape: {}".format(chl_shape))
    logger.info("chl_dtype: {}".format(USE_DTYPE))

    if sst_shape[2:] != chl_shape[2:]:
        logger.error("sst and chlorophyll x,y array shapes do not match got:")
        logger.error("chlorophyll: {}".format(chl_shape[2:]))
        logger.error("sst: {}".format(sst_shape[2:]))
        logger.error("quitting!")
        sys.exit()
    
    #simple regridding
    create_phenology_netcdf(chl_lons, chl_lats, chl_shape, output.replace(".nc", "_by_maxval.nc"), median=chl_median, std=chl_std_dev)
    create_phenology_netcdf(chl_lons, chl_lats, chl_shape, output.replace(".nc", "_by_date.nc"), date=True, median=chl_median, std=chl_std_dev)
    logger.info("using start date {}".format(chunk_start_date))
    get_multi_year_two_blooms_output(numpy_storage, 
                                    chunk_idx,
                                    chl_shape,
                                    chl_dtype, 
                                    chl_array, 
                                    sst_shape, 
                                    sst_dtype, 
                                    date_seperation_per_year=date_seperation_per_year, 
                                    start_date=chunk_start_date, 
                                    reverse_search=reverse_search,
                                    reference_index=ref_index,
                                    out_date_netcdf=output.replace(".nc", "_by_date.nc"),
                                    out_netcdf=output.replace(".nc", "_by_maxval.nc"))
    
    logger.setLevel(default_logging)
    shutil.rmtree(os.path.join(numpy_storage,str(chunk_idx)))
    return True

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
    parser.add_argument("--median_threshold", default=MEDIAN_THRESHOLD_DEFAULT, help="change median threshold, set as percentage value e.g. 20 = 20", required=False)
    parser.add_argument("--modelled_median", action='store_true',default=False, help="test for solar zenith of inputs", required=False)
    #probably needs a better description!
    #give options for both since we might end up in a situation where sst is 3 years and chl is one year (or vice versa)
    parser.add_argument("--extend_chl_data", default=False, action="store_true", help="extends chlorophyll by copying the (central) chl array to the year previous and year next")
    parser.add_argument("--extend_sst_data", default=False, action="store_true", help="extends sea surfaace temperature by copying the (central) chl array to the year previous and year next")
    parser.add_argument("--start_year", default=0, help="What year to use as the start point for metadata output, if not specified will use 0, affects no processing.")
    parser.add_argument("--debug_pixel",  nargs='+', default=None, type=int, help="pixle in x, y (lat, lon), these entries are 0 indexed.")
    parser.add_argument("--reshape", default=False, action="store_true", help="reshape to be t, 1, x, y")
    parser.add_argument("--reshape_sst", default=False, action="store_true", help="reshape to be t, 1, x, y")
    parser.add_argument("--sst_start_index", type=int, default=0)
    parser.add_argument("--sst_end_index", type=int, default=-1)
    parser.add_argument("--chl_begin_date", default=None, help="date relating to index 0 of a file in format dd/mm/yyyy")
    parser.add_argument("--no_thread", action='store_true', default=False, help="Don't use threading")
    parser.add_argument("--no_logfile", action='store_true', default=False, help="Don't create a logfile")
    args = parser.parse_args()

    time_of_run = datetime.datetime.now().strftime("%H%M")
    if not args.no_logfile:
        if len(args.chl_location) == 1:
            handler = logging.FileHandler(args.chl_location[0].replace(".nc", "_phenology_{}.logfile".format(time_of_run)))
            logger.info("initilised logfile at {}".format(args.chl_location[0].replace(".nc", "_phenology_{}.logfile".format(time_of_run))))
        else:
            handler = logging.FileHandler("phenology_bulk_run_{}.logfile".format(time_of_run))
            logger.info("initilised logfile at {}".format("phenology_bulk_run_{}.logfile".format(time_of_run)))
    else:
        handler = logging.NullHandler()

    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.info(" ".join(sys.argv))

    med_thresh = 1+ (float(args.median_threshold) / 100)
    if not args.intermediate_file_store:
        numpy_storage = tempfile.mkdtemp(prefix="phenology_")
    else:
        numpy_storage = args.intermediate_file_store
        if os.path.exists(numpy_storage):
            shutil.rmtree(numpy_storage)
        os.mkdir(numpy_storage)
    #remember to change median threshold to percentage!
    #reverse search should be ratio of 100 days so 5 days * 20 = 100 or 8 days * 12.5 (so 13) = 100 days
    #as 100 days is 0.27 of 365 we can just do that
    date_seperation_per_year = int(args.date_seperation_per_year)
    if not args.reverse_search:
        reverse_search = int(round(int(args.date_seperation_per_year) * 0.28))
    logger.info("Reverse search:{}".format(reverse_search))

    start_date = None
    start_date = int(args.first_date_index) - 1 if not start_date else start_date
    ref_index = int(args.reference_index) - 1

    if (ref_index + start_date) > (date_seperation_per_year + start_date):
        logger.info("reference index is too great and would result in mostly or all negative values.")
        sys.exit()

    #TODO list of files or file specified (mid november)
    for chl_location in args.chl_location:
        chl_files = chl_location
        logger.info("calculating on {}".format(chl_files))
        chl_filename = chl_location
        chl_ds = nc.Dataset(chl_location)
        #test for more than 1 variable, if so quit out and complain that it doesn't know which to use
        chl_variable = [x for x in chl_ds.variables if args.chl_var in x.lower()][0]
        chl_lon_var = [x for x in chl_ds.variables if "lon" in x.lower()][0]
        chl_lat_var = [x for x in chl_ds.variables if "lat" in x.lower()][0]
        chl_time_var = [x for x in chl_ds.variables if "time" in x.lower()][0]
        try:
            if not args.chl_begin_date:
                dts = nc.num2date(chl_ds[chl_time_var][:], chl_ds[chl_time_var].units)
                ref_date = dts[start_date+ref_index]
                start_datetime = dts[start_date]
            else:
                init_date = datetime.datetime.strptime(args.chl_begin_date, '%d/%m/%Y')
                start_datetime = init_date + datetime.timedelta(days=(365 * (start_date / date_seperation_per_year)))
                ref_date = start_datetime + datetime.timedelta(days=(365 * (ref_index / date_seperation_per_year)))
            logger.info(ref_date)
            REF_MONTH = calendar.month_name[ref_date.month]
            START_YEAR = start_datetime.year
            logger.info("reference date selected: "+ref_date.strftime("%d/%m/%Y"))
            logger.info("start date selected: "+start_datetime.strftime("%d/%m/%Y"))
            logger.info("reference month: "+ REF_MONTH)
            logger.info("start year: "+str(START_YEAR))
        except Exception as e:
            logger.error(e)
            logger.error("could not identify chlorophyll start date, please specify the first date you expect in the file (time index 0, not first date index) with --chl_begin_date e.g. --chl_begin_date 01/01/1997")
            sys.exit()
        chl_lons = chl_ds.variables[chl_lon_var][:]
        chl_lats = chl_ds.variables[chl_lat_var][:]
        logger.info("lats shape {}".format(chl_lats.shape))
        logger.info("lons shape {}".format(chl_lons.shape))
        LAT_IDX = chl_ds.variables[chl_variable].dimensions.index(chl_lat_var)
        LON_IDX = chl_ds.variables[chl_variable].dimensions.index(chl_lon_var)
        if LAT_IDX > LON_IDX:
            DIM_ORDER = ['TIME', 'DEPTH', 'LONGITUDE', 'LATITUDE']
        else:
            DIM_ORDER = ['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE']
        logger.info("done with dimensions")
        
        logger.info("getting chunks")

        """
            debug_pixel = args.debug_pixel if args.debug_pixel else [int(chl_lats.shape[0] * 0.45), int(chl_lons.shape[0] * 0.45)]
            logger.info("using debug pixel:", debug_pixel, " which equates to:", chl_lats[debug_pixel[1]], "N", chl_lons[debug_pixel[0]], "E", "(zero indexed so you may need to add 1 to get reference in other software)")
            debug_pixel = debug_pixel if LON_IDX > LAT_IDX else [debug_pixel[1], debug_pixel[0]]
        """
        debug_pixel = (50,50)


        chunks = [(x, y) for x in zip(list(range(0, chl_lons.shape[0], 100)),list(range(100,chl_lons.shape[0], 100))+ [chl_lons.shape[0]]) 
                         for y in zip(list(range(0, chl_lats.shape[0], 100)), list(range(100, chl_lats.shape[0], 100)) + [chl_lats.shape[0]])]


        logger.info("done chunks")
        logger.info("processing will take {} chunks, at 3 minutes a chunk this will take approximately {} minutes".format(len(chunks), len(chunks) * 3))

        sizes = [item for chunk in chunks for item in [chunk[0][1] -chunk[0][0], chunk[1][1] -chunk[1][0]]]
        min_size = min(sizes)
        debug_pixel = (min_size // 2, min_size // 2)
    
        default_logging = logging.INFO

        if not args.no_thread:
            threads = multiprocessing.cpu_count()
        else:
            threads = 1
        if threads != 1:
            logger.info("threading enabled, revised estimate is {} minutes using the current number of threads ({})".format((len(chunks) *3) /threads, threads))
            logger.info("setting log level to error only")
            logger.setLevel(logging.ERROR)
            default_logging = logging.ERROR
        sst_lock_maj = multiprocessing.Lock()
        chl_lock_maj = multiprocessing.Lock()
        log_info_lock_maj = multiprocessing.Lock()
        debug_chunk = 1
        pool = multiprocessing.Pool(threads, initializer=multi_proc_init, initargs=(chl_lock_maj, sst_lock_maj, log_info_lock_maj))
        #set this to be an early one
        res = pool.starmap_async(chunk_by_chunk_handler, [(chunk_idx, chunk) for chunk_idx, chunk in enumerate(chunks)], chunksize=1)
        with tqdm.tqdm(total=len(chunks)) as pbar:
            running = True
            last_known = res._number_left
            while running:
                time.sleep(1)
                if res._number_left != last_known:
                    pbar.update(last_known - res._number_left)
                    last_known = res._number_left
                if res._number_left <= 0 or last_known <= 0:
                    pbar.close()
                    running = False
                pbar.refresh()
        logger.setLevel(logging.INFO)
        logger.info("chunk processing finished, returned log level to info messages")

        final_output_maxval = chl_filename.replace(".nc", "_phenology_{}_by_maxval.nc".format(time_of_run))
        final_output_dates = chl_filename.replace(".nc", "_phenology_{}_by_date.nc".format(time_of_run))
        logger.info("saving reconstructed files to {}, {}".format(final_output_dates, final_output_maxval))

        chl_lon_var = [x for x in chl_ds.variables if "lon" in x.lower()][0]
        chl_lat_var = [x for x in chl_ds.variables if "lat" in x.lower()][0]
        chl_lons = chl_ds.variables[chl_lon_var][:]
        chl_lats = chl_ds.variables[chl_lat_var][:]

        create_phenology_netcdf(chl_lons, chl_lats, [1,1], final_output_maxval)
        create_phenology_netcdf(chl_lons, chl_lats, [1,1], final_output_dates)
        files = glob.glob(chl_filename.replace(".nc","_phenology_{}_chunk*_*").format(time_of_run))
        for chunk_idx, chunk in enumerate(chunks):
            maxval_chunk_file = [x for x in files if "chunk{}_".format(chunk_idx) in x and 'maxval' in x]
            date_chunk_file = [x for x in files if "chunk{}_".format(chunk_idx) in x and 'date' in x]
            #if we find the tile then add it in
            if len(maxval_chunk_file) and len(date_chunk_file):
                logger.info("processing chunk files: {}, {}".format(os.path.basename(maxval_chunk_file[0]), os.path.basename(date_chunk_file[0])))
                maxval_ds = nc.Dataset(final_output_maxval, 'r+')
                date_ds = nc.Dataset(final_output_dates, 'r+')
                maxval_chunk = nc.Dataset(maxval_chunk_file[0], 'r')
                date_chunk = nc.Dataset(date_chunk_file[0], 'r')
                logger.info("Files read, preparing to write")
                slc = [slice(None)] * 4
                med_prob_slc = [slice(None)] * 3
                x,y = chunk
                y = y if y[0] != 1 else (0, y[1])
                x = x if x[0] != 1 else (0, x[1])
                if LAT_IDX > LON_IDX:
                    slc[2] = slice(x[0], x[1])
                    slc[3] = slice(y[0], y[1])
                    med_prob_slc[1] = slice(x[0], x[1])
                    med_prob_slc[2] = slice(y[0], y[1])
                else:
                    slc[3] = slice(x[0], x[1])
                    slc[2] = slice(y[0], y[1])
                    med_prob_slc[2] = slice(x[0], x[1])
                    med_prob_slc[1] = slice(y[0], y[1])
                for var in [x for x in maxval_ds.variables.keys() if not x in maxval_ds.dimensions.keys() and x not in ['median_chlorophyll', 'probability']]:
                    try:
                        maxval_ds[var][slc] = maxval_chunk[var][:]
                        date_ds[var][slc] = date_chunk[var][:]
                    except ValueError:
                        logger.warning("Encountered potentially empty file during stitching, this might not be a problem but if you see a large square area of empty data this is likely the cause. Found in chunk {}".format(chunk_idx))
                for var in [x for x in maxval_ds.variables.keys() if x in ['median_chlorophyll', 'probability']]:
                    maxval_ds[var][med_prob_slc] = maxval_chunk[var][:]
                    date_ds[var][med_prob_slc] = date_chunk[var][:]
                logger.info("Written, closing.")
                maxval_ds.close()
                date_ds.close()
                maxval_chunk.close()
                date_chunk.close()
                logger.info("Removing chunk files.")
                os.remove(maxval_chunk_file[0])
                os.remove(date_chunk_file[0])

        logger.info("complete!")


                

