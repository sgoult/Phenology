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
import astropy.convolution
from fast_polarity.polarity import polarity_edge_finder_optimised as polaritiser
from scipy.signal import argrelextrema, argrelmin

#TODO Set dynamically from the input netcdf or user specified from command line
FILL_VAL = -9.999999999999998e+33

#TODO set this from the user inputs
MEDIAN_THRESHOLD_DEFAULT = 20
STD_DEV_THRESHOLD_DEFAULT = 2
MAX_MEANS_THRESHOLD_DEFAULT = 20
LAT_IDX = 2
LON_IDX = 3
REF_MONTH = 'January'
END_MONTH = 'December'
USE_DTYPE = "float64"
DEFAULT_CHUNKS = 200

output_location = None
output_location_date = None
median_output_name = None

logging.basicConfig()
logger = logging.getLogger("phenology_logger")
logger.setLevel(logging.DEBUG)


def find_maxes(array_like):
    max_indexes = argrelextrema(array_like, numpy.greater)
    min_indexes = argrelextrema(array_like, numpy.less)
    maxes = array_like[max_indexes]
    minimums = array_like[min_indexes]
    min_avg = numpy.ma.median(minimums)
    max_avg = numpy.ma.median(maxes)
    return max_avg, min_avg, (max_avg - min_avg) / max_avg


def reversedEnumerate(l):
    return zip(range(len(l)-1, -1, -1), l)

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
    cf = numpy.ma.convolve(squeezed, [1,0,-1],'same') / 1
    return cf

def get_start_index_and_duration(array_like,chl_values,date_offset,depth=5, pad_values=False, verbose=False):
    """
    takes a list of values and searches for the max, then the associated sign changes to indicate bloom start/end
    set depth to change limit of number of identifications
    
!comment[A] my approach was to save the  max Chl value found between each of the start and end time, and also the duration, i.e., estimated as number of steps between start and end times
! with the first derivative light availability or SST is increasing when PAR or SST first derivative is positive, and vice versa
 
    in a run using global data that took 30 minutes, this function made up 513 seconds of the processing time
    """
    global duration_minimum_value
    global start_minimum_seperation_value
    #before we do anything else, because the boxcar method gives us an output that is 0.0 in place of nans we should increment those by 0.0001 so that the polariser can find them
    array_like[array_like==0.0] = 0.0001
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
    #find out which way we're going
    starts = []
    ends = []
    #works much the same qas the SST one below
    #TODO work out why we have dropped the first index (1) from starts - theres no reason it should be doing this!!
    poss_starts = [x for x in zero_crossings if array_like[x] < 0]
    poss_ends = [x for x in zero_crossings if array_like[x] > 0]
    logger.debug("starts: {}".format(poss_starts))
    logger.debug("ends: {}".format(poss_ends))
    starts = poss_starts
    #flip this around!
    already_tested = False
    ends = []
    for start in starts:
        try:
            ends.append(next(x for x in poss_ends if x > start))
        except:
            continue
    #remove me
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
        logger.debug((idx, start))
        try:
            if durations[idx] <= duration_minimum_value:
                continue
        except IndexError:
            continue
        try:
           end = next(x for x in ends if x > start)
           max_idx = numpy.nanargmax(chl_values[start:end + 1])
           #we find the "max" value from the chl boxcar
           chl_val = chl_values[max_idx + date_offset + start] if not chl_values.mask[max_idx + date_offset + start] else numpy.nan
           if numpy.isnan(chl_val):
               continue
           dates.append([start + date_offset,end + date_offset,end-start,max_idx + date_offset + start - 1,chl_val])
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


def phen_records_to_one_val_on_max(records, year_start_index=0, index=4, verbose=False):
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
        if date_seperation_per_year:
            output_record[0] = output_record[0] - year_start_index
            output_record[1] = output_record[1] - year_start_index
            output_record[3] = output_record[3] - year_start_index
        return output_record
    else:
        return [None,None,None,None,None]

def polarity_edge_finder(input_array):
    return (numpy.diff(numpy.sign(input_array)) != 0)*1

def match_start_end_to_solar_cycle(array_like, chl_sbx_slice, chl_slice, date_seperation_per_year, reverse_search, start_date=0, reference_date=0, verbose=False, skip_end_year=False, processing_end_date=None, one_year_climato=False):
    """
    Attributes the start and end times in relation to the SST or solar cycle, takes an sst array (array_like), smoothed chlorophyll derivative slice (chl_sbx_slice) and the original chlorophyll data.

    Slices up the data based on high/low periods of SST (or otherwise), then feeds each period into get_start_index_and_duration, once finished it will output an array of shape (x, y, time, 2, 5)
    verbose will spame the terminal with information about whats going on, best to establish a few pixels you want to inspect rather than having this on all the time.

    in a run using global data that too 30 minutes, this function made up 703 seconds of the processing time
    I would guess that 500 of those seconds can be attributed to get_start_index_and_duration
    """ 
    if not processing_end_date:
        processing_end_date = chl_sbx_slice.shape[0]
        if skip_end_year:
            processing_end_date = processing_end_date - date_seperation_per_year

    logger.debug(f"using start date of {start_date} and end date of {processing_end_date}")

    if array_like.mask.all():
        logger.info("array all -32616.30273438")
        logger.info(array_like)
        #we can stop here, there's no point continuing with an empty array
        return [[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, processing_end_date, date_seperation_per_year) if not x + date_seperation_per_year > processing_end_date],[[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, processing_end_date, date_seperation_per_year) if not x + date_seperation_per_year > processing_end_date]

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
            return [[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, processing_end_date, date_seperation_per_year) if not x + date_seperation_per_year > processing_end_date],[[None,None,None,None,None], [None,None,None,None,None]], [None for x in range(start_date, processing_end_date, date_seperation_per_year) if not x + date_seperation_per_year > processing_end_date]
    except Exception as e:
        #triggered once, but was helpful to know what the contents were
        logger.error(highs)
        logger.error(lows)
        logger.error(e)

    #chuck out the results
    logger.debug("updated highs and lows after filtering and correcting for missing dates")
    logger.debug("highs: {}".format(numpy.array2string(numpy.array(highs))))
    logger.debug("lows: {}".format(numpy.array2string(numpy.array(lows))))
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

        logger.debug("blooms in this period")
        if activity_period:
            logger.debug([x for x in period_chl_phenology if x[3] >= index and x[3] < end_date and not all(p is None for p in x)])
            high_records.extend([x for x in period_chl_phenology if x[3] >= index and x[3] < end_date and not all(p is None for p in x)])
        else:
            logger.debug([x for x in period_chl_phenology if x[3] >= index and x[3] < end_date and not all(p is None for p in x)])
            low_records.extend([x for x in period_chl_phenology if x[3] >= index and x[3] < end_date and not all(p is None for p in x)])
        
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
    #remove any year that doesn't have full 12 months of coverage
    logger.debug("list of years being considered:")
    logger.debug(list(enumerate(range(start_date, processing_end_date, date_seperation_per_year))))
    for year, year_start_index in enumerate(range(start_date, processing_end_date, date_seperation_per_year)):
        logger.debug("doing year {}".format(year))
        logger.debug("this year starts on {}".format(year_start_index))
        logger.debug("this year ends on {}".format(year_start_index+date_seperation_per_year))
        logger.debug(f"our end date is {processing_end_date}")
        if ((year_start_index + date_seperation_per_year) > processing_end_date) or (year_start_index >= processing_end_date):
            logger.debug("skipping year as start or end are greater than end date")
            break


        sst_highs= [high for high in highs if high > year_start_index and high < (year_start_index + date_seperation_per_year)]
        sst_lows = [low for low in lows if low > year_start_index and low < (year_start_index + date_seperation_per_year)]

        try:
            logger.debug("doing high start end detection")
            last_high = max(sst_highs)
            logger.debug(last_high)
            first_high = min(sst_highs)
            logger.debug(first_high)
            last_high_end = next(x for x in lows if x > last_high)
            logger.debug(last_high_end)
            pre_high_start = next(x for x in reversed(lows) if x < first_high) 
            logger.debug(pre_high_start)
        except Exception as e:
            logger.debug(e)
            first_high = -1000
            last_high = -1000
            last_high_end = -1000
            pre_high_start = -1000
        
        try:
            logger.debug("doing low start end detection")
            first_low = min(sst_lows)
            logger.debug(first_low)
            last_low = max(sst_lows)
            logger.debug(last_low)
            logger.debug(highs)
            last_low_end = next(x for x in highs if x > last_low and x != last_low)
            logger.debug(last_low_end)
            pre_low_start = next(x for x in reversed(highs) if x < first_low)
            logger.debug(pre_low_start)
        except Exception as e:
            logger.debug(e)
            first_low = -1000
            last_low = -1000
            pre_low_start= -1000
            last_low_end = -1000


        if last_high_end > last_low_end:
            logger.debug("selected high end")
            forward_search = (last_high_end - (year_start_index + date_seperation_per_year)) + 1
        else:
            logger.debug("selected low end")
            forward_search = (last_low_end - (year_start_index + date_seperation_per_year)) + 1

        if forward_search == -1000:
            forward_search = 0

        logger.debug("forward search is: {}".format(forward_search))

        logger.debug("lows for this year: {}".format(numpy.array2string(numpy.array(sst_lows))))
        logger.debug("highs for this year: {}".format(numpy.array2string(numpy.array(sst_highs))))

        logger.debug("first low start: {}".format(first_low))
        logger.debug("high start before first low start: {}".format(pre_low_start))
        logger.debug("first high start: {}".format(first_high))
        logger.debug("low start before first high start: {}".format(pre_high_start))
        logger.debug("last high end: {}".format(last_high_end))
        logger.debug("last low end: {}".format(last_low_end))
        #find blooms that start after the year - our reverse search, end before the end of the year, and end during the current year
        #this is where I have concerns that there are problems with the selection of blooms, if we're doing a year with reference month of June
        #then this currently selects june - june, rather than january - december, is this correct? The central month 
        possible_high_blooms = [x for x in high_records if x[3] > year_start_index and x[3] < (year_start_index + date_seperation_per_year + forward_search) and not x[0] < (year_start_index - reverse_search)]
        possible_low_blooms = [x for x in low_records if x[3] > year_start_index and x[3] < (year_start_index + date_seperation_per_year + forward_search) and not x[0] < (year_start_index - reverse_search)]

        logger.debug("possible_high_blooms pre filtering")
        logger.debug(possible_high_blooms)
        logger.debug("possible_low_blooms pre filtering")
        logger.debug(possible_low_blooms)


        #if there is a high that straddles the years
        if pre_low_start:
            if pre_low_start < year_start_index:
                if abs(pre_low_start - year_start_index) > abs(first_low - year_start_index):
                    # if the high "belongs" to the previous year (ie the majority of it resides in the last year) then we should discount any highs that occur in there
                    possible_high_blooms = [x for x in possible_high_blooms if x[3] > first_low]

        #if there is a low that straddles the years
        if pre_high_start:
            if pre_high_start < year_start_index:
                if abs(pre_high_start - year_start_index) > abs(first_high - year_start_index):
                    # if the low "belongs" to the previous year (ie the majority of it resides in the last year) then we should discount any highs that occur in there
                    possible_low_blooms = [x for x in possible_low_blooms if x[3] > first_high]


        #test the start date, and the end date, if the start and end date are roughly equal
        #select the next max bloom

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
        #CORRECT VALUES
        low = phen_records_to_one_val_on_max(low_blooms, year_start_index=year_start_index)
        high = phen_records_to_one_val_on_max(high_blooms, year_start_index=year_start_index)
        if not low or not high:
            logger.debug("***************")
            logger.debug(high_records)
            logger.debug(low_records)
            logger.debug(year_start_index - reverse_search)
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

        
        #Alternative date sorting, based on start date
        high_start = high[0] if high[0] else -1000
        low_start = low[0] if low[0] else -1000
        if low_start > high_start:
            blooms_by_date.append([low, high])
        else:
            blooms_by_date.append([high, low])
        

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

    return blooms, ngds, blooms_by_date, ngds_date, total_blooms

def prepare_sst_variables(sst_array, chunk, skip=False, chunk_idx=None, output_name=None):
    """
    Creates smoothed sst, currently has a large portion commented out as source file
     is already the centered derivative diff data.
    """
    global debug_pixel_main
    ds = nc.Dataset(output_name.replace(".nc", "_intermediate_products.nc".format(chunk_idx)), 'r+', format='NETCDF4_CLASSIC')
    logger.info("sst_array shape before prep:{}".format(sst_array.shape))
    #smoothed sst
    if not skip:
        logger.info("smoothing sst")
        #logger.info(sst_array[:20,:,debug_pixel_main[0],debug_pixel_main[1]])
        if debug_pixel_main[2] == chunk_idx:
            logger.info(sst_array[:5,:,debug_pixel_main[0],debug_pixel_main[1]])
        sst_boxcar = numpy.apply_along_axis(lambda m: numpy.ma.convolve(m, numpy.ones(8)/8, mode='valid'), axis=0, arr=sst_array)
        logger.info("shape after sst sbx")
        logger.info(sst_boxcar.shape)
        fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])), mask=numpy.ones((1,1,sst_boxcar.shape[2],sst_boxcar.shape[3])))
        #sst_boxcar_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "sst_sbx"), mode="w+", shape=sst_boxcar.shape, dtype=USE_DTYPE)
        logger.debug("masking sst boxcar with sst_array mask")

        ds.createVariable('sst_boxcar', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)
        ds.variables['sst_boxcar'].setncattr("units", "degrees celsius")
        ds.variables['sst_boxcar'][:] = sst_boxcar[:]

        ds.close()
        #get sst derivative
        logger.info("doing sst_derivative")
        sst_der = numpy.apply_along_axis(centered_diff_derivative, 0, sst_boxcar[:,:,:,:])
        logger.info("shape after sst der")
        logger.info(sst_der.shape)
    else:
        logger.info("Skipped sst data preparation")
        sst_der = sst_array
        numpy.ma.set_fill_value(sst_array, fill_value=FILL_VAL)
        logger.info(sst_array.shape)

    new_timesteps = sst_der.shape[0] - sst_array.shape[0]
    logger.info("after sst preparation timestep difference is {}, original shape: {}, new shape: {}".format(new_timesteps, sst_array.shape[0], sst_der.shape[0]))
    additional_steps = new_timesteps / 2
    logger.info("this means there are {} new timesteps that must be created at the beginning and end of the sst derivitive array".format(abs(additional_steps)))

    
    start_fill_arrays = [fill_arr for i in range(0, missing_sst_dates_at_start + math.floor(abs(additional_steps)))]

    if not additional_steps.is_integer():
        logger.warning("as the number of new timesteps is not a whole number will pad one additional step to the beginning of the sst derivitive array")

    end_fill_arrays = [fill_arr for i in range(0, missing_sst_dates_at_end + math.ceil(abs(additional_steps)))]

    logger.info((sst_der.shape[0], len(start_fill_arrays), len(end_fill_arrays)))
    if not (sst_der.shape[0] + len(start_fill_arrays) + len(end_fill_arrays)) % date_seperation_per_year == 0:
        logger.error("After padding, the sst derivitive end product did not divide equally! shape is {} and division is {}".format((sst_der.shape[0] + len(start_fill_arrays) + len(end_fill_arrays)), (sst_der.shape[0] + len(start_fill_arrays) + len(end_fill_arrays)) % date_seperation_per_year))
        raise Exception("After padding, the sst derivitive end product did not divide equally!")
    

    sst_der = numpy.ma.concatenate(start_fill_arrays + [sst_der] + end_fill_arrays)
    sst_der = numpy.ma.masked_where(sst_der == 0, sst_der)
    sst_der = numpy.ma.filled(sst_der, fill_value=FILL_VAL)
    ds = nc.Dataset(output_name.replace(".nc", "_intermediate_products.nc".format(chunk_idx)), 'r+', format='NETCDF4_CLASSIC')
    ds.createVariable('sst_der', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)
    ds.variables['sst_der'].setncattr("units", "deg C")
    ds.variables['sst_der'][:] = sst_der[:]
    ds.close()


    logger.info("sst prep complete")
    return sst_der.shape, 'float64'

def prepare_chl_variables(chl_array, chunk, date_seperation, chl_lats, chl_lons, modelled_threshold=False, median_threshold=1.2, relative_max_anomaly=False, relative_median_anomaly=False, dynamic_max_means_threshold=0.25, max_means_threshold=0.2, median_filename="median_output.nc", do_fill_chl=True, std_dev_anomaly=False,  max_means_anomaly=False, std_dev_threshold=STD_DEV_THRESHOLD_DEFAULT, chunk_idx=None, output_name=None, date_seperation_per_year=0):
    """
    Creates the smoothed anomaly chlorophyll data, saves a file to the temporary directory that is read as a mem map later to conserve resources.
    """
    #median * 1.05
    global debug_pixel_main
    #! here i would prfer we output the media value (without adding 5% to it)

    ds = nc.Dataset(output_name.replace(".nc", "_intermediate_products.nc".format(chunk_idx)),'w',format='NETCDF4_CLASSIC')
    ds.createDimension('LONGITUDE', chl_lons.shape[0])
    ds.createDimension('LATITUDE', chl_lats.shape[0])
    ds.createDimension('DEPTH', 1)
    ds.createDimension('TIME', None)
    ds.createVariable('LATITUDE', 'float64', dimensions=['LATITUDE'], zlib=True)
    ds.variables['LATITUDE'].setncattr("units", "degrees_north")
    ds.variables['LATITUDE'][:] = chl_lats
    ds.createVariable('LONGITUDE', 'float64', dimensions=['LONGITUDE'], zlib=True)
    ds.variables['LONGITUDE'].setncattr("units", "degrees_east")
    ds.variables['LONGITUDE'][:] = chl_lons
    ds.createVariable('DEPTH', 'float32', dimensions=['DEPTH'], zlib=True)
    ds.variables['DEPTH'].setncattr("units", "meters")
    ds.variables['DEPTH'].setncattr("positive", "down")
    ds.variables['DEPTH'][:] = [0.1]
    ds.createVariable('TIME', 'float32', dimensions=['TIME'], zlib=True)
    ds.variables['TIME'].setncattr("units", "years")
    #switch back to LATITIUDE and LONGITUDE, establish why the flipping of the axis makes everything go screwey

    chl_array_fill_val = chl_array.fill_value
    #TODO insert the solar zenith angle establishment, only if handed a modelled input (so add flag for --modelled-input)
    #Here we make a very big assumption that index 0 is january of the year being considered in the chlorophyll array.
    if modelled_threshold:
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
                if zen >= 75:
                    date_zeniths.append(0)
                else:
                    date_zeniths.append(1)
                true_zen.append(zen)
            true_zens.append(true_zen)
            date_masks.append(date_zeniths)
        
        temp_chl_array = chl_array


        ds.createVariable('zen', 'float32', dimensions=['TIME', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)
        ds.variables['zen'].setncattr("units", "degrees")
        for year in range(0,  chl_array.shape[0], date_seperation):
            for index, date_mask in enumerate(date_masks):
                for row, row_mask in enumerate(date_mask):
                    if not row_mask:
                        if LAT_IDX == 2:
                            temp_chl_array.mask[year + index,0,row,:] = True
                        if LAT_IDX == 3:
                            temp_chl_array.mask[year + index,0,:,row] = True
                    ds.variables['zen'][index,row,:] = true_zens[index][row]

        #TODO add gap filling, --gap-filling with a few choices for interpolation options, if not specified then don't do it at all


    #create the threshold var, this should be done on the original chlorophyll data, not the filled variable
    if modelled_threshold:
        logger.info("chl_median")
        temp_chl_array = numpy.ma.masked_where(temp_chl_array == chl_array_fill_val, temp_chl_array)
        chl_median = numpy.ma.median(temp_chl_array,axis=0)
        logger.info("median value: {}".format(chl_median))
        logger.info("median threshold: {}".format(median_threshold))
    else:
        logger.info("chl_median")
        chl_median = numpy.ma.median(chl_array,axis = 0)
        logger.info("median value: {}".format(chl_median))
        logger.info("median threshold: {}".format(median_threshold))
    
    if modelled_threshold:
        std_dev = numpy.std(temp_chl_array, axis=0)
    else:
        std_dev = numpy.std(chl_array, axis=0)
    
    if modelled_threshold:
        max_means, min_means, perc_val_change = numpy.apply_along_axis(find_maxes, 0,temp_chl_array)
    else:
        max_means, min_means, perc_val_change = numpy.apply_along_axis(find_maxes, 0,chl_array)
        

    if modelled_threshold:
        chl_array = temp_chl_array

    #! after estimating the median, i apply a filling to the chl time-series, e.g., window of 8 time steps, but it would be good to be able to choose other width of the window, e.g., 5 time-steps or 3 time-steps...
    if do_fill_chl:
        #create a 9 x 9 x 9 filter such that central value will be average over time and space
        k = numpy.ones((5,5,5))
        logger.info("pre filling shape {}".format(chl_array.shape))
        chl_reshape = chl_array.copy()
        chl_reshape.shape = (chl_array.shape[0],chl_array.shape[2],chl_array.shape[3])
        #extend values beyond array edge, interpolate over nan values as opposed to filling them with the default fill value (-999999 or something similar)
        filled_chl = astropy.convolution.convolve(chl_reshape, k, boundary='extend', nan_treatment='interpolate', mask=chl_reshape.mask)/k.sum()
        logger.info("post filling shape {}".format(filled_chl.shape))
        filled_chl.shape = chl_array.shape
        logger.info("post filling and reshape shape {}".format(filled_chl.shape))

        filled_chl = numpy.ma.masked_invalid(filled_chl)
        filled_chl_temp = chl_array.copy()
        filled_chl_temp[filled_chl_temp.mask == True] = filled_chl[filled_chl_temp.mask == True] 
        filled_chl = filled_chl_temp
        filled_chl = numpy.ma.masked_invalid(filled_chl)
    else:
        filled_chl = chl_array

    logger.info("filled_chl has mask? {}".format(filled_chl.mask.any()))

    ds.createVariable('filled_chl', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)

    #! here it would be good to give the option to select the value of the median threshold, e.g., median plus 5%, 10%, 15%...



    #get anomalies
    #anomaly needs to be changed here to reflect what we are actually working with, which is normalised trend I guess?
    logger.info("anomaly")

    if std_dev_anomaly:
        anomaly = filled_chl - (std_dev*float(std_dev_threshold))
    elif max_means_anomaly:
        anomaly = filled_chl - (max_means - (max_means * max_means_threshold))
    elif relative_max_anomaly:
        #dynamic threshold calc
        anomaly = filled_chl - (max_means - (max_means * (perc_val_change * max_means_threshold)))
    elif relative_median_anomaly:
        #dynamic threshold calc
        #set 0.25 to be parameter
        anomaly = filled_chl - (chl_median + (chl_median * (perc_val_change * 0.25)))
    else:
        anomaly = filled_chl - (chl_median*median_threshold)

    if debug_pixel_main[2] == chunk_idx:
        logger.info(anomaly[:5,:,debug_pixel_main[0],debug_pixel_main[1]])
    #anomaly_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_anomaly"), mode="w+", shape=anomaly.shape, dtype=USE_DTYPE)
    #anomaly_map[:] = anomaly[:]

    #anomaly = anomaly_map
    ds.createVariable('anomaly', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)
    ds.variables['anomaly'].setncattr("units", "mg chl m^3")
    ds.variables['anomaly'][:] = anomaly[:]

    anomaly_map = None
    #need to ditch any empty entries here as they interefere with the cumsum

    #get cumsum of anomalies
    logger.info("chl cumsum")
    chl_cumsum = numpy.ma.cumsum(anomaly,axis=0)
    if debug_pixel_main[2] == chunk_idx:
        logger.info(chl_cumsum[:5,:,debug_pixel_main[0],debug_pixel_main[1]])
    #chl_cumsum_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_cumsum"), mode="w+", shape=chl_cumsum.shape, dtype=USE_DTYPE)
    #chl_cumsum_map[:] = chl_cumsum[:]
    #chl_cumsum = chl_cumsum_map
    anomaly = None
    chl_cumsum_map = None

    ds.createVariable('chl_cumsum', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)
    ds.variables['chl_cumsum'].setncattr("units", "mg chl m^3")
    ds.variables['chl_cumsum'][:] = chl_cumsum[:]

    #get centered derivative        
    logger.info("chl der")
    chl_der = numpy.apply_along_axis(centered_diff_derivative, 0, chl_cumsum)
    if debug_pixel_main[2] == chunk_idx:
        logger.info(chl_der[:5,:,debug_pixel_main[0],debug_pixel_main[1]])
    #chl_der_map = numpy.memmap(os.path.join(numpy_storage, str(chunk), "chl_der"), mode="w+", shape=chl_der.shape, dtype=USE_DTYPE)
    #chl_der_map[:] = chl_der[:]
    #chl_der = chl_der_map
    chl_cumsum = None
    chl_der_map = None
    
    ds.createVariable('chl_der', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)
    ds.variables['chl_der'].setncattr("units", "mg chl m^3")
    ds.variables['chl_der'][:] = chl_der[:]

    #boxcar filter with width of 3 (sbx) should be something like this:
    logger.info("chl sbx")
    chl_boxcar = numpy.apply_along_axis(lambda m: numpy.convolve(m, numpy.ones(3)/3, mode='full'), axis=0, arr=chl_der)
    if debug_pixel_main[2] == chunk_idx:
        logger.info(chl_boxcar[:5,:,debug_pixel_main[0],debug_pixel_main[1]])
    fill_arr = numpy.ma.masked_array(numpy.zeros((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])),
                                     mask=numpy.ones((1,1,chl_boxcar.shape[2],chl_boxcar.shape[3])))
    logger.info("chl shape after boxcar")
    logger.info(chl_boxcar.shape)

    #pad the boxcar back to january
    new_timesteps = chl_boxcar.shape[0] - chl_array.shape[0]
    additional_steps = new_timesteps / 2

    logger.info("after chl preparation {} new timesteps were created, original shape: {}, new shape: {}".format(new_timesteps, chl_array.shape[0], chl_boxcar.shape[0]))
    logger.info("this means there are {} new timesteps at the beginning and end of the chl boxcar array".format(additional_steps))

    if not missing_chl_dates_at_start - math.floor(additional_steps) < 0:
        logger.info("padding {} time steps from january {} in chl boxcar".format(missing_chl_dates_at_start - math.floor(additional_steps), START_YEAR))
        start_fill_arrays = [fill_arr for i in range(math.floor(additional_steps), missing_chl_dates_at_start)]
    else:
        start_fill_arrays = missing_chl_dates_at_start - math.floor(additional_steps)

    if not additional_steps.is_integer():
        logger.warning("as the number of new timesteps is not a whole number will pad one additional step to the beginning of the chl boxcar array")

    if not missing_chl_dates_at_end - math.ceil(additional_steps) < 0:
        logger.info("padding {} time steps to december in chl boxcar".format(missing_chl_dates_at_end - math.ceil(additional_steps), START_YEAR))
        end_fill_arrays = [fill_arr for i in range(math.ceil(additional_steps), missing_chl_dates_at_end)]
    else:
        end_fill_arrays = missing_chl_dates_at_end - math.ceil(additional_steps)

    if not isinstance(end_fill_arrays, list):
        logger.info("removing {} values from boxcar end".format(end_fill_arrays))
        chl_boxcar = chl_boxcar[0:chl_boxcar.shape[0] - abs(end_fill_arrays)]
        end_fill_arrays=[]
        logger.info(chl_boxcar.shape)

    if not isinstance(start_fill_arrays, list):
        logger.info("removing {} values from boxcar start".format(start_fill_arrays))
        chl_boxcar = chl_boxcar[0+abs(start_fill_arrays):chl_boxcar.shape[0]]
        start_fill_arrays=[] 
        logger.info(chl_boxcar.shape)

    

    if not (chl_boxcar.shape[0] + len(start_fill_arrays) + len(end_fill_arrays)) % date_seperation_per_year == 0:
        logger.error("After padding, the chlorophyll boxcar end product did not divide equally!")
        logger.error((chl_boxcar.shape[0] + len(start_fill_arrays) + len(end_fill_arrays), date_seperation_per_year, (chl_boxcar.shape[0] + len(start_fill_arrays) + len(end_fill_arrays)) % date_seperation_per_year))
        raise Exception("After padding, the chlorophyll boxcar end product did not divide equally!")
    
    #pad the boxcar back to january    
    chl_boxcar = numpy.ma.concatenate(start_fill_arrays + [chl_boxcar] + end_fill_arrays, axis=0)

    #now do the original chlorophyll array
    logger.info("padding {} time steps from january {} in chl array".format(missing_chl_dates_at_start, START_YEAR))
    start_fill_arrays = [fill_arr for i in range(0, missing_chl_dates_at_start)]
    logger.info("padding {} time steps to december in chl boxcar".format(missing_chl_dates_at_end, START_YEAR))
    end_fill_arrays = [fill_arr for i in range(0, missing_chl_dates_at_end)]

    if not (filled_chl.shape[0] + len(start_fill_arrays) + len(end_fill_arrays)) % date_seperation_per_year == 0:
        logger.error("After padding, the chlorophyll boxcar end product did not divide equally!")
        raise Exception("After padding, the chlorophyll boxcar end product did not divide equally!")

    logger.info("mask pre concatenation")
    logger.info(filled_chl.mask.any())
    filled_chl = numpy.ma.concatenate(start_fill_arrays + [filled_chl] + end_fill_arrays)
    filled_chl = numpy.ma.masked_where(filled_chl == chl_array_fill_val, filled_chl)

    logger.info("mask post concatenation")
    logger.info(filled_chl.mask.any())

    ds.variables['filled_chl'].setncattr("units", "mg chl m^3")
    ds.variables['filled_chl'][:] = filled_chl[:]
    logger.info("chl array shape after padding")
    logger.info(filled_chl.shape)

    logger.info("chl boxcar shape after padding")
    logger.info(chl_boxcar.shape)
    
    ds.createVariable('chl_boxcar', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=True)
    ds.variables['chl_boxcar'].setncattr("units", "mg chl m^3")
    ds.variables['chl_boxcar'][:] = chl_boxcar[:]

    ds.close()

    logger.info("chl prep complete")
    return chl_boxcar.shape, USE_DTYPE, chl_median, std_dev, max_means, perc_val_change, filled_chl



def create_phenology_netcdf(chl_lons, chl_lats, output_shape=None,name="phenology_{}.nc".format(datetime.datetime.now().strftime("%H%M")), date=False, median=None, std=None, max_means=None, perc_val_change=None):
    """
    Creates the skeleton of the netcdf file to be used by write_to_output_netcdf, all of this is metadata.
    """
    global chunk_size
    global zlib_compression
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
    ds.createVariable('date_start1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=True)
    week_descriptor = 'weeks from {}'.format(REF_MONTH)
    #description = ' data between {} and {} for years {} to {}'.format(REF_MONTH, END_MONTH, START_YEAR, END_YEAR)
    ds.variables['date_start1'].setncattr("units", week_descriptor)
    ds.createVariable('date_max1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['date_max1'].setncattr("units", week_descriptor)
    ds.createVariable('date_end1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['date_end1'].setncattr("units", week_descriptor)
    ds.createVariable('duration1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['duration1'].setncattr("units", week_descriptor)
    ds.createVariable('date_start2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['date_start2'].setncattr("units", week_descriptor)
    ds.createVariable('date_max2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['date_max2'].setncattr("units", week_descriptor)
    ds.createVariable('date_end2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['date_end2'].setncattr("units", week_descriptor)
    ds.createVariable('duration2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['duration2'].setncattr("units", week_descriptor)
    ds.createVariable('max_val1', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.createVariable('max_val2', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['max_val2'].setncattr("units", "mgChl/m3")
    ds.variables['max_val1'].setncattr("units", "mgChl/m3")
    ds.createVariable('total_blooms', 'float32', dimensions=DIM_ORDER,fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['total_blooms'].setncattr("units", "observations")
    ds.createVariable('probability', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['probability'].setncattr("units", "likelihood")
    ds.createVariable('median_chlorophyll', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.createVariable('chlorophyll_std_dev', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.createVariable('max_mean', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.createVariable('perc_val_change', 'float32', dimensions=DIM_ORDER[1:4],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['median_chlorophyll'].setncattr("units", "mg chl m^3")
    ds.variables['chlorophyll_std_dev'].setncattr("units", "mg chl m^3")
    ds.variables['max_mean'].setncattr("units", "mg chl m^3")
    ds.variables['perc_val_change'].setncattr("units", "percent")
    ds.setncattr("generation command", str(" ".join(sys.argv)))
    ds.setncattr("run location", str(os.getcwd()))
    ds.setncattr("reverse_search", reverse_search)
    ds.setncattr("date_zero", date_zero_datetime.strftime("%Y/%m/%d"))
    ds.setncattr("steps_per_year", date_seperation_per_year)
    if isinstance(median, numpy.ndarray):
        ds.variables['median_chlorophyll'][:] = median
    if isinstance(std, numpy.ndarray):
        ds.variables['chlorophyll_std_dev'][:] = std
    ds.variables['max_mean'][:] = max_means
    ds.variables['perc_val_change'][:] = perc_val_change
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

def create_intermediate_netcdf(output_name, chl_lons, chl_lats):
    global chunk_size
    global zlib_compression
    ds = nc.Dataset(output_name,'w',format='NETCDF4_CLASSIC')
    ds.createDimension('LONGITUDE', chl_lons.shape[0])
    ds.createDimension('LATITUDE', chl_lats.shape[0])
    ds.createDimension('DEPTH', 1)
    ds.createDimension('TIME', None)
    ds.createVariable('LATITUDE', 'float64', dimensions=['LATITUDE'], zlib=zlib_compression,)
    ds.variables['LATITUDE'].setncattr("units", "degrees_north")
    ds.variables['LATITUDE'][:] = chl_lats
    ds.createVariable('LONGITUDE', 'float64', dimensions=['LONGITUDE'], zlib=zlib_compression,)
    ds.variables['LONGITUDE'].setncattr("units", "degrees_east")
    ds.variables['LONGITUDE'][:] = chl_lons
    ds.createVariable('DEPTH', 'float32', dimensions=['DEPTH'], zlib=zlib_compression,)
    ds.variables['DEPTH'].setncattr("units", "meters")
    ds.variables['DEPTH'].setncattr("positive", "down")
    ds.variables['DEPTH'][:] = [0.1]
    ds.createVariable('TIME', 'float32', dimensions=['TIME'], zlib=zlib_compression,)
    ds.variables['TIME'].setncattr("units", "years")
    ds.createVariable('zen', 'float32', dimensions=['TIME', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['zen'].setncattr("units", "degrees")
    ds.createVariable('filled_chl', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['filled_chl'].setncattr("units", "mg chl m^3")
    ds.createVariable('anomaly', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['anomaly'].setncattr("units", "mg chl m^3")
    ds.createVariable('chl_cumsum', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['chl_cumsum'].setncattr("units", "mg chl m^3")
    ds.createVariable('chl_der', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['chl_der'].setncattr("units", "mg chl m^3")
    ds.createVariable('chl_boxcar', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['chl_boxcar'].setncattr("units", "mg chl m^3")
    ds.createVariable('sst_boxcar', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['sst_boxcar'].setncattr("units", "degrees celsius")
    ds.createVariable('sst_der', 'float32', dimensions=['TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'],fill_value=FILL_VAL, zlib=zlib_compression)
    ds.variables['sst_der'].setncattr("units", "degrees celsius")
    ds.close()


def total_blooms_to_probability(array_like):
    return numpy.count_nonzero(array_like == 2) / array_like.size

def get_multi_year_two_blooms_output(output_name, chunk, chl_shape, chl_dtype, chl_data, sst_shape, sst_dtype, date_seperation_per_year=47, start_date=0, reverse_search=20, reference_index=0, out_netcdf=output_location, out_date_netcdf=output_location_date, chunk_idx=None, end_date=False):
    #this all works on the assumption the axis 0 is time
    global debug_pixel_main
    global extend_array

    ds = nc.Dataset(output_name.replace(".nc", "_intermediate_products.nc".format(chunk)),'r',format='NETCDF4_CLASSIC')
    
    logger.info("reading variables")
    chl_boxcar = ds.variables["chl_boxcar"][:]

    if not end_date:
        end_date = chl_boxcar.shape[0]
        if extend_array:
            end_date = end_date - date_seperation_per_year
    logger.info(sst_dtype)
    sst_der = ds.variables["sst_der"][:]
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
    logger.info("output array shape expecting to write = {}".format((chl_data.shape[2],chl_data.shape[3], int((end_date - start_date) // date_seperation_per_year), 2,5)))
    year_true_start_end_array = numpy.ndarray((chl_data.shape[2],chl_data.shape[3], int((end_date - start_date) // date_seperation_per_year), 2,5))

    logger.info("number of years expecting to write = {}".format(int((end_date -start_date) // date_seperation_per_year)))
    total_blooms = numpy.ndarray((chl_data.shape[2],chl_data.shape[3], int((chl_data.shape[0] - start_date) // date_seperation_per_year)))
    year_true_start_end_array.fill(FILL_VAL)
    blooms_by_date = year_true_start_end_array.copy()
    total_blooms.fill(FILL_VAL)
    total_blooms_date = total_blooms.copy()
    logger.info("doing sst initiations and correction")
    logger.debug("start date : {}".format(start_date))
    for ix, iy in tqdm.tqdm(numpy.ndindex(chl_data.shape[2], chl_data.shape[3]), total=(chl_data.shape[2] * chl_data.shape[3]), disable=True):
        expected_exceptions = []
        try:
            verbose=False
            if iy == debug_pixel_main[0] and ix == debug_pixel_main[1] and chunk_idx == debug_pixel_main[2]:
                if logger.level <= logging.INFO:
                    logger.setLevel(logging.DEBUG)
                    logger.info("debug pixel {} {} encountered".format(ix,iy))
                    verbose = True
            else:
                logger.setLevel(default_logging)
                pass
            results = match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],
                                                     chl_boxcar[:,:,ix,iy],
                                                     chl_data[:,:,ix,iy], 
                                                     date_seperation_per_year, 
                                                     reverse_search, 
                                                     verbose=verbose, 
                                                     start_date=start_date, 
                                                     reference_date=reference_index,
                                                     processing_end_date=end_date)
            year_true_start_end_array[ix,iy] = results[0]
            total_blooms[ix,iy] = results[1]
            blooms_by_date[ix, iy] = results[2]
            total_blooms_date[ix,iy] = results[3]

            logger.debug("end duration array")
            logger.debug(year_true_start_end_array[ix,iy])
        except Exception as e:
            if e in expected_exceptions:
                pass
            expected_exceptions.append(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            logger.error(e)
            logger.error(repr(e))
            logger.error(results)
            logger.error((ix, iy))
            logger.error(exc_tb.tb_lineno)
            match_start_end_to_solar_cycle(sst_der[:,:,ix,iy],
                                           chl_boxcar[:,:,ix,iy], 
                                           chl_data[:,:,ix,iy], 
                                           date_seperation_per_year, 
                                           reverse_search, 
                                           verbose=False, 
                                           start_date=start_date, 
                                           reference_date=reference_index,
                                           end_date=end_date)

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

def chunk_by_chunk_handler(chunk_idx, chunk):
    global default_logging
    global debug_pixel_only
    global debug_chunk
    if (not debug_chunk == chunk_idx and debug_pixel_only):
        return
    try:
        log = log_lock.acquire(block=False)
        if log:
            handler.setLevel(logging.DEBUG)
            logger.debug("set logging to debug on filehandler")
        logger.debug("processing chunk {} with bounds {}".format(chunk_idx, str(chunk)))
        chunk_by_chunk(chunk_idx, chunk)
        if log:
            try:
                log_lock.release()
                handler.setLevel(default_logging)
            except ValueError:
                pass

    except Exception as e:
        traceback.print_exc()
        log_lock.release()
        handler.setLevel(default_logging)
        logger.error(e)
    raise e

def get_chl_ngd(start_date, date_seperation_per_year, chl_array):
    ngd_arr = []
    for year in range(start_date, chl_array.shape[0], date_seperation_per_year):
        year_chl = chl_array[year:year+date_seperation_per_year]
        year_ngd = numpy.ma.count_masked(year_chl, axis=0)
        ngd_arr.append(year_ngd)
    ngd_arr = numpy.dstack(ngd_arr)
    return ngd_arr

def find_nearest(array, value):
    array = numpy.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
    return idx


def chunk_by_chunk(chunk_idx, chunk):
    global chl_ds
    global args
    global start_date
    global debug_chunk
    global do_only_debug_chunk

    chunk_start_date = start_date

    #if we have been instructed to skip everything then do so
    if do_only_debug_chunk and not chunk_idx == debug_chunk:
        return True

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
    chl_array = numpy.ma.masked_array(chl_array, numpy.isnan(chl_array))
    chl_array = numpy.ma.masked_where(chl_array == chl_array.fill_value, chl_array)
    chl_array = numpy.ma.masked_where(chl_array <= 0, chl_array)
    chl_array = numpy.ma.masked_invalid(chl_array)
    if chl_array.mask.all():
        logger.info(numpy.isnan(chl_array).all())
        logger.info("skipping as empty")
        chl_lock.release()
        #output empty netcdf
        empty_med = numpy.empty((chl_lats.shape[0], chl_lons.shape[0]), dtype=chl_array.dtype) if LAT_IDX >LON_IDX else numpy.empty(shape=(chl_lons.shape[0], chl_lats.shape[0]), dtype=chl_array.dtype)
        create_phenology_netcdf(chl_lons, chl_lats, [1,1], output.replace(".nc", "_by_maxval.nc"), median=empty_med, std=empty_med, max_means=empty_med)
        create_phenology_netcdf(chl_lons, chl_lats, [1,1], output.replace(".nc", "_by_date.nc"), date=True, median=empty_med, std=empty_med, max_means=empty_med)
        return True
    chl_lock.release()
    if len(chl_array.shape) == 3:
        logger.info("reshaping chl to {}".format((chl_array.shape[0], 1, chl_array.shape[1], chl_array.shape[2])))
        chl_array.shape = (chl_array.shape[0], 1, chl_array.shape[1], chl_array.shape[2])
    
    logger.info("making temp storage")

    """
    if not (chl_array.shape[2] == chl_lats.shape[0] and chl_array.shape[3] == chl_lons.shape[0]):
        logger.info("adjusting to flip lat and lon")
        chl_array.shape = (chl_array.shape[0], chl_array.shape[1], chl_array.shape[3], chl_array.shape[2])
    """
    if args.extend_chl_data:
        chl_array, chunk_start_date = extend_array(chl_array, date_seperation_per_year, chunk_start_date)
        logger.info("start date after extension: {}".format(chunk_start_date))

    chl_shape, chl_dtype, chl_median, chl_std_dev, chl_max_means, perc_val_change, filled_chl = prepare_chl_variables(chl_array, 
                                                                          chunk_idx, 
                                                                          date_seperation_per_year, 
                                                                          chl_lats, 
                                                                          chl_lons, 
                                                                          modelled_threshold=args.modelled_threshold, 
                                                                          median_threshold=med_thresh, 
                                                                          median_filename=chl_filename.replace(".nc", 
                                                                                                               "_median_{}_chunk{}.nc".format(time_of_run,
                                                                                                                                              chunk_idx)
                                                                                                              ),
                                                                          chunk_idx=chunk_idx,
                                                                          std_dev_anomaly=do_std_dev,
                                                                          std_dev_threshold=std_dev_threshold,
                                                                          relative_max_anomaly=do_rel_max,
                                                                          max_means_threshold=max_means_threshold,
                                                                          max_means_anomaly=do_max_means,
                                                                          relative_median_anomaly=do_rel_med_anomaly,
                                                                          do_fill_chl=fill_chl,
                                                                          output_name=output,
                                                                          date_seperation_per_year=date_seperation_per_year)
    logger.info("chl_shape: {}".format(chl_shape))
    logger.info("chl_dtype: {}".format(USE_DTYPE))

    chl_ngd = get_chl_ngd(chunk_start_date, date_seperation_per_year, filled_chl)

    if args.sst_location:
        sst_lock.acquire()
        logger.info("sst file provided, reading array")
        logger.info("only one file found, assuming full stack of observations")
        sst_ds = nc.Dataset(args.sst_location)
        sst_variable = [x for x in sst_ds.variables if args.sst_var in x.lower()][0]
        sst_lon_var, sst_lat_var, sst_time_var = ds_to_dim_vars(sst_ds)
        sst_lons = sst_ds.variables[sst_lon_var][:]
        sst_lats = sst_ds.variables[sst_lat_var][:]
        sst_time = sst_ds.variables[sst_time_var][:]
        SST_LAT_IDX = sst_ds.variables[sst_variable].dimensions.index(sst_lat_var)
        SST_LON_IDX = sst_ds.variables[sst_variable].dimensions.index(sst_lon_var)
        SST_TIME_IDX = sst_ds.variables[sst_variable].dimensions.index(sst_time_var)
        sst_slc = [slice(None)] * len(sst_ds.variables[sst_variable].shape)
        sst_slc[SST_LON_IDX] = slice(x[0], x[1])
        sst_slc[SST_LAT_IDX] = slice(y[0], y[1])
        sst_slc[SST_TIME_IDX] = slice(args.sst_start_index, sst_time.shape[0])
        sst_array = sst_ds.variables[sst_variable][sst_slc]
        logger.info("sst array shape at read:{}".format(sst_array.shape))
        try:
            sst_fill_val =  sst_ds.variables[sst_variable]._FillValue
        except:
            sst_fill_val = FILL_VAL
        
        if not numpy.ma.is_masked(sst_array):
            logger.debug("sst is not a masked array, masking it")
            sst_array = numpy.ma.masked_array(sst_array, numpy.isnan(sst_array))
            sst_array = numpy.ma.masked_where(sst_array == sst_fill_val, sst_array)
        sst_lock.release()

        if args.extend_sst_data:
            logger.debug("extending sst array by repeating year at start and end")
            sst_array, _ = extend_array(sst_array, date_seperation_per_year, start_date)

        if len(sst_array.shape) == 3:
            logger.info("reshaping sst to {}".format((sst_array.shape[0], 1, sst_array.shape[1], sst_array.shape[2])))
            sst_array.shape = (sst_array.shape[0], 1, sst_array.shape[1], sst_array.shape[2])

        sst_shape, sst_dtype = prepare_sst_variables(sst_array, chunk_idx, skip=args.skip_sst_prep,chunk_idx=chunk_idx, output_name=output)
        logger.info("sst_shape: {}".format(sst_shape))
        logger.info("sst_dtype: {}".format(USE_DTYPE))
        sst_array = None

    if sst_shape[2:] != chl_shape[2:]:
        logger.error("sst and chlorophyll x,y array shapes do not match got:")
        logger.error("chlorophyll: {}".format(chl_shape[2:]))
        logger.error("sst: {}".format(sst_shape[2:]))
        logger.error("quitting!")
        sys.exit()
    if not isinstance(filled_chl.mask, numpy.ndarray):
        logger.error(filled_chl.mask)
        logger.error("chl data has no mask for its variable data, this will completely break our logic")
        sys.exit()
        return
    #simple regridding
    create_phenology_netcdf(chl_lons, chl_lats, chl_shape, output.replace(".nc", "_by_maxval.nc"), median=chl_median, std=chl_std_dev, max_means=chl_max_means, perc_val_change=perc_val_change)
    create_phenology_netcdf(chl_lons, chl_lats, chl_shape, output.replace(".nc", "_by_date.nc"), date=True, median=chl_median, std=chl_std_dev, max_means=chl_max_means, perc_val_change=perc_val_change)
    logger.info("using start date {}".format(chunk_start_date))
    get_multi_year_two_blooms_output(output,
                                    chunk_idx,
                                    chl_shape,
                                    chl_dtype, 
                                    filled_chl, 
                                    sst_shape, 
                                    sst_dtype, 
                                    date_seperation_per_year=date_seperation_per_year, 
                                    start_date=chunk_start_date, 
                                    reverse_search=reverse_search,
                                    reference_index=ref_index,
                                    out_date_netcdf=output.replace(".nc", "_by_date.nc"),
                                    out_netcdf=output.replace(".nc", "_by_maxval.nc"),
                                    chunk_idx=chunk_idx)
    
    logger.setLevel(default_logging)
    return True

def ds_to_dim_vars(ds):
    lon_var, lat_var, time_var = (None,)*3
    try:
        try:
            lon_var = [x for x in ds.variables if "lon" in x.lower()][0]
        except:
            logger.info("trying to get lat/lon from standard name")
            for va in ds.variables:
                if "standard_name" in ds.variables[va].__dict__.keys():
                    if ds.variables[va].__dict__["standard_name"] == "longitude":
                        lon_var= va

        try:
            lat_var = [x for x in ds.variables if "lat" in x.lower()][0]
        except:
            logger.info("trying to get lat/lon from standard name")
            for va in ds.variables:
                if "standard_name" in ds.variables[va].__dict__.keys():
                    if ds.variables[va].__dict__["standard_name"] == "latitude":
                        lat_var= va
        try:
            time_var = [x for x in ds.variables if "time" in x.lower()][0]
        except:
            logger.info("trying to get lat/lon from standard name")
            for va in ds.variables:
                if "standard_name" in ds.variables[va].__dict__.keys():
                    if ds.variables[va].__dict__["standard_name"] == "time":
                        time_var= va
    except:
        pass

    if not all([lon_var,lat_var,time_var]):
        logger.error("Unable to identify matching variables for longitude, latitude or time, to enable this functionality you must update the contents of your input netcdf files. Dimensions should have a matching variable to describe them, and each of those most either have a standard_name attribute or name that contains 'lat' 'lon' or 'time'. Variables in this file are: {}".format(", ".join(ds.variables)))
        raise ValueError("could not find one of lat, lon or time")
    return lon_var,lat_var,time_var

def get_ds_time_data(ds, time_var, begin_date=False, var="chl"):
    global date_seperation_per_year
    logger.info("ds time var shape: {}".format(ds[time_var].shape[0]))
    try:
        if not begin_date:
            dts = nc.num2date(ds[time_var][:], ds[time_var].units)
            if not date_seperation_per_year:
                date_seperation_per_year = math.ceil(365 / round(((dts[-1] - dts[0]).days) / ds[time_var][:].shape[0]))
                logger.info("estimated seperation per year as {}".format(date_seperation_per_year))
            ref_date = datetime.datetime(year=dts[0].year, month=1, day=1) + datetime.timedelta(days=(365 * (start_date +ref_index / date_seperation_per_year)))
            start_datetime = datetime.datetime(year=dts[0].year, month=1, day=1) + datetime.timedelta(days=(365 * (start_date / date_seperation_per_year)))
            init_date = dts[0]
            end_date = dts[-1]
        else:
            if not date_seperation_per_year:
                #we cant guess
                logger.error("please specify number of indexes per year with --date_seperation_per_year (e.g. --date_seperation_per_year 46 would mean 46 steps per year")
                sys.exit()
            init_date = datetime.datetime.strptime(begin_date, '%d/%m/%Y')
            start_datetime = datetime.datetime(year=init_date.year,month=1,day=1) + datetime.timedelta(days=(365 * (start_date / date_seperation_per_year)))
            ref_date =  datetime.datetime(year=init_date.year,month=1,day=1) + datetime.timedelta(days=(365 * (start_date +ref_index / date_seperation_per_year)))
            logger.info("number of days between start and end")
            logger.info(ds[time_var].shape[0] * math.ceil(365 / date_seperation_per_year))
            end_date = init_date + datetime.timedelta(days=(ds[time_var].shape[0] * math.ceil(365 / date_seperation_per_year)))
        logger.info(ref_date)
        REF_MONTH = ref_date.strftime("%d %B")
        START_YEAR = start_datetime.year
        logger.info("reference date selected: "+ref_date.strftime("%d/%m/%Y"))
        logger.info("start date selected: "+start_datetime.strftime("%d/%m/%Y"))
        logger.info("reference month: "+ REF_MONTH)
        logger.info("start year: "+str(START_YEAR))


        logger.info("getting {} missing timesteps".format(var))
        start_difference = init_date - datetime.datetime(year=init_date.year, month=1, day=1)
        logger.info(end_date)
        end_difference =  datetime.datetime(year=end_date.year, month=12, day=31) - end_date
        missing_dates_at_start = start_difference.days// math.ceil(365 / date_seperation_per_year)
        missing_dates_at_end = end_difference.days// math.ceil(365 / date_seperation_per_year)
        logger.info("missing steps to january at start of file: {}".format(missing_dates_at_start))
        logger.info("missing steps to december at end of file: {}".format(missing_dates_at_end))
    except Exception as e:
        logger.error(e)
        logger.error("could not identify {} start date, please specify the first date you expect in the file (time index 0, not first date index) with --{}_begin_date e.g. --chl_begin_date 01/01/1997".format(var, var))
        sys.exit()
    
    return REF_MONTH, START_YEAR, missing_dates_at_start, missing_dates_at_end, ds[time_var].shape[0], start_datetime

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("chl_location", nargs='+', help="A chlorophyll file or list of chlorophyll files. These must have a variable called chl, depth, lat and lon (or contain those strings) which will be used to define the array shape")
    parser.add_argument("--chl_var", help="specify chl variable, otherwise is guessed based on variables that contain'chl' in the name", default="chl", required=False)
    parser.add_argument("--sst_location", help="An sst file, or glob of files (e.g. sst*.nc) that matches the chlorophyll observations. if it does not match some interpolation can be attempted with --sst_date_interp", required=False)
    parser.add_argument("--sst_var", help="specify the sst variable name, otherwise is guessed based on variables containing 'sst'", default="sst", required=False)
    parser.add_argument("--output_folder", help="output folder, if not specified defaults to (chlorophyll_filename)_phenology.nc in the current folder", default=None, required=False)
    parser.add_argument("--date_seperation_per_year", help="how many temporal observations we have in a year, if not specified will be guessed", default=0, required=False)
    parser.add_argument("--first_date_index", help="Specify if the first date you want to include is not the first date present in the date stack.", default=1, required=False)
    parser.add_argument("--reference_index", help="Date index in relation to first_date_index to use as the reference from which weeks are measured in the phenology output. If specified is used as first_date_index + reference_index, if not set will measure from first_date_index (default 1)", default=1, required=False)
    parser.add_argument("--chl_climatology", action="store_true", help="extend the input chlorophyll array by creatingidentical copied for the year previous and year next", default=0, required=False)
    parser.add_argument("--reverse_search", default=False, help="specify the number of observations to search in the previous year, if not specified will be calculated as a representation of 100 days (date_seperation_per_year / 0.27).", required=False)
    parser.add_argument("--skip_sst_prep", action="store_true", default=False, help="skip sst preparation, instead the program will assume that the sst input is already in a state where low/high period fluctuation can be identified", required=False)
    parser.add_argument("--median_threshold", default=MEDIAN_THRESHOLD_DEFAULT, help="change median threshold, set as percentage value e.g. 20 = 20", required=False)
    parser.add_argument("--modelled_threshold", action='store_true',default=False, help="test for solar zenith of inputs", required=False)
    parser.add_argument("--std_deviation", action='store_true',default=False, help="Use standard deviation instead of median median chlorophyll to generate anomalies", required=False)
    parser.add_argument("--std_deviation_threshold", default=STD_DEV_THRESHOLD_DEFAULT, help="change std deviation overage, default is 2", required=False)
    parser.add_argument("--max_means", action='store_true',default=False, help="Use mean of maximums instead of median chlorophyll to generate anomalies", required=False)
    parser.add_argument("--rel_max_means", action='store_true',default=False, help="Use mean of maximums instead of median chlorophyll to generate anomalies", required=False)
    parser.add_argument("--rel_median", action='store_true',default=False, help="use a threshold generated by scaling from min/max instead of a fixed value as with median threshold", required=False)
    parser.add_argument("--max_means_threshold", default=MAX_MEANS_THRESHOLD_DEFAULT, help="change amount subtracted from max means for anomaly thresholding, default is 20\% of max means", required=False)
    #probably needs a better description!
    #give options for both since we might end up in a situation where sst is 3 years and chl is one year (or vice versa)
    parser.add_argument("--extend_chl_data", default=False, action="store_true", help="extends chlorophyll by copying the (central) chl array to the year previous and year next")
    parser.add_argument("--extend_sst_data", default=False, action="store_true", help="extends sea surfaace temperature by copying the (central) chl array to the year previous and year next")
    parser.add_argument("--start_year", default=0, help="What year to use as the start point for metadata output, if not specified will use 0, affects no processing.")
    parser.add_argument("--debug_pixel",  nargs='+', default=None, type=float, help="pixle in x, y (lat, lon), these entries are 0 indexed.")
    parser.add_argument("--debug_pixel_only",  default=False, action="store_true", help="don't process anything that isnt the debug pixel")
    parser.add_argument("--sst_start_index", type=int, default=0)
    parser.add_argument("--sst_end_index", type=int, default=-1)
    parser.add_argument("--chl_begin_date", default=None, help="date relating to index 0 of a file in format dd/mm/yyyy")
    parser.add_argument("--sst_begin_date", default=None, help="date relating to index 0 of a file in format dd/mm/yyyy")
    parser.add_argument("--no_thread", action='store_true', default=False, help="Don't use threading")
    parser.add_argument("--no_compress", action='store_true', default=False, help="Don't compress the data - for large files this can be the difference of a days processing :(")
    parser.add_argument("--no_logfile", action='store_true', default=False, help="Don't create a logfile")
    parser.add_argument("--stitch_only", action='store_true', default=False, help="Don't create a logfile")
    parser.add_argument("--no_delete", action='store_true', default=False, help="Don't delete chunks")
    parser.add_argument("--no_chl_fill", action='store_true', default=False, help="Don't fill the chlorophyll variable")
    parser.add_argument("--debug_chunk_only", action='store_true', default=False, help="Only process our debug pixel and its relative chunk")
    parser.add_argument("--minimum_bloom_duration", type=int, default=15, help="the minimum duration in days for a bloom to be retained, this is converted to timesteps by: timesteps_per_year * (days / 365). Default 15.")
    parser.add_argument("--chunk_size", type=int, default=DEFAULT_CHUNKS, help="size of chunks along one side to process (e.g. value of 200 becomes 200x200 chunk")
    parser.add_argument("--minimum_bloom_seperation_duration", type=int, default=10, help="the minimum days between a two initiation dates, this is converted to timesteps by: timesteps_per_year * (days / 365). Default 15.")
    args = parser.parse_args()

    output_folder = args.output_folder if args.output_folder else os.path.dirname(args.chl_location[0])
    zlib_compression = not args.no_compress
    chunk_size = args.chunk_size
    do_only_debug_chunk = args.debug_chunk_only
    time_of_run = datetime.datetime.now().strftime("%Y%m%d%H%M")
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

    fill_chl = not args.no_chl_fill
    med_thresh = 1+ (float(args.median_threshold) / 100)
    do_std_dev = args.std_deviation
    std_dev_threshold = args.std_deviation_threshold
    do_max_means = args.max_means
    max_means_threshold = float(args.max_means_threshold) / 100
    do_rel_max = args.rel_max_means
    do_rel_med_anomaly = args.rel_median
    #remember to change median threshold to percentage!
    #reverse search should be ratio of 100 days so 5 days * 20 = 100 or 8 days * 12.5 (so 13) = 100 days
    #as 100 days is 0.27 of 365 we can just do that
    date_seperation_per_year = int(args.date_seperation_per_year)
    duration_minimum_value = int(round(date_seperation_per_year * (args.minimum_bloom_duration / 365)))
    start_minimum_seperation_value = int(round(date_seperation_per_year * (args.minimum_bloom_seperation_duration / 365)))
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
        chl_lon_var,chl_lat_var,chl_time_var = ds_to_dim_vars(chl_ds)
        REF_MONTH, START_YEAR, missing_chl_dates_at_start, missing_chl_dates_at_end, chl_time_len, date_zero_datetime = get_ds_time_data(chl_ds, chl_time_var, begin_date=args.chl_begin_date, var="chl")
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
        
        logger.info("getting sst missing timesteps")
        if args.sst_location:
            logger.info("sst file provided, reading array")
            logger.info("only one file found, assuming full stack of observations")
            sst_ds = nc.Dataset(args.sst_location)
            sst_lon_var,sst_lat_var,sst_time_var = ds_to_dim_vars(sst_ds)
            sst_ref_month, sst_ref_year, missing_sst_dates_at_start, missing_sst_dates_at_end, sst_time_len, sst_zero_datetime = get_ds_time_data(sst_ds, sst_time_var, begin_date=args.sst_begin_date, var="sst")
            if not (missing_sst_dates_at_start + missing_sst_dates_at_end + sst_time_len) == (missing_chl_dates_at_start + missing_chl_dates_at_end + chl_time_len):
                logger.error("sst time shape : {} chl time shape : {}".format(sst_time_len, chl_time_len))
                logger.error("sst time after padding would not equal chl after padding, got {} for sst and {} for chl, this might be a bug or a problem with your file".format((missing_sst_dates_at_start + missing_sst_dates_at_end + sst_time_len), (missing_chl_dates_at_start + missing_chl_dates_at_end + chl_time_len)))
                sys.exit()
            if not sst_ref_year == START_YEAR:
                logger.error("sst start year and chl start year do not match, got {} for sst and {} for chl".format(sst_ref_year, START_YEAR))
                sys.exit()
        logger.info("getting chunks")

        if args.debug_pixel:
            debug_lat = find_nearest(chl_lats, args.debug_pixel[0])
            debug_lon = find_nearest(chl_lons, args.debug_pixel[1])
            debug_pixel = [debug_lat, debug_lon]
        else:
            debug_pixel = [int(chl_lats.shape[0] * 0.45), int(chl_lons.shape[0] * 0.45)]
        debug_pixel_main = debug_pixel if LON_IDX < LAT_IDX else [debug_pixel[1], debug_pixel[0]]
        logger.info(debug_pixel_main)


        chunks = [(x, y) for x in zip(list(range(0, chl_lons.shape[0], args.chunk_size)),list(range(args.chunk_size,chl_lons.shape[0], args.chunk_size))+ [chl_lons.shape[0]]) 
                         for y in zip(list(range(0, chl_lats.shape[0], args.chunk_size)), list(range(args.chunk_size, chl_lats.shape[0], args.chunk_size)) + [chl_lats.shape[0]])]


        debug_pixel_main = [[debug_pixel_main[0] - chunk[0][0], debug_pixel_main[1] - chunk[1][0], chunk_idx] for chunk_idx, chunk in enumerate(chunks) if chunk[0][0] < debug_pixel_main[0] and chunk[0][1] > debug_pixel_main[0] and chunk[1][0] < debug_pixel_main[1] and chunk[1][1] > debug_pixel_main[1]][0]
        logger.info(debug_pixel_main)
        logger.info("using debug pixel: {} which equates to: {} N {} E (zero indexed so you may need to add 1 to get reference in other software)".format(debug_pixel, chl_lats[debug_pixel_main[0]], chl_lons[debug_pixel_main[1]]))
        
        debug_pixel_only = args.debug_pixel_only


        logger.info("done chunks")
        logger.info("processing will take {} chunks, at 3 minutes a chunk this will take approximately {} minutes".format(len(chunks), len(chunks) * 3))

        sizes = [item for chunk in chunks for item in [chunk[0][1] -chunk[0][0], chunk[1][1] -chunk[1][0]]]
        min_size = min(sizes)
    
        default_logging = logging.INFO
        if not args.stitch_only:
            if not args.no_thread:
                #shut down unhelpful warnings
                logger.info("disabling numpy warnings to reduce spam during threaded processing")
                numpy.warnings.filterwarnings('ignore')
                threads = multiprocessing.cpu_count()
            else:
                threads = 1
            if threads != 1:
                logger.info("threading enabled, revised estimate is {} minutes using the current number of threads ({})".format((len(chunks) * (3 * (args.chunk_size / 100))) /threads, threads))
                logger.info("setting log level to error only")
                logger.setLevel(logging.ERROR)
                default_logging = logging.ERROR
            sst_lock_maj = multiprocessing.Lock()
            chl_lock_maj = multiprocessing.Lock()
            log_info_lock_maj = multiprocessing.Lock()
            debug_chunk = debug_pixel_main[2]
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
            pool.close()
            logger.info("pool closed")
            logger.setLevel(logging.INFO)
            logger.info("chunk processing finished, returned log level to info messages")

            final_output_maxval = os.path.join(output_folder, os.path.basename(chl_filename.replace(".nc", "_phenology_{}_by_maxval.nc".format(time_of_run))))
            final_output_dates = os.path.join(output_folder, os.path.basename(chl_filename.replace(".nc", "_phenology_{}_by_date.nc".format(time_of_run))))
            final_output_intermediate = os.path.join(output_folder, os.path.basename(chl_filename.replace(".nc", "_phenology_{}_intermediate_products.nc".format(time_of_run))))
            logger.info("saving reconstructed files to {}, {}".format(final_output_dates, final_output_maxval))

            create_phenology_netcdf(chl_lons, chl_lats, [1,1], final_output_maxval)
            create_phenology_netcdf(chl_lons, chl_lats, [1,1], final_output_dates)
            create_intermediate_netcdf(final_output_intermediate, chl_lons, chl_lats)
        else:
            logger.info("skipped file creation")
            final_output_maxval = glob.glob(chl_filename.replace(".nc", "_phenology_{}_by_maxval.nc".format(args.time_glob)))[0]
            final_output_dates = glob.glob(chl_filename.replace(".nc", "_phenology_{}_by_date.nc".format(args.time_glob)))[0]
            final_output_intermediate = glob.glob(chl_filename.replace(".nc", "_phenology_{}_intermediate_products.nc".format(args.time_glob)))
            logger.info(final_output_maxval,final_output_dates,final_output_intermediate)

        chl_lons = chl_ds.variables[chl_lon_var][:]
        chl_lats = chl_ds.variables[chl_lat_var][:]

        intermediate_files = glob.glob(chl_filename.replace(".nc","_phenology_{}_chunk*_intermediate*.nc").format(time_of_run))
        files = glob.glob(chl_filename.replace(".nc","_phenology_{}_chunk*_*by*.nc").format(time_of_run))
        logger.info("stitching chunk files")
        try:
            maxval_ds = nc.Dataset(final_output_maxval, 'r+')
            date_ds = nc.Dataset(final_output_dates, 'r+')
            intermediate_ds = nc.Dataset(final_output_intermediate, 'r+')
        except Exception as e:
            logger.error("problem found in files {} {} {}".format(maxval_ds,date_ds,intermediate_ds))
        
        for chunk_idx, chunk in enumerate(tqdm.tqdm(chunks)):
            maxval_chunk_file = [x for x in files if "chunk{}_".format(chunk_idx) in x and 'maxval' in x]
            date_chunk_file = [x for x in files if "chunk{}_".format(chunk_idx) in x and 'date' in x]
            intermediate_chunk_file = [x for x in intermediate_files if "chunk{}_".format(chunk_idx) in x]
            #if we find the tile then add it in
            if len(maxval_chunk_file) and len(date_chunk_file) and len(intermediate_chunk_file):
                logger.debug("processing chunk files: {}, {}".format(os.path.basename(maxval_chunk_file[0]), os.path.basename(date_chunk_file[0])))
                maxval_chunk = nc.Dataset(maxval_chunk_file[0], 'r')
                date_chunk = nc.Dataset(date_chunk_file[0], 'r')
                intermediate_chunk = nc.Dataset(intermediate_chunk_file[0], 'r')
                logger.debug("Files read, preparing to write")
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
                for var in [x for x in maxval_ds.variables.keys() if not x in maxval_ds.dimensions.keys() and x not in ['median_chlorophyll', 'probability', 'chlorophyll_std_dev', 'max_mean', 'perc_val_change']]:
                    try:
                        maxval_ds[var][slc] = maxval_chunk[var][:]
                        date_ds[var][slc] = date_chunk[var][:]
                    except ValueError:
                        logger.warning("Encountered potentially empty file during stitching, this might not be a problem but if you see a large square area of empty data this is likely the cause. Found in chunk {} var {}".format(chunk_idx, var))

                    try:
                        #get the shape of time here
                        time_shape = date_chunk[var].shape[0] + 1 
                    except:
                        continue

                for var in [x for x in maxval_ds.variables.keys() if x in ['median_chlorophyll', 'probability', 'chlorophyll_std_dev', 'max_mean', 'perc_val_change']]:
                    maxval_ds[var][med_prob_slc] = maxval_chunk[var][:]
                    date_ds[var][med_prob_slc] = date_chunk[var][:]

                """
                for var in [x for x in intermediate_ds.variables.keys() if not x in intermediate_ds.dimensions.keys() and x not in ['zen']]:
                    intermediate_ds[var][slc] = intermediate_chunk[var][:]
                """
                logger.debug("Written, closing.")
                maxval_chunk.close()
                date_chunk.close()
                intermediate_chunk.close()
                logger.debug("Removing chunk files.")
                if not args.no_delete:
                    os.remove(maxval_chunk_file[0])
                    os.remove(date_chunk_file[0])

        maxval_ds['TIME'][:] = range(1, time_shape)
        date_ds['TIME'][:] = range(1, time_shape)
        intermediate_ds.close()
        maxval_ds.close()
        date_ds.close()

        logger.info("completed stitching of phenology products, files have been output, now creating intermediate products netcdf. This will take a long time (80 seconds per iteration).")
        try:
            intermediate_ds = nc.Dataset(final_output_intermediate, 'r+')
        except Exception as e:
            logger.error("problem found in files {}".format(intermediate_ds))
            continue
        for chunk_idx, chunk in enumerate(tqdm.tqdm(chunks)):
            logger.info("on chunk")
            logger.info(chunk_idx)
            intermediate_chunk_file = [x for x in intermediate_files if "chunk{}_".format(chunk_idx) in x]
            if len(intermediate_chunk_file):
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
                intermediate_chunk = nc.Dataset(intermediate_chunk_file[0], 'r')
                #only output chl_sbx variabe, not all
                for var in [x for x in intermediate_ds.variables.keys() if not x in intermediate_ds.dimensions.keys() and x not in ['zen']]:
                    logger.info(var)
                    intermediate_ds[var][slc] = intermediate_chunk[var][:]
                    try:
                        time_inter_shape = intermediate_chunk["chl_boxcar"].shape[0]
                    except:
                        continue
                intermediate_chunk.close()
                os.remove(intermediate_chunk_file[0])

        intermediate_ds["TIME"][:] = range(1, time_inter_shape)
        intermediate_ds.close()
        
        logger.info("complete!")


                

