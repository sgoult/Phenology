
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os

def find_nearest(array, value):
    array = numpy.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
    return idx


def gen_plots(phen_file, lat, lon, start, stop, pdf, xsize=35, ysize=3, start_index=1, date_seperation_per_year=46, intermediate_file=None):
    t = nc.Dataset(phen_file.replace("_by_date","_intermediate_products").replace("_by_maxval","_intermediate_products"))
    timings = nc.Dataset(phen_file)
    lon_var = [x for x in t.variables if "lon" in x.lower()][0]
    lat_var = [x for x in t.variables if "lat" in x.lower()][0]
    lats = t.variables[lat_var][:]
    lons = t.variables[lon_var][:]
    print "lat min max:", lats.min(), lats.max()
    print "lon min max:", lons.min(), lons.max()
    lat_idx = find_nearest(lats, lat)
    lon_idx = find_nearest(lons, lon)
    ser =  t.variables["filled_chl"][start:stop,:,lat_idx,lon_idx]
    x = range(start, stop)
    fig1 = plt.figure(figsize=(xsize,ysize), dpi=300)
    plt.title("{}, {}".format(lat, lon))
    plt.plot(x, ser, color='blue', label="filled chl")
    ser =  t.variables["chl_boxcar"][start:stop,:,lat_idx,lon_idx]
    x = range(start, stop)
    plt.plot(x, ser, color='red', label="smoothed filled chl")
    #add vertical lines for start (red) and end (green) of each year
    last_start = None
    last_end = None
    last_date_start = None
    last_date_start_two = None
    for year_indx, year in enumerate(range(start_index, stop, date_seperation_per_year)):
        if not (year >= start and year < stop):
            continue
        last_start = plt.axvline(x=year, color='red')
        last_end = plt.axvline(x=year +  date_seperation_per_year -1, color='green')
        date_start =  timings.variables["date_start1"][year_indx,:,lat_idx,lon_idx]
        date_max =  timings.variables["date_max1"][year_indx,:,lat_idx,lon_idx]
        date_end =  timings.variables["date_end1"][year_indx,:,lat_idx,lon_idx]
        duration =  timings.variables["duration1"][year_indx,:,lat_idx,lon_idx]
        if not numpy.ma.is_masked(date_start) and not (date_start+year >stop and date_stop < start):
            last_date_start = plt.axvline(x=year + date_start, color='yellow')
            plt.axvline(x=year + date_max - 1, color='yellow')
            print "primary"
            print date_start, date_max, date_end
            plt.axvline(x=year + date_end, color='yellow')
        date_start =  timings.variables["date_start2"][year_indx,:,lat_idx,lon_idx]
        date_max =  timings.variables["date_max2"][year_indx,:,lat_idx,lon_idx]
        date_end =  timings.variables["date_end2"][year_indx,:,lat_idx,lon_idx]
        duration =  timings.variables["duration2"][year_indx,:,lat_idx,lon_idx]
        if not numpy.ma.is_masked(date_start)  and not (date_start+year >stop and date_stop < start):
            last_date_start_two = plt.axvline(x=year + date_start, color='purple')
            plt.axvline(x=year + date_max -1, color='purple')
            print "secondary"
            print date_start, date_max, date_end
            plt.axvline(x=year + date_end, color='purple')

    if last_start:
        last_start.set_label('year start')
    if last_end:
        last_end.set_label('year end')
    if last_date_start:
        last_date_start.set_label('date_start1 metrics')
    if last_date_start_two:
        last_date_start_two.set_label('date_start2 metrics')
    plt.legend()
    pdf.savefig(fig1)

def main(args):
    points = [p.split(",") for p in args.points]
    if not args.output_folder:
        output_file = args.phen_file.replace("nc", "pdf")
    else:
        output_file = os.path.join(args.output_folder, 
                                   os.path.basename(args.phen_file.replace("nc", "pdf")))

    with PdfPages(output_file) as pdf:
        for point in points:
            print point
            gen_plots(args.phen_file, float(point[0]), float(point[1]), 0, 1012, pdf)
            for year in range(0, 1012, 46):
                print year
                gen_plots(args.phen_file, float(point[0]), float(point[1]), year, year+46, pdf, xsize=6, ysize=6)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("phen_file", help="A processed phenology file")
    parser.add_argument("--output_folder", default=None, help="a location to put things, if not specified default to chlorophyll file location")
    parser.add_argument("points", nargs="+", help="A point or list of points in format lat,lon lat,lon etc. E.g. 70.93,-13.11")
    parser.add_argument("--intermediate_file", default=None, help="associated intermediate file, if not given will default to replacing by_data/by_maxval with intermediate_products")
    args = parser.parse_args()
    main(args)