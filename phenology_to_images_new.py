from __future__ import print_function, division
import netCDF4 as nc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import argparse
import os
import calendar
from matplotlib.backends.backend_pdf import PdfPages
import datetime

def output_variables_as_pngs(in_file, reference_date, variable_names=None, specified_year=False):
    """
    Too big of a function but basically just churns through all of our variables and applies the same colour scales to everything. 
    
    Outputs an average where possible.
    """
    data_source = nc.Dataset(in_file)
    lon_var = [x for x in data_source.variables if "lon" in x.lower()][0]
    lat_var = [x for x in data_source.variables if "lat" in x.lower()][0]
    lats = data_source.variables[lat_var]
    lons = data_source.variables[lon_var]
    reverse_search = data_source.getncattr("reverse_search")
    date_zero = datetime.datetime.strptime(data_source.getncattr("date_zero"), "%Y/%m/%d")
    date_seperation_per_year = data_source.getncattr("steps_per_year")
    days_per_step = round(365./date_seperation_per_year)
    central_longitude = lons[-1] - lons[0]
    in_file_bname = in_file.replace(".nc","")
    with PdfPages('{}_output_images_aveonly.pdf'.format(in_file_bname)) as pdf:
        #to months
        #data_source.variables[variable] / (date_seperation / 12)
        #if you just want one or two variables:
        #var_list = ['date_start1', 'date_start2']
        var_list = [v for v in data_source.variables if not v in data_source.dimensions] if not variable_names else variable_names
        for variable in var_list:
            if 'total_blooms' in variable:
                continue
            print(variable) 
            if 'probability' in variable:
                fig1 = plt.figure(figsize=(8, 4), dpi=100)
                m = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
                f1 = plt.pcolormesh(lons, lats, np.ma.masked_invalid(data_source.variables[variable][0]), shading='flat', vmin=0.0, vmax=1.0, cmap=plt.cm.Blues)
                cbar = plt.colorbar(f1, orientation="horizontal", fraction=0.2, pad=0.07, cmap=plt.cm.Blues, ticks=np.arange(0,1,0.1)) 
                cbar.ax.set_xticklabels(['0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'], fontsize=15) 
                #cbar.ax.set_xticklabels(np.arange(0,1,0.1), fontsize=20) 
                cbar.set_label('probability', fontsize=15)
                m.coastlines(resolution='50m', color='black', linewidth=1)
                m.add_feature(cfeature.LAND, facecolor='0.75')
                g1 = m.gridlines(draw_labels = True)
                g1.xlabels_top = False
                g1.xlabel_style = {'size': 2, 'color': 'gray'}
                g1.ylabel_style = {'size': 2, 'color': 'gray'}
                plt.title(variable, fontsize=20)
                #plt.savefig('{}_{}.png'.format(in_file_bname,variable), dpi=100)
                pdf.savefig(fig1)
                plt.close()
            elif 'median_chlorophyll' in variable or "chlorophyll_std_dev" in variable:
                fig1 = plt.figure(figsize=(8, 4), dpi=100)
                m = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
                f1 = plt.pcolormesh(lons, lats, np.ma.masked_invalid(data_source.variables[variable][0]), shading='flat', norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=5), cmap=plt.cm.viridis)
                cbar = plt.colorbar(f1, orientation="horizontal", pad=0.07, ticks=[0.1,0.3,1,5], extend="max", norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=3)) 
                cbar.ax.set_xticklabels(['0.1','0.3','1','3','5'], fontsize=15) 
                cbar.set_label('Chlorophyll, mg m$^{-3}$', fontsize=15)
                m.coastlines(resolution='50m', color='black', linewidth=1)
                m.add_feature(cfeature.LAND, facecolor='0.75')
                g1 = m.gridlines(draw_labels = True)
                g1.xlabels_top = False
                g1.xlabel_style = {'size': 2, 'color': 'gray'}
                g1.ylabel_style = {'size': 2, 'color': 'gray'}
                plt.title(variable, fontsize=20)
                pdf.savefig(fig1)
                plt.close()    
            else:
                #average
                if specified_year:
                    print("year specified, no averaging")
                    plot_vals = data_source.variables[variable][specified_year]
                else:
                    plot_vals = np.ma.mean(data_source.variables[variable][:],axis=0)
                
                fig1 = plt.figure(figsize=(8, 4), dpi=100)
                m = plt.axes(projection=ccrs.PlateCarree(central_longitude=central_longitude))
                if 'max_val' in variable:
                    f1 = plt.pcolormesh(lons, lats, np.ma.masked_invalid(plot_vals[0]), shading='flat', norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=5), cmap=plt.cm.viridis)
                    #if you want numbers:
                    cbar = plt.colorbar(f1, orientation="horizontal", pad=0.07, ticks=[0.1,0.3,1,5], extend="max", norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=3)) 
                    cbar.ax.set_xticklabels(['0.1','0.3','1','3','5'], fontsize=15) 
                    cbar.set_label('Chlorophyll, mg m$^{-3}$', fontsize=15)
                elif 'duration' in variable:
                    f1 = plt.pcolormesh(lons, lats, np.ma.masked_invalid(plot_vals[0]), shading='flat', vmin=1, vmax=30, cmap=plt.cm.gnuplot)
                    cbar = plt.colorbar(f1, orientation="horizontal", pad=0.07, ticks=[1,5,10,15,20,25,30], cmap=plt.cm.gnuplot) 
                    cbar.ax.set_xticklabels(['1','5','10','15','20','25','30'], fontsize=15) 
                    cbar.set_label('Weeks', fontsize=15)    
                else:
                    f1 = plt.pcolormesh(lons, lats, np.ma.masked_invalid(plot_vals[0]), shading='flat', vmin=0 - reverse_search, vmax= date_seperation_per_year + reverse_search, cmap=plt.cm.cubehelix)
                    """
                    cbar = plt.colorbar(f1, orientation="horizontal", pad=0.07, cmap=plt.cm.cubehelix) 
                    cbar.ax.set_xticklabels(range(0 - reference_date, reference_date, 10), fontsize=15) 
                    cbar.set_label('Weeks', fontsize=15)
                    """
                    #uncomment for months as ticks:
                    ticks = range(0 - reverse_search, date_seperation_per_year + reverse_search, int(30 / days_per_step))
                    cbar = plt.colorbar(f1, orientation="horizontal", pad=0.07, cmap=plt.cm.cubehelix,  ticks=ticks) 
                    cbar.ax.set_xticklabels([(date_zero +datetime.timedelta(days=t*days_per_step)).strftime("%b") for t in ticks])
                    cbar.set_label('Months', fontsize=15)
                m.coastlines(resolution='50m', color='black', linewidth=1)
                m.add_feature(cfeature.LAND, facecolor='0.75')
                g1 = m.gridlines(draw_labels = True)
                g1.xlabels_top = False
                g1.xlabel_style = {'size': 2, 'color': 'gray'}
                g1.ylabel_style = {'size': 2, 'color': 'gray'}
                plt.title("{} average".format(variable), fontsize=20)
                #plt.savefig('{}_{}_plot_vals.png'.format(in_file_bname,variable), dpi=100)
                pdf.savefig(fig1)
                plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("in_files", nargs='+', help="A processed phenology file or list of files")
    parser.add_argument("--reference_date", help="e.g. 37, where the 'central' date is")
    parser.add_argument("--variables", default=None,nargs='+', help="list of vars to use")
    parser.add_argument("--specified_year", default=False, help="index of year to plot")
    args = parser.parse_args()

    for f in args.in_files:
        if args.variables:
            output_variables_as_pngs(f, args.reference_date, args.variables, specified_year=args.specified_year)
        else:
            output_variables_as_pngs(f, args.reference_date, specified_year=args.specified_year)
