import xarray
import pandas
import numpy
import datetime
import dask

netcdf_file = "/local0/sst_reshaped_new.nc"

output = "/local0/sst_reshaped_new_5d_comp.nc"

input_file = xarray.open_dataset(netcdf_file)
years = set([int(f.dt.year) for f in input_file.sst.TIME])
end_list=[]
ods = None
for year in years:
    print(year)
    yr_sst = input_file.sst.sel(TIME=input_file.sst.TIME.dt.year.isin([year-1, year, year+1]))
    if len(input_file.sst.TIME.dt.year.isin([year])) != 12:
        print([year-1, year, year+1])
        print(input_file.sst.TIME.dt.year.isin([year-1, year, year+1]))
        print("using min max")
        min_month = min([int(x.dt.month) for x in yr_sst.TIME]) 
        max_month = max([int(x.dt.month) for x in yr_sst.TIME]) 
        sst = yr_sst.interp(TIME=pandas.date_range(start=datetime.datetime(year=year, month=min_month, day=1), end= datetime.datetime(year=year, month=max_month, day=31), freq='5D'))  
    else:
        sst = yr_sst.interp(TIME=pandas.date_range(start=datetime.datetime(year=year, month=1, day=1), end= datetime.datetime(year=year, month=12, day=31), freq='5D'))
    sst.to_netcdf(path="preproc_interp_to_delete_5d_{}.nc".format(year))

sst=None

joined_arrays = xarray.open_mfdataset('preproc_interp_to_delete_5d_*.nc', concat_dim='TIME')
joined_arrays.sst.to_netcdf(output, encoding={'sst':{'zlib': True, 'complevel': 9}})
