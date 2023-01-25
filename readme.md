[![DOI](https://zenodo.org/badge/148331696.svg)](https://zenodo.org/badge/latestdoi/148331696)

Requirements

python 2.7 or 3.7

numpy=1.15.4
netCDF4=1.4.2

this can be installed with:

conda create --name phenology_env python=3 anaconda

conda install --name phenology_env numpy=1.15.4 netCDF4=1.4.2

once checked out and installed you need to cd into the directory of this repo and run

`python3 setup.py build_ext --inplace`

this will compile the c module for use with python3, if you are using python 2 use:

`python2 setup.py build_ext --inplace`

Example usage:
`python3 phenology_multibloom.py --sst_location sstder_98_06_global_nofixval_newg.nc chl_5day_MIT_global_98_06.nc --date_seperation_per_year=73`

for use with a pre prepared sst derivitive file and modelled chlorophyll input use:

`python phenology_multibloom.py --sst_location sstder_98_06_global_nofixval.nc fill_chlfillh2_nofixval_OCCCI_NOBM_5d_98_06_newg.nc --date_seperation_per_year 73 --skip_sst_prep --modelled_median --first_date_index 37`

This will produce a file named chl_5day_MIT_global_98_06_phenology.nc, using the two bloom method with time slices matched to every 73 entries (5 day composites).

For multiple files usage is as follows:
`python3 phenology_multibloom.py --sst_location sstder_98_06_global_nofixval_newg.nc chl_5day_MIT_global_98_06.nc chl_5day_MIT_global_98_06_2.nc --date_seperation_per_year=73`

Which creates file named chl_5day_MIT_global_98_06_phenology.nc and a file named chl_5day_MIT_global_98_06_2_phenology.nc

```
usage: phenology_multibloom.py [-h] [--sst_location SST_LOCATION]
                               [--output OUTPUT]
                               [--date_seperation_per_year DATE_SEPERATION_PER_YEAR]
                               [--first_date_index FIRST_DATE_INDEX]
                               [--intermediate_file_store INTERMEDIATE_FILE_STORE]
                               [--median_threshold MEDIAN_THRESHOLD]
                               chl_location

positional arguments:
  chl_location          A chlorophyll file or glob for files (e.g. chl*.nc)
                        that can be stacked into an array of shape
                        (time,depth,lat,lon). These must have a variable
                        called chl, lat and lon (or contain those strings)
                        which will be used to define the array shape

optional arguments:
  -h, --help            show this help message and exit
  --sst_location SST_LOCATION
                        An sst file, or glob of files (e.g. sst*.nc) that
                        matches the chlorophyll observations. if it does not
                        match some interpolation can be attempted with
                        --sst_date_interp
  --output OUTPUT       output location, default is phenology.nc in current
                        folder
  --date_seperation_per_year DATE_SEPERATION_PER_YEAR
                        how many temporal observations we have in a year, if
                        not specified will be guessed
  --first_date_index FIRST_DATE_INDEX
                        specify if the first date you want to include is not
                        the first date present in the date stack
  --intermediate_file_store INTERMEDIATE_FILE_STORE
                        change where intermediate numpy files are placed, if
                        not specified then /tmp is assumed - you should
                        specify somewhere else if your tmp cannot handle the
                        array sizes needed (and currently this program will
                        fill it until it cannot).
  --median_threshold MEDIAN_THRESHOLD
                        change median threshold
```
