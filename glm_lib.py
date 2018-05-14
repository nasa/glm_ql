#!/usr/bin/env python
"""
------------
glm_lib.py

Purpose: function library

------------
Usage:

Not standalone. For use with Geostationary Lightning Mapper Quick Look
Software

------------
Author: Matt Smith (Univ. of Alabama in Huntsville - Information
                    Technology and Systems Center)
                   (NASA SPoRT - MSFC - Huntsville, Alabama, USA)
        matthew.r.smith@nasa.gov

Date: November 2017

Version: 1.0

------------
Modifications:

------------
Description:

This script contains functions used within glm_cspp.py.

------------
Notes:

------------
Required Files:
glm_config.py

------------
License Info:

Geostationary Lightning Mapper Quick Look Software,
United States Government as represented by the Administrator of the
National Aeronautics and Space Administration. All rights reserved.

Geostationary Lightning Mapper Quick Look Software is licensed under the
NOSA, Version 1.3 (the "License"); you may not use this file except in
compliance with the License. You may obtain a copy of the License at:

* * * * https://github.com/nasa/glm_ql/blob/master/GLM_Quicklook_NOSA.txt

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

# Standard Libraries
import logging
import shutil
from datetime import datetime
import warnings
import os
import calendar
import sys

# Third Party Libraries
import matplotlib       # to avoid 'failed to get the current screen resources'
matplotlib.use('Agg')   # to avoid 'failed to get the current screen resources'
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
import numpy as np

# Local Libraries
import glm_config

# Conditionally import GDAL-related modules (since its import is sometimes
#  problematic)
GDAL_IMPORTED = True
try:
    import gdalconst
except ImportError:
    logging.warning('No GDAL found')
    GDAL_IMPORTED = False

try:
    from osgeo import osr, gdal
except ImportError:
    logging.warning('No osgeo found')
    GDAL_IMPORTED = False

# Conditionally import Cartopy-related modules (required for PNG creation)
try:
    import cartopy.crs as ccrs
except ImportError:
    logging.critical('No cartopy found')
    sys.exit(-1)

try:
    import cartopy.feature as cfeature
except ImportError:
    logging.critical('No cartopy.feature found')
    sys.exit(-1)

# Filter UserWarnings from matplotlib's LONGITUDE_FORMATTER/LATITUDE_FORMATTER
warnings.filterwarnings('ignore', ".*Steps argument.*", UserWarning)
try:
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
except ImportError:
    logging.critical('No cartopy.mpl.gridliner found')
    sys.exit(-1)

####################################################
def copyFileToLocalDir(fil, local_dir):
    """Copy file (fil) to local directory"""

    logging.info('Copying %s to %s', fil, local_dir)
    try:
        shutil.copy(fil, local_dir)
    except:
        logging.warning('Error copying file')
    else:
        logging.info('File copied')

####################################################
def setup_domain(args):
    """Setup domain grid for binning features."""

    global MAX_BINS_Y, MAX_BINS_X
    global PIX_PER_DEG_X, PIX_PER_DEG_Y, SIZ_LAT, SIZ_LON, DEG_PER_PIX_X, DEG_PER_PIX_Y

    logging.info('Start setup_domain')
    SIZ_LAT = glm_config.DOMAIN_OPTS[args.domain]['URlat'] \
            - glm_config.DOMAIN_OPTS[args.domain]['LLlat']
    SIZ_LON = glm_config.DOMAIN_OPTS[args.domain]['URlon'] \
            - glm_config.DOMAIN_OPTS[args.domain]['LLlon']
    DEG_PER_PIX_X = DEG_PER_PIX_Y = glm_config.DOMAIN_OPTS[args.domain]['resKm'] / 111.1
    PIX_PER_DEG_X = 1. / DEG_PER_PIX_X
    PIX_PER_DEG_Y = 1. / DEG_PER_PIX_Y
    MAX_BINS_Y = int((glm_config.DOMAIN_OPTS[args.domain]['URlat'] \
               - glm_config.DOMAIN_OPTS[args.domain]['LLlat']) * PIX_PER_DEG_Y)
    MAX_BINS_X = int((glm_config.DOMAIN_OPTS[args.domain]['URlon'] \
               - glm_config.DOMAIN_OPTS[args.domain]['LLlon']) * PIX_PER_DEG_X)

####################################################
def array_2_gtif(tif_array, tif_file, args):
    """
    Generate GeoTIFF using data in tif_array using GDAL.
    Copy to GTIF_DIR if directed.
    """

    logging.info('Start array_2_gtif')

    # Set up output to geoTIFF
    driver = gdal.GetDriverByName('GTiff')
    out_raster = driver.Create(tif_file, MAX_BINS_X, MAX_BINS_Y, 1,
                               gdal.GDT_Byte)

    # pixelHeight negative since we have upper left origin
    out_raster.SetGeoTransform((glm_config.DOMAIN_OPTS[args.domain]['LLlon'],
                                DEG_PER_PIX_X,
                                0,
                                glm_config.DOMAIN_OPTS[args.domain]['URlat'],
                                0,
                                -DEG_PER_PIX_Y))
    out_band = out_raster.GetRasterBand(1)
    out_band.WriteArray(tif_array)

    # Get the spatial reference information
    out_raster_srs = osr.SpatialReference()
    out_raster_srs.ImportFromEPSG(4326)
    out_raster.SetProjection(out_raster_srs.ExportToWkt())

    out_band.FlushCache()

    # Make sure the file is explicitly closed
    out_raster = None
    del out_band, out_raster

    # If sendGtifTo flag set, copy geoTIFF to GTIF_DIR.
    if glm_config.DOMAIN_OPTS[args.domain]['sendGtifTo']:
        copyFileToLocalDir(tif_file, glm_config.GTIF_DIR)

####################################################
def array_2_png(array, png_file, title, args):
    """
    Generate PNG image using array.
    This is where cartopy is sometimes hard to manage. Whitespace is rampant
    at times. Copy to PNG_DIR if directed.
    """

    logging.info('Start array_2_png')
    # Plot it and save it to a PNG
    # Determine latitude labeling increment
    lat_inc = 10.
    if SIZ_LAT <= 8:
        lat_inc = 1.
    elif SIZ_LAT <= 16:
        lat_inc = 2.
    elif SIZ_LAT <= 25:
        lat_inc = 5.

    # Determine longitude labeling increment
    lon_inc = 10.
    if SIZ_LON <= 8:
        lon_inc = 1.
    elif SIZ_LON <= 16:
        lon_inc = 2.
    elif SIZ_LON <= 25:
        lon_inc = 5.

    # Extent in lat/lon coords (left, right, bottom, top)
    geo_extent = [
        glm_config.DOMAIN_OPTS[args.domain]['LLlon'],
        glm_config.DOMAIN_OPTS[args.domain]['URlon'],
        glm_config.DOMAIN_OPTS[args.domain]['URlat'],
        glm_config.DOMAIN_OPTS[args.domain]['LLlat'],
        ]

    # Set up plot
    fig, ax = plt.subplots(1, 1,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(geo_extent, crs=ccrs.PlateCarree())
    ax.background_patch.set_facecolor('black')

    # Draw Coastlines, Countries, States, and Lakes
    coastlines = cfeature.NaturalEarthFeature(
        'physical',
        'coastline',
        glm_config.DOMAIN_OPTS[args.domain]['mapres'],
        edgecolor='0.8',
        facecolor='none'
        )
    ax.add_feature(coastlines, linewidth=0.2,)

    countries = cfeature.NaturalEarthFeature(
        'cultural',
        'admin_0_countries',
        glm_config.DOMAIN_OPTS[args.domain]['mapres'],
        edgecolor='0.5',
        facecolor='none'
        )
    ax.add_feature(countries, linewidth=0.2,)

    states = cfeature.NaturalEarthFeature(
        'cultural',
        'admin_1_states_provinces',
        glm_config.DOMAIN_OPTS[args.domain]['mapres'],
        edgecolor='0.4',
        facecolor='none'
        )
    ax.add_feature(states, linewidth=0.1,)

    lakes = cfeature.NaturalEarthFeature(
        'physical',
        'lakes',
        glm_config.DOMAIN_OPTS[args.domain]['mapres'],
        edgecolor='0.4',
        facecolor='none'
        )
    ax.add_feature(lakes, linewidth=0.1,)

    # Setup Gridlines (Lats and Lons)
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.4,
        color='0.6',
        alpha=0.3,
        linestyle='dashed')
    gl.xlabels_top = False
    lon = int(glm_config.DOMAIN_OPTS[args.domain]['LLlon']-lon_inc)
    lon_list = [lon]
    while lon <= glm_config.DOMAIN_OPTS[args.domain]['URlon']:
        lon += lon_inc
        lon_list.append(lon)
    gl.xlocator = mticker.FixedLocator(lon_list)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'color':'black', 'size':3}

    gl.ylabels_right = False
    lat = int(glm_config.DOMAIN_OPTS[args.domain]['LLlat']-lat_inc)
    lat_list = [lat]
    while lat <= glm_config.DOMAIN_OPTS[args.domain]['URlat']:
        lat += lat_inc
        lat_list.append(lat)
    gl.ylocator = mticker.FixedLocator(lat_list)
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'color':'black', 'size':3}

    # Plot the data!
    img = ax.imshow(
        array,
        extent=geo_extent,
        aspect=1,
        alpha=1.0,
        cmap=glm_config.GLM_COLORS,
        vmin=0,
        vmax=glm_config.NUM_COLORS
        )

    # Make adjustments to locations of tick mark labels
    tick_offsets = [0.0]
    for tck in range(0, glm_config.NUM_COLORS + 1):
        tick_offsets.append(0.5 + tck - (tck / 12.))

    # Specs: [left, bottom, width, height]
    cbar_ax = fig.add_axes([0.15, 0.06, 0.70, 0.015])
    clrbar = fig.colorbar(img, cax=cbar_ax, orientation='horizontal',
                          ticks=tick_offsets)
    cbar_ax.tick_params(labelsize=3, length=-0.1)
    cbar_ax.set(adjustable='datalim')

    # Generate list of colorbar labels
    # Legend is necessarily a part of the colorbar
    if args.feature == 'flash':
        # Generate 'X through Y-1' strings
        labels = ['# Flashes', '0']
        for val in range(0, glm_config.NUM_COLORS - 1):
            labels.append(str(flash_bins_temp[val]) + '-'
                          + str(flash_bins_temp[val + 1] - 1))
        # Add last value ('>N')
        labels.append('>' + str(flash_bins_temp[glm_config.NUM_COLORS -1]))

    elif args.feature == 'group':
        # Generate 'X through Y-1' strings
        labels = ['# Groups', '0']
        for val in range(0, glm_config.NUM_COLORS -1):
            labels.append(str(group_bins_temp[val]) + '-'
                          + str(group_bins_temp[val + 1] - 1))
        # Add last value ('>N')
        labels.append('>' + str(group_bins_temp[glm_config.NUM_COLORS - 1]))

    # Plot colorbar ticks and title
    clrbar.ax.set_xticklabels(labels)
    clrbar.ax.sep = 1
    clrbar.set_label(title, size=7, labelpad=3)

    # Save to file for web graphics
    # specify dpi - otherwise it defaults to 100 dpi
    plt.savefig(png_file, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close()

    # If sendPngTo flag set, copy PNG to PNG_DIR.
    if glm_config.DOMAIN_OPTS[args.domain]['sendPngTo']:
        copyFileToLocalDir(png_file, glm_config.PNG_DIR)

####################################################
def process_minute(file_list, args):
    """
    Process the specified list of files. Generate PNG and optional GeoTIFF.
    """

    logging.info('Start process_minute')

    # Get and calculate all domain-related variables
    setup_domain(args)
    lit_array = np.zeros((MAX_BINS_Y, MAX_BINS_X), dtype='int')
    glm_file_list = []  # start file list (without paths)

    first_file = True
    for abs_file in file_list:
        glm_file = os.path.basename(abs_file)
        glm_file_list.append(glm_file)

        if first_file:
            # Get start date & time of the first file (hh:mm:00)
            syyyy = glm_file[20:24]
            sddd = glm_file[24:27]
            shh = glm_file[27:29]
            smm = glm_file[29:31]
            sss = glm_file[31:33]
            smonthday = datetime.strptime(syyyy+sddd, '%Y%j').strftime('%m%d')
            smon = smonthday[:2]
            sdom = smonthday[2:]
            first_file = False

        logging.debug('glm_file = %s', glm_file)

        # Open GLM data file to read
        f_h = Dataset(abs_file, 'r')

        # Get number of Flashes and Groups
        num_flashes = len(f_h.dimensions['number_of_flashes'])
        num_groups = len(f_h.dimensions['number_of_groups'])

        # Read group data
        gpf_ids = f_h.variables['group_parent_flash_id'][:]

        # Process groups or flashes (not events yet)
        if args.feature == 'group':

            # Read lats and lons
            glons = f_h.variables['group_lon'][:]
            glats = f_h.variables['group_lat'][:]

            # Loop through all groups for *binning*
            for group in range(num_groups):

                # Skip group if its parent flash has fewer than 3 groups
                if (gpf_ids == gpf_ids[group]).sum() < 3:
                    continue

                # Determine x and y bin values
                y = int((glats[group] - \
                        glm_config.DOMAIN_OPTS[args.domain]['LLlat']) * PIX_PER_DEG_Y)
                x = int((glons[group] - \
                        glm_config.DOMAIN_OPTS[args.domain]['LLlon']) * PIX_PER_DEG_X)

                # Skip this group if off the grid
                if y < 0 or x < 0 or y >= MAX_BINS_Y or x >= MAX_BINS_X:
                    continue

                # Increment counter
                lit_array[y][x] += 1

        else: # Flashes

            flons = f_h.variables['flash_lon'][:]
            flats = f_h.variables['flash_lat'][:]

            # Loop through all flashes for *binning*
            for flash in range(num_flashes):

                # Determine x and y bin values
                y = int((flats[flash] - \
                        glm_config.DOMAIN_OPTS[args.domain]['LLlat']) * PIX_PER_DEG_Y)
                x = int((flons[flash] - \
                        glm_config.DOMAIN_OPTS[args.domain]['LLlon']) * PIX_PER_DEG_X)

                # Skip this flash if off the grid
                if y < 0 or x < 0 or y >= MAX_BINS_Y or x >= MAX_BINS_X:
                    continue

                # Increment counter
                lit_array[y][x] += 1

        f_h.close()

    # Adjust all bin values for >8km resolution
    global flash_bins_temp, group_bins_temp
    flash_bins_temp = list(glm_config.FLASH_BINS)
    group_bins_temp = list(glm_config.GROUP_BINS)
    if glm_config.DOMAIN_OPTS[args.domain]['resKm'] != 8:
        flash_bins_temp = [1]
        group_bins_temp = [1]
        for val in range(1, glm_config.NUM_COLORS):
            flash_bins_temp.append(
                int(round(glm_config.FLASH_BINS[val] \
                    * (glm_config.DOMAIN_OPTS[args.domain]['resKm']**2) / 64.)))
            group_bins_temp.append(
                int(round(glm_config.GROUP_BINS[val] \
                    * (glm_config.DOMAIN_OPTS[args.domain]['resKm']**2) / 64.)))

    # Re-Bin lit_array with digitized binnings
    if args.feature == 'group':
        lit_array_new = np.digitize(lit_array.flatten(), group_bins_temp)
    elif args.feature == 'flash':
        lit_array_new = np.digitize(lit_array.flatten(), flash_bins_temp)

    png_array = lit_array_new.reshape(MAX_BINS_Y, MAX_BINS_X).astype('int8')

    logging.debug('max of png_array=%d', np.max(png_array))

    png_file = os.path.join(glm_config.PROC_DIR, 'goes16_glm_' + args.domain
                            + '_' + args.feature + '_' + syyyy + smon + sdom
                            + '_' + shh + smm + '.png')
    logging.debug('png_file=%s', png_file)
    if args.feature == 'group':
        label = 'GOES-16 GLM Group Density ' + syyyy + ' ' \
                + calendar.month_abbr[int(smon)] + ' ' + sdom + ' ' + shh \
		+ ':' + smm + ' UTC (1 min)\nPreliminary Non-Operational'
    elif args.feature == 'flash':
        label = 'GOES-16 GLM Flash Centroid Density ' + syyyy + ' ' \
                + calendar.month_abbr[int(smon)] + ' ' + sdom + ' ' + shh \
                + ':' + smm + ' UTC (1 min)\nPreliminary Non-Operational'
    array_2_png(np.flipud(png_array), png_file, label, args)

    if args.g:
        tif_file = os.path.join(glm_config.PROC_DIR, 'goes16_glm_' \
                   + args.domain + '_' + args.feature + '_' + syyyy + smon \
                   + sdom + '_' + shh + smm + '.tif')
        logging.debug('tif_file=%s', tif_file)
        tif_array = np.flipud(lit_array_new.reshape(MAX_BINS_Y, MAX_BINS_X).\
                              astype('int8'))
        logging.debug('max of tif_array=%d', np.max(tif_array))
        array_2_gtif(tif_array, tif_file, args)

    # Cleanup the (current) processing directory
    os.remove(png_file)

####################################################
def get_args(sys_args):
    """
    Parse arguments supplied with command. Available FEATURES and DOMAINS
    are specified in glm_config.py.
    """

    import argparse

    log = logging.getLogger(__file__)
    logging.debug('Start get_args')

    parser = argparse.ArgumentParser(
        formatter_class=lambda prog: argparse.HelpFormatter(prog,
                                                            max_help_position=30),
        prog=sys.argv[0],
        usage='%(prog)s [options]   Generate GLM quicklook\
 for 1 minute of GLM files (3)'
        )
    parser.add_argument(
        'glmDir',
        help='local directory where GLM "OR" files are located'
        )
    parser.add_argument(
        'feature',
        help='lightning feature to display\n',
        choices=glm_config.FEATURES
        )
    parser.add_argument(
        'domain',
        help='domain to display\n',
        choices=glm_config.DOMAINS
        )
    if GDAL_IMPORTED:
        parser.add_argument(
            '-g',
            action='store_true',
            help='generate geoTIFF'
            )
    parser.add_argument(
        '-time',
        type=int,
        help='yyyyMMddhhmm or yyyydddhhmm',
        metavar="[2017-2025]"
        )
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='diagnostic output'
        )
    parser.add_argument(
        '-r',
        action='store_true',
        help='reverse process (oldest 1st) directory of files'
        )

    args = parser.parse_args()

   # Check validity of supplied directory arg.
    if not os.path.isdir(args.glmDir):
        logging.critical('%s is not a valid directory', args.glmDir)
        sys.exit(-1)

    # set to DEBUG - otherwise defaults to WARNING
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    return args

####################################################
def check_time(t_0):
    """Perform sanity check on specific supplied date/time."""

    days_in_month = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, \
                   7:31, 8:31, 9:30, 10:31, 11:30, 12:31}
    days_in_year = 365
    if calendar.isleap(int(str(t_0)[:4])):
        days_in_month[2] = 29
        days_in_year = 366

    # dayOfYear designation
    if t_0 < 201701010000:
        # Check date range for dayOfYear formatted arg.
        if t_0 < 20170010000 or t_0 > 20253652359:
            logging.critical('Time (year) out of range 2017-2025')
            sys.exit(-1)
        ddd = int(str(t_0)[4:7])
        if ddd > days_in_year:
            logging.critical('Time (day-of-year) out of range 001-%03d',
                             days_in_year)
            sys.exit(-1)
        return datetime.strptime(str(t_0), '%Y%j%H%M')

    # monthDay designation
    else:
        # Check date range monthDay formatted arg.
        if t_0 < 201701010000 or t_0 > 202512312359:
            logging.critical('Time (year) out of range 2017-2025')
            sys.exit(-1)
        month = int(str(t_0)[4:6])
        if month > 12:
            logging.critical('Time (month) out of range 01-12')
            sys.exit(-1)
        day = int(str(t_0)[6:8])
        if day > days_in_month[month]:
            logging.critical('Time (day) out of range 01-%02d for month %02d',
                             days_in_month[month], month)
            sys.exit(-1)
        hour = int(str(t_0)[8:10])
        if hour > 23:
            logging.critical('Time (hour) out of range 00-23')
            sys.exit(-1)
        minute = int(str(t_0)[10:12])
        if minute > 59:
            logging.critical('Time (minute) out of range 00-59')
            sys.exit(-1)
        return datetime.strptime(str(t_0)[:12], '%Y%m%d%H%M')
