#!/usr/bin/env python
"""
------------
glm_cspp.py

------------
Purpose:

Create GLM quick look images and optional geoTIFFs

------------
Usage:

glm_cspp.py glmDir feature domain [-t yyyyMMddhhmm|yyyydddhhmm] [-g][-d][-r]

where:
    glmDir:      local directory where 20-sec GLM 'OR' files are located
    feature:     flash | group (or 'event' someday)
    domain:      fulldisk|southamer|conus|carib|seus|swus|neus|nwus
    [-t minute]: option - date/time (month/day or day-of-year format)
    [-g]:        option - also generate geoTIFF
    [-r]:        option - process in reverse order (oldest first)
    [-d]:        option - debug mode

------------
Author:

Matt Smith (Univ. of Alabama in Huntsville - Information Technology and
            Systems Center)
           (NASA SPoRT - MSFC)
    Huntsville, Alabama, USA)

    matthew.r.smith@nasa.gov

Date: November 2017

Version: 1.0

------------
Notes:

This code will process every minute of GLM data (3 20-sec files) in 'glmDir',
generating PNGs/GeoTIFFs using the specified domain, for the specified
lighting feature.

* Constants should be set in glm_config.py.

If an optional 'minute' is specified, only that minute's PNG/GeoTIFF will be
generated.

If any of the 3 files (:00, :20, :40) from any minute are missing - a
PNG/GeoTIFF will NOT be generated. The first are referred to as 00-sec.

Events are not handled yet. As of February 2018, the geolocation for events is
specious.

The origin of the data array is the Upper Left corner. The data array first
dimension is y (lat).

Domains with resolutions (resKm) higher than 8km are binned differently -
to account for the larger gridbox area (in function process_minute). I.e., the
digitizations are intended for 8km (64 sq. km) boxes. The flash_bin and group_bin
lists are modified if not an 8km box.

This code was written for python3 and using cartopy as the cartographic
package. The bulk of the import issues with GDAL and cartopy are located in
glm_lib.py. Cartopy has occasional trouble with unwanted whitespace.

GLM Level-2 file names (example):
OR_GLM-L2-LCFA_G16_s20171170850000_e20171170850200_c20171170850224.nc

OR:     Operational System - Realtime Data
GLM:    Geostationary Lightning Mapper (on GOES-16 and -17)
L2:     Level 2
LCFA:   Lightning Cluster Filter Algorithm
G16|17: GOES-16 and -17
s...:   start date/time (yyyydddhhmmsst) yr, dayOfYr, hr, min, sec, tenthOfSec
e...:   end date/time
c...:   creation date/time
nc:     NetCDF(4)

------------
Required files:

-glm_cspp.py (this file)
    main
-glm_config.py (CONSTANTS)
-glm_lib.py (Functions)
    copyFileToLocalDir
    setup_domain
    array_2_gtif
    array_2_png
    process_minute
    get_args
    check_time

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
import glob
from datetime import datetime
import os
import sys

# Third Party Libraries

# Local Libraries
import glm_lib
import glm_config

# Initiate logging, specifying log report format
LOG = logging.getLogger(__file__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', \
    datefmt='%Y/%m/%d %H:%M:%S')
##################################
#
#  MAIN PROGRAM
#
##################################
def main():
    """
    Main program
    """

    # Parse arguments
    args = glm_lib.get_args(sys.argv)

    logging.info('Start main')

    # Generate PNGs for all minutes in directory glmDir.
    #
    if args.time is None:

        # Make list of all files in 'glmDir' that have :00 seconds start time.
        glm_list_00 = glob.glob(os.path.join(
            args.glmDir,
            glm_config.GLM_FILE_PREFIX + '20?????????00?_e*_c*.nc'
            ))
        if not glm_list_00: # check if no files found
            logging.error('no 00-sec GLM files found')
            sys.exit(-1)
        logging.info('Found %d 00-sec GLM files in %s', len(glm_list_00), args.glmDir)

        # Loop through all :00 GLM files (most recent first)
        for glm_file_00 in sorted(glm_list_00, reverse=(not args.r)):

            # Begin list of 3 files for this minute
            minute_list = [glm_file_00]

            # string for this minute
            this_min = os.path.basename(glm_file_00)[20:31]

            # Check the existence of the xx:20 file
            glm_file_20 = glob.glob(os.path.join(
                args.glmDir, glm_config.GLM_FILE_PREFIX + this_min + '200*'))

            # Check the existence of the xx:40 file
            glm_file_40 = glob.glob(os.path.join(
                args.glmDir, glm_config.GLM_FILE_PREFIX + this_min + '400*'))

            # Verify one (and only one) of each file
            if len(glm_file_20) == 1 and len(glm_file_40) == 1:
                minute_list.append(glm_file_20[0])
                minute_list.append(glm_file_40[0])

            # Otherwise, complain and skip it
            else:
                logging.warning('Cannot generate quicklook for %s', this_min)
                continue

            logging.info('\nProcessing %s', this_min)

            # Generate PNG (and GeoTIFF) for this minute
            glm_lib.process_minute(minute_list, args)

    # Generate PNGs for specific minute in directory glmDir.
    #
    else:

        # Perform sanity check of specified time
        d_t = glm_lib.check_time(args.time)

        # Get list of files for this minute
        glm_list = glob.glob(os.path.join(
            args.glmDir, glm_config.GLM_FILE_PREFIX
            + d_t.strftime('%Y%j%H%M') # Converts to day-of-year
            + '?0?_e*_c*.nc'
            ))

        if len(glm_list) != 3:
            logging.critical('Could not find 3 files for minute: %s',\
                  datetime.strftime(d_t, '%Y%m%d%H%M'))
            sys.exit(-1)

        logging.info('Processing %s', d_t.strftime('%Y%m%d %H:%M'))

        # Generate PNG (and GeoTIFF) for this minute
        glm_lib.process_minute(glm_list, args)

    logging.info('End main')

######################################
#if being called as main.py - execute
if __name__ == "__main__":
    main()
