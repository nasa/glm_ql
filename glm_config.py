#!/usr/bin/env python
"""
------------
glm_config.py

Purpose: to designate constants

------------
Usage:

Not standalone. Only for use with Geostationary Lightning Mapper Quick
Look Software

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

This script sets several constants, such as directory locations, a color
table, and domain extents.

------------
Notes:

------------
Required Files:

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
import os

# Third Party Libraries
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# Local Libraries

# GLM lightning features allowed (event lat/lon not yet stable - 2/20/2018)
FEATURES = ['flash', 'group']

# Domains allowed
DOMAINS = ['fulldisk', 'southamer', 'conus']

# Set processing directory
# (set to desired location)
PROC_DIR = '/data1/msmith/glm/processing'

# Set PNG destination directory
# (set to desired location)
PNG_DIR = os.path.join(os.environ['HOME'], 'glm/png')

# Set GeoTIFF destination directory
# (set to desired location)
GTIF_DIR = os.path.join(os.environ['HOME'], 'glm/gtiff')

# Operational GLM filename prefix
GLM_FILE_PREFIX = 'OR_GLM-L2-LCFA_G16_s'

# Digitizations for each lightning feature (based on 8km resolution binning)
GROUP_BINS = np.array([1, 3, 10, 20, 40, 70, 100, 150, 200, 300, 500])
# color #             0 1  2   3   4   5   6    7    8    9   10   11

FLASH_BINS = np.array([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
# color #             0 1   2   3   4   5   6   7   8   9  10   11

# Note - we are assuming the length of 'FLASH_BINS' is the same as 'GROUP_BINS'
if len(FLASH_BINS) == len(GROUP_BINS):
    NUM_COLORS = len(FLASH_BINS)
else:
    logging.critical('Error with sizes of FLASH_BINS and GROUP_BINS')
    sys.exit(-1)

# Domain options: extent, map & data res, flags
DOMAIN_OPTS = {
    'fulldisk': {
        'LLlat':    -55.0,	# Lower-Left latitude
        'URlat':     55.0,	# Upper-Right Latitude
        'LLlon':   -130.2,	# Lower-Left Longitude
        'URlon':    -20.2,	# Upper-Right Longitude
        'mapres':  '110m',	# cartopy map resolution
        'resKm':     16,	# resolution of output PNG and GeoTIFF
        'sendPngTo':  1,	# send PNG to PNG_DIR
        'sendGtifTo': 1		# send GeoTIFF to GTIF_DIR
        },
    'conus': {
        'LLlat':     20.0,	# Lower-Left latitude
        'URlat':     55.0,	# Upper-Right Latitude
        'LLlon':   -130.0,	# Lower-Left Longitude
        'URlon':    -60.0,	# Upper-Right Longitude
        'mapres':   '50m',	# cartopy map resolution
        'resKm':      8,	# resolution of output PNG and GeoTIFF
        'sendPngTo':  1,	# send PNG to PNG_DIR
        'sendGtifTo': 1		# send GeoTIFF to GTIF_DIR
        },
    'southamer': {
        'LLlat':    -55.0,	# Lower-Left latitude
        'URlat':     20.0,	# Upper-Right Latitude
        'LLlon':    -90.0,	# Lower-Left Longitude
        'URlon':    -30.0,	# Upper-Right Longitude
        'mapres':  '110m',	# cartopy map resolution
        'resKm':     16,	# resolution of output PNG and GeoTIFF
        'sendPngTo':  0,	# send PNG to PNG_DIR
        'sendGtifTo': 0		# send GeoTIFF to GTIF_DIR
        },
    }

COLOR_DICT = {
    'red': [(0.0, 0.0, 0.0),      #1
            (0.083, 0.0, 1.0),    #2
            (0.167, 1.0, 1.0),    #3
            (0.25, 1.0, 1.0),     #4
            (0.33, 1.0, 1.0),     #5
            (0.42, 1.0, 1.0),     #6
            (0.5, 1.0, 1.0),      #7
            (0.58, 1.0, 0.734),   #8
            (0.667, 0.734, 0.797),#9
            (0.75, 0.797, 1.0),   #10
            (0.833, 1.0, 0.5),    #11
            (0.917, 0.5, 0.25),   #12
            (1.0, 0.25, 0.0)      #13
           ],
    'green': [(0.0, 0.0, 0.0),     #1
              (0.083, 0.0, 1.0),   #2
              (0.167, 1.0, 0.836), #3
              (0.25, 0.836, 0.641),#4
              (0.33, 0.641, 0.469),#5
              (0.42, 0.469, 0.27), #6
              (0.5, 0.27, 0.0),    #7
              (0.58, 0.0, 0.0),    #8
              (0.667, 0.0, 0.0),   #9
              (0.75, 0.0, 0.0),    #10
              (0.833, 0.0, 0.0),   #11
              (0.917, 0.0, 0.0),   #12
              (1.0, 0.0, 0.0),     #13
             ],
    'blue': [(0.0, 0.0, 0.0),      #1
             (0.083, 0.0, 0.0),    #2
             (0.167, 0.0, 0.0),    #3
             (0.25, 0.0, 0.0),     #4
             (0.33, 0.0, 0.0),     #5
             (0.42, 0.0, 0.004),   #6
             (0.5, 0.004, 0.0),    #7
             (0.58, 0.0, 0.145),   #8
             (0.667, 0.145, 0.399),#9
             (0.75, 0.399, 1.0),   #10
             (0.833, 1.0, 0.5),    #11
             (0.917, 0.5, 0.75),   #12
             (1.0, 0.75, 0.0)      #13
            ]
    }
# Generate colormap
GLM_COLORS = LinearSegmentedColormap('colors', COLOR_DICT)
