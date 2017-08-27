#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 21:29:42 2017

@author: Bodo Bookhagen, August 2017
This code requires a series of modules. If you have conda/miniconda installed, use:
conda create -n py35 python=3.5 scipy pandas numpy matplotlib gdal scikit-image gdal ipython python=3.5.0 spyder scikit-fmm

#run with:
source activate py35
python3 /home/bodo/Dropbox/soft/ChanGeom/run_ChanGeom.py --i ChanGeom_Nahan_salient.ini
"""

import pandas as pd
import numpy as np
import os, time
import configparser
import argparse
import ChanGeom

parser = argparse.ArgumentParser(description='Convert KML/KMZ to channel width and store the channel width and distance along channel. Requires a driver file with parameters - if not given, assumes ini file is called "ChanGeom.ini". See https://github.com/UP-RS-ESP/ChanGeom for more information (B. Bookhagen, bodo.bookhagen@uni-potsdam.de, v0.1, 25-Aug-2017).')
parser.add_argument('--i', dest="ChanGeom_ini_fn", default = "ChanGeom.ini", help='Definition/Driver/Initialization file to be processed (e.g., /raid/bodo/data/Himalaya/CentralNepal/ChanGeom_CentralNepal.ini')

args = parser.parse_args()
ChanGeom_ini_fn = args.ChanGeom_ini_fn
#ChanGeom_ini_fn='/home/bodo/Dropbox/Himalaya/ChanGeom_Himalaya/Nahan_salient/ChanGeom_Nahan_salient.ini'
#
### Load configuration from ChanGeom.ini
config = configparser.ConfigParser()

if os.path.exists(ChanGeom_ini_fn) == True:
    config.read(ChanGeom_ini_fn)
else:
    print('Config File: Can not find %s (wrong filename or path?)' %(ChanGeom_ini_fn))

#loading parameters and setting parameters for ChanGeom run
#adding IN_PATH to filenames to make sure we have the absolute paths
IN_KML_PATH                 = os.path.join(config['PATH']['IN_PATH'], config['PATH']['IN_KML_PATH'])
OUT_SHAPEFILE_PATH          = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OUT_SHAPEFILE_PATH'])
OUT_TIF_PATH                = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OUT_TIF_PATH'])
OUT_CENTERLINE_TIF_PATH     = os.path.join(config['PATH']['IN_PATH'],config['PATH']['OUT_CENTERLINE_TIF_PATH'])
OUT_CENTERLINE_SHAPE_PATH   = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OUT_CENTERLINE_SHAPE_PATH'])
CTRLINE_SHP_FN              = os.path.join(config['PATH']['IN_PATH'], config['PATH']['REGION_LABEL']+config['PATH']['CTRLINE_SHP_FN'])
HDF_FNAME                   = os.path.join(config['PATH']['IN_PATH'], config['PATH']['REGION_LABEL']+config['PATH']['HDF_FNAME'])
ASCII_FNAME                 = os.path.join(config['PATH']['IN_PATH'], config['PATH']['REGION_LABEL']+config['PATH']['ASCII_FNAME'])
DATA_PATH2SAVE              = os.path.join(config['PATH']['IN_PATH'], config['PATH']['DATA_PATH2SAVE'])
PLOT_PATH2SAVE              = os.path.join(config['PATH']['IN_PATH'], config['PATH']['PLOT_PATH2SAVE'])
OVERVIEW_FIG_FN             = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OVERVIEW_FIG_FN'])
shapemerger_command         = os.path.join(config['PATH']['CODE_PATH'], config['COMMAND_FNAMES']['shapemerger_command'])

### Main Run starts here
t_start = time.time()
outSpatialRef = ChanGeom.convert_KML(OUT_SHAPEFILE_PATH, OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, 
                OUT_CENTERLINE_SHAPE_PATH, IN_KML_PATH, 
                config['PARAMETERS']['PROJ4_STRING'], 
                config['COMMAND_FNAMES']['ogr2ogr_command'])

ChanGeom.polygon_to_raster(OUT_SHAPEFILE_PATH, OUT_TIF_PATH, 
                                         config['COMMAND_FNAMES']['gdal_rasterize_command'],
                                         float(config['PARAMETERS']['PIXEL_SIZE']))
                
ChanGeom.find_centerline_PYTHON(OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, 
                         config['PATH']['CODE_PATH'], 
                         , 
                         OUT_CENTERLINE_SHAPE_PATH, 
                         config['COMMAND_FNAMES']['gdal_translate_command'], 
                         outSpatialRef)

ChanGeom.find_centerline_MATLAB(OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, 
                         config['PATH']['CODE_PATH'], 
                         config['COMMAND_FNAMES']['matlab_command'], 
                         OUT_CENTERLINE_SHAPE_PATH, 
                         config['COMMAND_FNAMES']['gdal_translate_command'], 
                         outSpatialRef)

ChanGeom.prepare_DEM(config['PATH']['CODE_PATH'], 
                     config['COMMAND_FNAMES']['matlab_command'], 
                     config['COMMAND_FNAMES']['ogr2ogr_command'], 
                     outSpatialRef,
                     config['PATH']['DEM_FNAME'], 
                     config['PATH']['DEM_SHAPE_OUT_FNAME'],
                     config['PATH']['TOPOTOOLBOX_PATH'],
                     config['PATH']['DEM_AREA_THRESHOLD'])

ChanGeom.merge_CTRL_shapefiles(CTRLINE_SHP_FN, shapemerger_command, 
                          config['PATH']['IN_PATH'])


ChanGeom.postprocess_shapefiles(config['PATH']['IN_PATH'], 
                                config['PATH']['DEM_SHAPE_OUT_FNAME'], 
                                CTRLINE_SHP_FN, 
                                float(config['PARAMETERS']['DISTANCE_THRESHOLD']), 
                                float(config['PARAMETERS']['SEARCH_RADIUS']), 
                                DATA_PATH2SAVE, 
                                PLOT_PATH2SAVE, OVERVIEW_FIG_FN, HDF_FNAME, 
                                ASCII_FNAME)

print("Processing time %0.2f s" % (time.time() - t_start) )
print("Processing time %0.2f m" % ((time.time() - t_start)/60.) )
print("Processing time %0.2f h" % ((time.time() - t_start)/(60.*60.)) )
