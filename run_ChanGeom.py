#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 21:29:42 2017

@author: Bodo Bookhagen, August 2017
"""

import pandas as pd
import numpy as np
import os, time
import configparser
import ChanGeom

### Load configuration from ChanGeom.ini
config = configparser.ConfigParser()
ChanGeom_ini_fn = "ChanGeom.ini"
#default init filename - create separate filenames for different regions or parameters
if os.path.exists(ChanGeom_ini_fn) == True:
    config.read(ChanGeom_ini_fn)
else:
    print('Config File: %s does not exist (wrong path?)' %(ChanGeom_ini_fn))

#loading parameters and setting parameters for ChanGeom run
#adding IN_PATH to filenames to make sure we have the absolute paths
IN_KML_PATH                 = os.path.join(config['PATH']['IN_PATH'], config['PATH']['IN_KML_PATH'])
OUT_SHAPEFILE_PATH          = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OUT_SHAPEFILE_PATH'])
OUT_TIF_PATH                = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OUT_TIF_PATH'])
OUT_CENTERLINE_TIF_PATH     = os.path.join(config['PATH']['IN_PATH'],config['PATH']['OUT_CENTERLINE_TIF_PATH'])
OUT_CENTERLINE_SHAPE_PATH   = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OUT_CENTERLINE_SHAPE_PATH'])
CTRLINE_SHP_FN              = os.path.join(config['PATH']['IN_PATH'], config['PATH']['REGION_LABEL']+config['PATH']['CTRLINE_SHP_FN'])
SRTM1_SHP_FN                = os.path.join(config['PATH']['IN_PATH'], config['PATH']['REGION_LABEL']+config['PATH']['SRTM1_SHP_FN'])
SHAPEFILE_MERGED_OUT        = SRTM1_SHP_FN
HDF_FNAME                   = os.path.join(config['PATH']['IN_PATH'], config['PATH']['REGION_LABEL']+config['PATH']['HDF_FNAME'])
ASCII_FNAME                 = os.path.join(config['PATH']['IN_PATH'], config['PATH']['REGION_LABEL']+config['PATH']['ASCII_FNAME'])
DATA_PATH2SAVE              = os.path.join(config['PATH']['IN_PATH'], config['PATH']['DATA_PATH2SAVE'])
PLOT_PATH2SAVE              = os.path.join(config['PATH']['IN_PATH'], config['PATH']['PLOT_PATH2SAVE'])
OVERVIEW_FIG_FN             = os.path.join(config['PATH']['IN_PATH'], config['PATH']['OVERVIEW_FIG_FN'])

### Main Run starts here
t_start = time.time()
ChanGeom.convert_KML(OUT_SHAPEFILE_PATH, OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, 
                OUT_CENTERLINE_SHAPE_PATH, IN_KML_PATH, config['PARAMETERS']['PROJ4_STRING'], config['COMMAND_FNAMES']['ogr2ogr_command'])

outSpatialRef = ChanGeom.polygon_to_raster(OUT_SHAPEFILE_PATH, OUT_TIF_PATH, 
                                         config['COMMAND_FNAMES']['gdal_rasterize_command'],
                                         float(config['PARAMETERS']['PIXEL_SIZE']))
ChanGeom.find_centerline(OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, 
                         config['COMMAND_FNAMES']['matlab_command'], OUT_CENTERLINE_SHAPE_PATH,
                         config['COMMAND_FNAMES']['gdal_translate_command'], outSpatialRef)
ChanGeom.merge_shapefiles(SHAPEFILE_MERGED_OUT, config['COMMAND_FNAMES']['shapemerger_command'], 
                          config['PATH']['IN_PATH'])

#STILL need to add:
# Call matlab <label>__DEM_to_FAC_STR.m
#ADD projection information to output from Matlab
#ogr2ogr -s_srs epsg:32644 -t_srs epsg:32644 WHimalaya_srtm1_1km2_UTM44N_WGS84.shp WHimalaya_srtm1_1km2.shp

ChanGeom.postprocess_shapefiles(config['PATH']['IN_PATH'], SRTM1_SHP_FN, CTRLINE_SHP_FN, 
                                float(config['PARAMETERS']['DISTANCE_THRESHOLD']), 
                                float(config['PARAMETERS']['SEARCH_RADIUS']), DATA_PATH2SAVE, 
                                PLOT_PATH2SAVE, OVERVIEW_FIG_FN, HDF_FNAME, ASCII_FNAME)

print("Processing time %s s" % str(ChanGeom.prettyfloat(time.time() - t_start)) )
print("Processing time %s m" % str(ChanGeom.prettyfloat((time.time() - t_start)/60.)) )
print("Processing time %s h" % str(ChanGeom.prettyfloat((time.time() - t_start)/(60.*60.))) )
