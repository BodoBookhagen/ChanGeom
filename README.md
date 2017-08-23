# ChanGeom
ChanGeom - Channel Geometry, River Width, and along-stream distance extraction from KML files. The KML files can be created in Google Earth (manual clicking along river paths) or by automatic extraction from airphotos or satellite imagery. This is a modified and enhanced version originally described in:

Fisher, G.B., et al., Channel planform geometry and slopes from freely available high-spatial resolution imagery and DEM fusion: Implications for channel width scalings, erosion proxies, and fluvial signatures in tectonically active landscapes, Geomorphology (2013), http://dx.doi.org/10.1016/j.geomorph.2013.04.011

Please cite that manuscript when using these scripts.

Contact Bodo Bookhagen (bodo.bookhagen@uni-potsdam.de) for questions and comments pertinent to the code.

The following major changes have been done to the original code:
- Full Python 3.6 implementation
- gdal and ogr workflows to read and write large number of shapefiles and tif files
- automatic (re-) projection of input files and fully automatized workflows
- Scipy-based channel extraction (no Matlab required) *_WILL BE UPDATED to allow CUDA and multi-core processing_*
- merging all channel width shapefiles into one coherent shapefile for each region for easier analysis
- linking channel width file with topographic information from DEM file (drainage area, flow distance, slope) through Kd trees (no ESRI arcpy or ArcMap necessary). This is a significant speed improvement over arcpy.
- generating standard output plots through matplotlib to visualize and analyze results
- Generating a Python-Pandas DataFrame with all data for efficient and fast analysis
- Saving all data to ASCII and HDF files for further analysis in Python-Pandas or Matlab


## Installation
The code will run on any operating system. You will need to install a number of packages for python, but these should be standard on any Python 3.X installation (e.g., Anaconda, WinPython): gdal, ogr, numpy, scipy, matplotlib, pandas

## Running the Code
You will only need to change/set parameters in the _ChanGeom.ini_ file. These will need to be modifed for each region. See the [ChanGeom.ini](ChanGeom.ini) file for more information. You should create a separate ChanGeom.ini file for each region.

## Example
Will follow.
