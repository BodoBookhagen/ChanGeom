# -*- coding: utf-8 -*-
#Prepare data for ChanGeom

# ChanGeom
# Script prepares KML files and performs the following steps:
# 1. Take KML and re-project (e.g., into UTM or equal area projection)
# 2. store as polygon shapefile 
# 3. Rasterize polygon shapefile
# 4. Run Matlab-based RiverWidth algorithm
# 5. combine outputs into a point Shapefile with two Fields: RWidth_m and CumDist_m
#
# Edit the variables below to match your system and set pixel_size
#
# The script will determine if the output file(s) exists and will not delete these files.
# That is, if a file has been previously processed, you want to delete it if
# you want to re-run the script.
#
# Created May 2014 by Bodo Bookhagen
# Modified Dec 2015 by Bodo Bookhagen
# Modified Aug 2017 by Bodo Bookhagen


# Import gdal and ogr modules
try:
    from osgeo import gdal, ogr, osr
except ImportError:
    import gdal, ogr, osr
import os, subprocess, sys, csv

import glob
import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import pandas as pd
from skimage.morphology import skeletonize, medial_axis
from skimage.segmentation import find_boundaries
from scipy.ndimage.morphology import distance_transform_edt
from scipy.ndimage import convolve, generate_binary_structure, binary_erosion
#import skfmm

def mk_LRP_DA_subplots(SRTM1_merged_data_sorted, RWidth_data_sorted, 
                       Ctrline_fname_list, Shp_fname_txt, path2save, i):
# Plot LRP profile and width
    if os.path.exists(path2save) == False:
        os.mkdir(path2save)
    plt.clf()
    # create all axes we need
    ax1a = plt.subplot(311)
    ax1b = ax1a.twinx()
    ax2a = plt.subplot(312)
    ax2b = ax2a.twinx()
    ax3a = plt.subplot(313)
    ax3b = ax3a.twinx()
#    ax4a = plt.subplot(224)
#    ax4b = ax4a.twinx()
    
    # SRTM1: Distance
    ax1a.get_shared_x_axes().join(ax1a, ax1b)
    ax1a.set_title('SRTM1: ' + Ctrline_fname_list[i], fontsize=16)
    ax1a.plot(SRTM1_merged_data_sorted[i].RWidth_CumDist_m, 
             SRTM1_merged_data_sorted[i].SRTM1_elevation, 'k-',
                     label='LRP')
    ax1a.set_ylabel('Elevation asl [m]', fontsize=12, color='k')
    ax1b.errorbar(SRTM1_merged_data_sorted[i].RWidth_CumDist_m, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean, color='b',
                     yerr=SRTM1_merged_data_sorted[i].RWidth_m_std, 
                     label='RWidth')
    ax1b.set_ylabel('Riverwidth mean and std. dev [m]', fontsize=12, color='b')
    ax1a.set_xlabel('Cumulative distance from channel top [m]', fontsize=12, color='k')
    ax1a.grid()
    lines, labels = ax1a.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    ax1a.legend(lines + lines2, labels + labels2, loc=0)
    
    # SRTM1: DA
    ax2a.get_shared_x_axes().join(ax2a, ax2b)
    ax2a.set_title('SRTM1: ' + Ctrline_fname_list[i], fontsize=16)
    #ax2a.set_xscale('log')
    #ax2a.set_yscale('log')
    ax2a.plot(SRTM1_merged_data_sorted[i].SRTM1_DA_km2, 
             SRTM1_merged_data_sorted[i].SRTM1_gradient, 'k-', marker='+',
                     label='DA vs slope')
    ax2a.set_xlabel('Drainage Area [km^2]', fontsize=12, color='k')
    ax2a.set_ylabel('Gradient [m/m]', fontsize=12, color='k')
    #ax2b.set_yscale('log')
    ax2b.errorbar(SRTM1_merged_data_sorted[i].SRTM1_DA_km2, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean, color='b',
                     marker='+', yerr=SRTM1_merged_data_sorted[i].RWidth_m_std, 
                     label='RWidth')
    ax2b.set_ylabel('Riverwidth mean and std. dev [m]', fontsize=12, color='b')
    ax2a.grid()
    lines, labels = ax2a.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2a.legend(lines + lines2, labels + labels2, loc=0)

    # RivWidth
    #ax3a.get_shared_x_axes().join(ax3a, ax3b)
    ax3a.set_title('RivWidth: ' + Shp_fname_txt[i], fontsize=16)
    ax3a.plot(np.array(RWidth_data_sorted[i].CumDist_m), 
             np.array(RWidth_data_sorted[i].RWidth_m), '-', color='0.75', 
                     label='RWidth')
    ax3a.set_ylabel('Riverwidth [m]', fontsize=12, color='k')
    ax3a.set_xlabel('Cumulative distance from channel top [m]', fontsize=12, color='k')
    ax3a.grid()
    ax3a.legend()

    f = plt.gcf()  # f = figure(n) if you know the figure number
    #f.set_size_inches(11.69,8.27) #A4
    f.set_size_inches(16.53,11.69) #A3
    fname2save=os.path.join(path2save, Ctrline_fname_list[i] + '_plt.png')
    plt.savefig(fname2save, papertype='a3', orientation='landscape')

def mk_single_plots(SRTM1_merged_data_sorted, SRTM1_merged_data_sorted_df, 
                    Ctrline_fname_list, Ctrline_fname_list_nr ):
    #Plot single plot to see labels
    plt.clf()
    for i in range(Ctrline_fname_list_nr):
        plt.errorbar(SRTM1_merged_data_sorted[i].SRTM1_distance_from_outlet, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean,
                     yerr=SRTM1_merged_data_sorted[i].RWidth_m_std,
                     label=Ctrline_fname_list[i])
    plt.xlabel('Distance from outlet (km)', fontsize=12)
    plt.ylabel('Riverwidth mean and std. dev (m)',fontsize=12)
    plt.yscale('log')
    plt.grid()
    plt.title('Distance from outlet vs. Riverwidth')
    plt.legend()
    f = plt.gcf()  # f = figure(n) if you know the figure number
    f.set_size_inches(16.53,11.69) #A3
    plt.savefig('Centerline_SRTM1_DA_labels.png', papertype='a3', orientation='landscape')
    
    plt.clf()
    for i in range(Ctrline_fname_list_nr):
        plt.errorbar(SRTM1_merged_data_sorted[i].SRTM1_elevation, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean,
                     yerr=SRTM1_merged_data_sorted[i].RWidth_m_std,
                     label=Ctrline_fname_list[i])
    plt.xlabel('Elevation asl (m)', fontsize=12)
    plt.ylabel('Riverwidth mean and std. dev (m)',fontsize=12)
    plt.yscale('log')
    plt.grid()
    plt.title('Distance from outlet vs. Riverwidth')
    plt.legend()
    f = plt.gcf()  # f = figure(n) if you know the figure number
    f.set_size_inches(16.53,11.69) #A3
    plt.savefig('Centerline_SRTM1_elevation_river_width_with_labels.png', papertype='a3', orientation='landscape')
    
    plt.clf()
    for i in range(Ctrline_fname_list_nr):
        plt.errorbar(SRTM1_merged_data_sorted[i].RWidth_CumDist_m, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean,
                     yerr=SRTM1_merged_data_sorted[i].RWidth_m_std,
                     label=Ctrline_fname_list[i])
    plt.xlabel('Cumulative downstream distance (mm)', fontsize=12)
    plt.ylabel('Riverwidth mean and std. dev (m)',fontsize=12)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.title('Cumulative downstream distance vs. Riverwidth')
    plt.legend()
    f = plt.gcf()  # f = figure(n) if you know the figure number
    f.set_size_inches(16.53,11.69) #A3
    plt.savefig('Centerline_SRTM1_CumDistance_downstream_with_labels.png', papertype='a3', orientation='landscape')
    
    plt.clf()
    for i in range(Ctrline_fname_list_nr):
        plt.errorbar(SRTM1_merged_data_sorted[i].SRTM1_maxdistance_from_channelhead, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean,
                     yerr=SRTM1_merged_data_sorted[i].RWidth_m_std,
                     label=Ctrline_fname_list[i])
    plt.xlabel('Max. Distance from channel head (km)', fontsize=12)
    plt.ylabel('Riverwidth mean and std. dev (m)',fontsize=12)
    plt.yscale('log')
    plt.xscale('log')
    plt.grid()
    plt.title('Max. Distance from channel head vs. Riverwidth')
    plt.legend()
    f = plt.gcf()  # f = figure(n) if you know the figure number
    f.set_size_inches(16.53,11.69) #A3
    plt.savefig('Centerline_SRTM1_Distance_max_from_channelhead_with_labels.png', papertype='a3', orientation='landscape')

    #2x2 subplot
    plt.clf()
    plt.subplot(2,2,1)
    plt.scatter(SRTM1_merged_data_sorted_df.SRTM1_UTM_x, SRTM1_merged_data_sorted_df.SRTM1_UTM_y,
                s=SRTM1_merged_data_sorted_df.SRTM1_DA_km2, 
                c=SRTM1_merged_data_sorted_df.SRTM1_DA_km2, label='merged data (DA in km2)')
    plt.grid()
    plt.colorbar()
    plt.title('Drainage Area Map', fontsize=16)
    
    plt.subplot(2,2,2)
    plt.scatter(SRTM1_merged_data_sorted_df.SRTM1_UTM_x, SRTM1_merged_data_sorted_df.SRTM1_UTM_y,
                s=SRTM1_merged_data_sorted_df.RWidth_m_mean, 
                c=SRTM1_merged_data_sorted_df.RWidth_m_mean, label='merged data (RWdith mean)')
    plt.colorbar()
    plt.grid()
    plt.title('Riverwidth Map (mean)', fontsize=16)
    
    plt.subplot(2,2,3)
    for i in range(Ctrline_fname_list_nr):
        plt.errorbar(SRTM1_merged_data_sorted[i].SRTM1_distance_from_outlet, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean, 
                    yerr=SRTM1_merged_data_sorted[i].RWidth_m_std, label=Ctrline_fname_list[i])
    plt.xlabel('Distance from outlet (km)', fontsize=12)
    plt.ylabel('Riverwidth mean and std. dev (m)',fontsize=12)
    plt.yscale('log')
    plt.grid()
    plt.title('Distance from outlet vs. Riverwidth')
    
    plt.subplot(2,2,4)
    for i in range(Ctrline_fname_list_nr):
        plt.errorbar(SRTM1_merged_data_sorted[i].SRTM1_DA_km2, 
                     SRTM1_merged_data_sorted[i].RWidth_m_mean, 
                    yerr=SRTM1_merged_data_sorted[i].RWidth_m_std,
                    label=Ctrline_fname_list[i])
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.xlabel('DA (km2)', fontsize=12)
    plt.ylabel('Riverwidth mean and std. dev (m)',fontsize=12)
    plt.grid()
    plt.title('DA vs. Riverwidth')
    f = plt.gcf()  # f = figure(n) if you know the figure number
    #f.set_size_inches(11.69,8.27) #A4
    f.set_size_inches(16.53,11.69) #A3
    plt.savefig('Centerline_SRTM1_data_plots.png', papertype='a3', orientation='landscape')

def load_SRTM_shp_data(SRTM1_shp):
    dataSource = ogr.Open(SRTM1_shp)
    daLayer = dataSource.GetLayer(0)
    
    SRTM1_DA_km2=[]
    SRTM1_distance=[]
    SRTM1_d_maxdd = []
    SRTM1_gradient=[]
    SRTM1_ksn=[]
    SRTM1_UTM_x=[]
    SRTM1_UTM_y=[]
    SRTM1_elevation=[]
    for i in range(daLayer.GetFeatureCount()):
        feature = daLayer.GetFeature(i)
        #list feature
        #feature.items()   
        SRTM1_elevation.append(feature.GetField("elev"))
        SRTM1_distance.append(feature.GetField("d_out"))
        SRTM1_d_maxdd.append(feature.GetField("d_maxdd"))
        SRTM1_DA_km2.append(feature.GetField("A_km2"))
        SRTM1_gradient.append(feature.GetField("gradient"))
        SRTM1_ksn.append(feature.GetField("ksn"))
        geom = feature.GetGeometryRef()
        SRTM1_UTM_x.append(geom.GetX())
        SRTM1_UTM_y.append(geom.GetY())
    return SRTM1_UTM_x, SRTM1_UTM_y, SRTM1_elevation, SRTM1_d_maxdd, SRTM1_DA_km2, SRTM1_gradient, SRTM1_ksn, SRTM1_distance

def load_Ctrl_shp_data(CtrLine_shp):
    #Get Centerline data
    dataSource = ogr.Open(CtrLine_shp)
    daLayer = dataSource.GetLayer(0)
    
    Shp_fname=[]
    Shp_src_ID=[]
    RWidth_m=[]
    CumDist_m=[]
    CtrLine_UTM_x=[]
    CtrLine_UTM_y=[]
    for i in range(daLayer.GetFeatureCount()):
        feature = daLayer.GetFeature(i)
        #list feature
        #feature.items()   
        Shp_fname.append(feature.GetField("SOURCESHP"))
        Shp_src_ID.append(feature.GetField("SOURCEFID"))
        RWidth_m.append(feature.GetField("RWidth_m"))
        CumDist_m.append(feature.GetField("CumDist_m"))
        geom = feature.GetGeometryRef()
        CtrLine_UTM_x.append(geom.GetX())
        CtrLine_UTM_y.append(geom.GetY())
    return CtrLine_UTM_x, CtrLine_UTM_y, RWidth_m, CumDist_m, Shp_fname, Shp_src_ID
###
###
# Code starts here - no editing below this line

def raster2point(srcfile1, srcfile2, dstfile): #Write tif to CSV, ignoring 0s
    # based on gdal2xyz.py
    # reads two tif files and collects data from non-0 values 
    srcwin = None
    skip = 1    
    delim = ','

    # Open source file. 
    srcds1 = gdal.Open( srcfile1 )
    if srcds1 is None:
        print('Could not open %s.' % srcfile1)
        sys.exit( 1 )

    band = srcds1.GetRasterBand(1)
    if band is None:
        print('Could not get band %d' % band_num)
        sys.exit( 1 )
    bands = []
    bands.append(band)
    gt = srcds1.GetGeoTransform()

    srcds2 = gdal.Open( srcfile2 )
    if srcds2 is None:
        print('Could not open %s.' % srcfile2)
        sys.exit( 1 )

    band2 = srcds2.GetRasterBand(1)
    if band2 is None:
        print('Could not get band %d' % band_num)
        sys.exit( 1 )
    bands2 = []
    bands2.append(band)
    gt = srcds2.GetGeoTransform()
  
    # Collect information on all the source files.
    srcwin1 = (0,0,srcds1.RasterXSize,srcds1.RasterYSize)

    # Open the output file.
    if dstfile is not None:
        dst_fh = open(dstfile,'wt')
    else:
        print("Can not open destination file: %s" % dstfile)
        sys.exit( 1 )

    band_format = (("%g" + delim + "%g" + delim) * len(bands)).rstrip(delim) + '\n'

    # Setup an appropriate print format.
    if abs(gt[0]) < 180 and abs(gt[3]) < 180 \
       and abs(srcds1.RasterXSize * gt[1]) < 180 \
       and abs(srcds1.RasterYSize * gt[5]) < 180:
        format = '%.10g' + delim + '%.10g' + delim + '%s' + delim + '%s'
    else:
        format = '%.3f' + delim + '%.3f' + delim + '%s' + delim + '%s'
    dst_fh.write( 'X,Y,RWidth_m,CumDist_m\n' )

    # Loop emitting data.
    for y in range(srcwin1[1],srcwin1[1]+srcwin1[3],skip):
        data1 = []
        data2 = []
        for band in bands:
            band_data1 = band.ReadAsArray( srcwin1[0], y, srcwin1[2], 1 )    
            band_data1 = np.reshape( band_data1, (srcwin1[2],) )
            band_data2 = band2.ReadAsArray( srcwin1[0], y, srcwin1[2], 1 )    
            band_data2 = np.reshape( band_data2, (srcwin1[2],) )
            
            data1.append(band_data1)
            data2.append(band_data2)

        for x_i in range(0,srcwin1[2],skip):
            x = x_i + srcwin1[0]
            if data1[0][x_i] == 0:
                continue
            geo_x = gt[0] + (x+0.5) * gt[1] + (y+0.5) * gt[2]
            geo_y = gt[3] + (x+0.5) * gt[4] + (y+0.5) * gt[5]

            x_i_data1 = []
            x_i_data2 = []
            for i in range(len(bands)):
                x_i_data1.append(data1[i][x_i])
                x_i_data2.append(data2[i][x_i])
            band_str1 = "%g" % tuple(x_i_data1)
            band_str2 = "%g\n" % tuple(x_i_data2)
            line = format % (float(geo_x),float(geo_y), band_str1, band_str2)
            dst_fh.write( line )

def csv2shapefile( csvfile, shapefile, spatialreference ):
    # use a dictionary reader so we can access by field name
    reader = csv.DictReader(open(csvfile,"rt"),
        delimiter=',',
        quoting=csv.QUOTE_NONE)
    
    # set up the shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(shapefile)
    # create the layer
    layer = data_source.CreateLayer("RWidth", spatialreference, ogr.wkbPoint)
    
    # Add the fields we're interested in
    layer.CreateField(ogr.FieldDefn("Y", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("X", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("RWidth_m", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("CumDist_m", ogr.OFTReal))
    
    # Process the text file and add the attributes and features to the shapefile
    for row in reader:
      feature = ogr.Feature(layer.GetLayerDefn())
      # Set the attributes using the values from the delimited text file
      feature.SetField("X", row['X'])
      feature.SetField("Y", row['Y'])
      feature.SetField("RWidth_m", row['RWidth_m'])
      feature.SetField("CumDist_m", row['CumDist_m'])
    
      # create the WKT for the feature using Python string formatting
      wkt = "POINT(%f %f)" %  (float(row['X']) , float(row['Y']))
    
      # Create the point from the Well Known Txt
      point = ogr.CreateGeometryFromWkt(wkt)
    
      # Set the feature geometry using the point
      feature.SetGeometry(point)
      # Create the feature in the layer (shapefile)
      layer.CreateFeature(feature)
      # Destroy the feature to free resources
      feature.Destroy()
    
    # Destroy the data source to free resources
    data_source.Destroy()
        
def convert_KML(OUT_SHAPEFILE_PATH, OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, 
                OUT_CENTERLINE_SHAPE_PATH, IN_KML_PATH, PROJ4_STRING, ogr2ogr_command):
    if not os.path.exists(OUT_SHAPEFILE_PATH):
        os.makedirs(OUT_SHAPEFILE_PATH)
    if not os.path.exists(OUT_TIF_PATH):
        os.makedirs(OUT_TIF_PATH)
    if not os.path.exists(OUT_CENTERLINE_TIF_PATH):
        os.makedirs(OUT_CENTERLINE_TIF_PATH)
    if not os.path.exists(OUT_CENTERLINE_SHAPE_PATH):
        os.makedirs(OUT_CENTERLINE_SHAPE_PATH)
    
    for file in os.listdir(IN_KML_PATH):
         if file.endswith('.kml'):
              vector_file = os.path.join(IN_KML_PATH, file)
              print("Converting and Reprojecting %s" %(vector_file))
              
              #project KML into equal-area projection (UTM, LETIBET, LESAMRR)
              #define output projection
              outSpatialRef = osr.SpatialReference()
              #if there is a *.prj file in the PATH directory, take the projection information from that file
              #Can read ESRI and PROJ4 projection information
              #A .prj file can be created with: gdalsrsinfo -o proj4 landsat.tif >UTM19S_WGS84.prj
              #gdalsrsinfo -o proj4 "EPSG:32719" >UTM19S_WGS84.prj
              #this will need to be stored in the IN_KML_PATH file
              if len(glob.glob(os.path.join(IN_KML_PATH, '*.prj'))) > 0:
                  #there is a .prj file
                  for projfile in os.listdir(IN_KML_PATH):
                      if projfile.endswith('.prj'):
                          prj_file = os.path.join(IN_KML_PATH, projfile)
                          #import projection file. FIRST try if this is an ESRI prj file
                          prj_text = open(prj_file, 'r').read()
                          try:
                              err = outSpatialRef.ImportFromESRI([prj_text])
                          except err != 0:
                              err = outSpatialRef.ImportFromProj4(prj_text)
                              if err != 0:
                                  raise ValueError("Error importing ESRI or PROJ4 projection information from: %s" % prj_file)
                                  break
                          if err == 0:
                              #print outSpatialRef.ExportToProj4()
                              print
                  # Set Projection name to name of prj file
                  outSpatialRef.SetProjCS(prj_file.split('/')[-1][0:-4]);
              else:
                  err = outSpatialRef.ImportFromProj4(PROJ4_STRING)
                  if err != 0:
                      raise ValueError("Error importing ESRI or PROJ4 projection information from: %s" % prj_file)
                      break
                  
              #verify if OUT_SHAPEFILE_PATH exists
              if not os.path.exists(OUT_SHAPEFILE_PATH):
                  os.makedirs(OUT_SHAPEFILE_PATH)
              if os.path.exists(OUT_SHAPEFILE_PATH) == False:
                  print("Can not create directory: %s" %OUT_SHAPEFILE_PATH)
    
              out_shapefile = os.path.join(OUT_SHAPEFILE_PATH, file.split('.')[0] + "_projected.shp")
              if os.path.exists(out_shapefile):   
                  print("%s exists, skipping to next file..." %out_shapefile)
              else:    
                  #call ogr2ogr to convert KML
                  out_shapefile_txt = out_shapefile + '.ogr2ogr.out'
                  ogr2ogr_subprocess_command = ogr2ogr_command + ' -a_srs EPSG:4326 -t_srs ' + '"' + outSpatialRef.ExportToProj4() + '"' + ' -f "ESRI Shapefile" ' + out_shapefile + ' ' + vector_file + '>' + out_shapefile_txt
                  os.system(ogr2ogr_subprocess_command)

    return outSpatialRef
    
def polygon_to_raster(OUT_SHAPEFILE_PATH, OUT_TIF_PATH, gdal_rasterize_command, pixel_size):
    ### Process: Polygon to Raster
    # Load projected shapefile and convert to grid
    for file in os.listdir(OUT_SHAPEFILE_PATH):
         if file.endswith('.shp'):
              vector_file = os.path.join(OUT_SHAPEFILE_PATH, file)
              print("Rasterizing %s" % vector_file )
    
              #verify if OUT_TIF_PATH exists
              if not os.path.exists(OUT_TIF_PATH):
                  os.makedirs(OUT_TIF_PATH)
              if os.path.exists(OUT_TIF_PATH) == False:
                  print("Can not create directory: %s" %OUT_TIF_PATH)
    
              out_tiffile = os.path.join(OUT_TIF_PATH, file.split('.')[0] + "_r" + str(pixel_size) + "m.tif")
    
              if os.path.exists(out_tiffile):   
                  print("%s exists, skipping to next file..." %out_tiffile)
              else:    
                  #call gdal_rasterize to rasterize shapefile
                  out_tiffile_txt = out_tiffile.split('.')[0] + '.gdal_rasterize.out'
                  gdal_rasterize_subprocess_command = gdal_rasterize_command + ' -co COMPRESS=DEFLATE -co ZLEVEL=9 -co "NBITS=1" -burn 1 -tap -a_nodata 0 -ot Byte -tr ' + str(pixel_size) + ' ' + str(pixel_size) + ' ' + vector_file + ' ' + out_tiffile + '>'  + out_tiffile_txt
                  os.system(gdal_rasterize_subprocess_command)

def prepare_DEM(CODE_PATH, matlab_command, ogr2ogr_command, 
                outSpatialRef, DEM_FNAME, DEM_SHAPE_OUT_FNAME, 
                TOPOTOOLBOX_PATH, DEM_AREA_THRESHOLD):
    ### start Matlab script to generate shapefile and KML for river channels from DEM
    cwd = os.getcwd()
    os.chdir(CODE_PATH)
    print("Prepare DEM: Derive SHP and KML from DEM: %s" %DEM_FNAME.split('/')[-1])
    if os.path.exists(DEM_SHAPE_OUT_FNAME):   
        print("%s exist, skipping to file conversion..." %(DEM_SHAPE_OUT_FNAME.split('/')[-1]) )
    else:    
        #call Matlab to create centerline and cumulative distance TIFs
        #function prepare_DEM(DEM_fname, SHAPE_OUT_FNAME, TOPOTOOLBOX_PATH, Area_threshold)
        subprocess.call([matlab_command+" -nosplash -nodisplay -r \"prepare_DEM(\'%s\',\'%s\',\'%s\',%s)\"" % (DEM_FNAME, DEM_SHAPE_OUT_FNAME, TOPOTOOLBOX_PATH, DEM_AREA_THRESHOLD)], shell=True);
    
    print("Converting Shapefile to KML: %s" %DEM_SHAPE_OUT_FNAME.split('/')[-1])
    #adding .prj file for shapefile 
    if os.path.exists(DEM_SHAPE_OUT_FNAME[:-3] + 'prj') == False:
        DEM_SHAPE_OUT_FNAME_PRJ = DEM_SHAPE_OUT_FNAME[:-3] + 'prj'
        out_fid = open(DEM_SHAPE_OUT_FNAME_PRJ, 'w')
        out_fid.write(outSpatialRef.ExportToWkt())
        out_fid.close()
    DEM_KML_OUT_FNAME = DEM_SHAPE_OUT_FNAME
    DEM_KML_OUT_FNAME = DEM_KML_OUT_FNAME[:-3] + 'kml'
    if os.path.exists(DEM_KML_OUT_FNAME):   
        print("%s exist, skipping to next file..." %(DEM_KML_OUT_FNAME) )
    else:    
        #call ogr2ogr to add spatial information/projection to shapefile
        #Matlab's shapewrite doesn't add spatial information
        ogr2ogr_subprocess_command = ogr2ogr_command  + ' -f KML ' + DEM_KML_OUT_FNAME + ' ' + DEM_SHAPE_OUT_FNAME
        os.system(ogr2ogr_subprocess_command)

    os.chdir(cwd)

def find_centerline_PYTHON(OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, OUT_CENTERLINE_SHAPE_PATH,
                    outSpatialRef):
    #use SKIMAGE to skeletonize image
    for file in os.listdir(OUT_TIF_PATH):
         if file.endswith('.tif'):
              in_tiffile = os.path.join(OUT_TIF_PATH, file)
              print("Finding Centerline %s" %in_tiffile)
    
              #verify if OUT_CENTERLINE_TIF_PATH exists
              if not os.path.exists(OUT_CENTERLINE_TIF_PATH):
                  os.makedirs(OUT_CENTERLINE_TIF_PATH)
              if os.path.exists(OUT_CENTERLINE_TIF_PATH) == False:
                  print("Can not create directory: %s" % OUT_CENTERLINE_TIF_PATH)
    
              out_tiffile1 = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + "_centerline.tif")
              out_tiffile2 = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + "_centerline_chwidth.tif")
              out_tiffile3 = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + "_centerline_cumdist.tif")              
    
              if os.path.exists(out_tiffile1) and os.path.exists(out_tiffile2) and os.path.exists(out_tiffile3):   
                  print("%s  and %s exist, skipping to next file..." %(out_tiffile1.split('/')[-1], out_tiffile2.split('/')[-1]))
              else:
                  #load image
                  ds = gdal.Open(in_tiffile)
                  in_array = np.array(ds.GetRasterBand(1).ReadAsArray())
                  gt = ds.GetGeoTransform()
                  cs = ds.GetProjection()
                  cs_sr = osr.SpatialReference()
                  cs_sr.ImportFromWkt(cs)
                 
                  #Skeletonize, i.e. thin image
                  skeleton = skeletonize(in_array)
                  #medial axis results appear to be less useful
                  #medial_axis_array, medial_distance = medial_axis(skeleton, mask=None, return_distance=True)
                  #Get boundary or outline (original polygon)
                  boundary = find_boundaries(in_array, connectivity=1, mode='inner', background=0)
                  skeleton_and_boundary = skeleton+boundary
                  #find endpoints on merged skeleton and boundary binary images. 
                  #All 'real' channel endpoints of the stream are connected to the boundary. 
                  #All non-connected endpoints are spurious results and will be pruned.
                  endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy = find_network_endpoints(skeleton_and_boundary)
                  #endpoints not connected to boundary will be removed - these are spurs
                  while len(endpoints_idx) > 0:
                      #print('Numer of endpoints: %d' %len(endpoints_idx))
                      skeleton.ravel()[endpoints_idx] = 0
                      skeleton_and_boundary = skeleton + boundary
                      endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy = find_network_endpoints(skeleton_and_boundary)
                      
                  #remove pixels with more than 2 neighbors - e.g. endpoints of T-like structures
                  skeleton = bw_remove_3nn(skeleton)

                  #Make sure there are no further endpoints
                  endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy = find_network_endpoints(skeleton)
                  while len(endpoints_idx) > 2:
                      #print('Numer of endpoints: %d' %len(endpoints_idx))
                      skeleton.ravel()[endpoints_idx] = 0
                      endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy = find_network_endpoints(skeleton)
                  #Now find endpoints of pruned binary skeleton
                  #endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy = find_network_endpoints(skeleton)
                  
                  #remove mean channel width number of pixels from endpoints
#                  mask = ~in_array.astype(bool)
#                  distance_from_centerline_to_boundary = np.ma.masked_array(distance_transform_edt(in_array), mask)
#                  distance_to_remove = np.int(np.round(np.max(np.ma.masked_array(distance_from_centerline_to_boundary, skeleton != 1))))
#                  for i in range(distance_to_remove):
#                      skeleton.ravel()[endpoints_idx] = 0
#                      endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy = find_network_endpoints(skeleton)
#
#                  endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy = find_network_endpoints(skeleton)
                  print('Pruned binary skeleton and found endpoints')
                  
                  #Alternative approach to filter skeleton: find branch points and go from there. Not very efficient
                  #branch_pts = find_branch_points(skeleton)
                  
                  #Calculate distance transform between centerline and boundaries
                  #mask = ~in_array.astype(bool)
                  #m = np.ones_like(in_array)
                  #m[skeleton==1] = 0
                  #distance_from_boundary_to_centerline = np.ma.masked_array(distance_transform_edt(m), mask)
                  
                  #distance from centerline
                  mask = ~in_array.astype(bool)
                  distance_from_centerline_to_boundary = np.ma.masked_array(distance_transform_edt(in_array), mask)
                  distance_from_centerline_to_boundary = np.ma.masked_array(distance_from_centerline_to_boundary, skeleton != 1)
                  #arbitrarily set start point to first endpoint - will determine direction of flow later
                  #the following definition corresponds to the northern point as start point
                  #note that x and y coordinates are reversed, because of image direction: Y starts counting at 0 at upper left corner, X at lower left corner
                  idx_start_point_x = endpoints_idxx[0]
                  idx_start_point_y = endpoints_idxy[0]
                  idx_end_point_x = endpoints_idxx[1]
                  idx_end_point_y = endpoints_idxy[1]
                  cell_distance_cs_m, cell_downstream_list2D, cell_downstream_list1D, sinuosity, cell_distance_grid = trace_path_distance(skeleton, idx_start_point_x, idx_start_point_y, idx_end_point_x, idx_end_point_y)
                  
                  centerline_chwidth = np.zeros_like(distance_from_centerline_to_boundary)
                  centerline_chwidth[skeleton==1] = distance_from_centerline_to_boundary[skeleton == 1]
                  
                  #write outpus to geotiff grid files
                  driver = gdal.GetDriverByName('GTiff')
                  driver.Register()
                  outRaster = driver.Create(out_tiffile1, skeleton.shape[1], skeleton.shape[0], 1, gdal.GDT_Byte)
                  outRaster.SetGeoTransform(gt)
                  outRaster.SetProjection(cs)
                  outband = outRaster.GetRasterBand(1)
                  outband.SetMetadataItem('Band', 'Centerline')
                  outband.SetNoDataValue(0)
                  outband.WriteArray(skeleton,0,0)
                  outband.FlushCache()
                  del driver, outRaster

                  driver = gdal.GetDriverByName('GTiff')
                  driver.Register()
                  outRaster = driver.Create(out_tiffile2, centerline_chwidth.shape[1], centerline_chwidth.shape[0], 1, gdal.GDT_Float32)
                  outRaster.SetGeoTransform(gt)
                  outRaster.SetProjection(cs)
                  outband = outRaster.GetRasterBand(1)
                  outband.SetMetadataItem('Band', 'Channel_width_m')
                  outband.SetNoDataValue(0)
                  outband.WriteArray(centerline_chwidth,0,0)
                  outband.FlushCache()
                  del driver, outRaster

                  driver = gdal.GetDriverByName('GTiff')
                  driver.Register()
                  outRaster = driver.Create(out_tiffile3, cell_distance_grid.shape[1], distance_from_centerline_to_boundary.shape[0], 1, gdal.GDT_Float32)
                  outRaster.SetGeoTransform(gt)
                  outRaster.SetProjection(cs)
                  outband = outRaster.GetRasterBand(1)
                  outband.SetMetadataItem('Band', 'CumDistance_channel_m')
                  outband.SetNoDataValue(0)
                  outband.WriteArray(cell_distance_grid,0,0)
                  outband.FlushCache()
                  del driver, outRaster

                  fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16.53,11.69),
                                           sharex=True, sharey=True,
                                           subplot_kw={'adjustable': 'box-forced'})
                  ax = axes.ravel()
                  colormap2use = plt.cm.cool
                  distance_map1 = ax[0].imshow(cell_distance_grid, cmap=colormap2use)
                  ax[0].axis('on')
                  ax[0].set_title('cell_distance_grid', fontsize=20)
                  cbar1 = fig.colorbar(distance_map1, orientation='horizontal', ax=ax[0])
                  cbar1.set_label('cumulative distance in # of cells')
                
                  ax[1].imshow(skeleton+boundary, cmap=plt.cm.gray)
                  ax[1].axis('on')
                  ax[1].set_title('skeleton + boundary', fontsize=20)
                
                  ax[2].plot(cell_distance_cs_m, channel_width_m, '+')
                  ax[2].axis('on')
                  ax[2].set_xlabel('Cumulative distance along channel (m)')
                  ax[2].set_ylabel('Channel width (m)')
                  ax[2].set_title('Cell distance vs. channel width', fontsize=20)
                
                  distance_map2 = ax[3].imshow(distance_from_centerline_to_boundary, cmap=colormap2use)
                  ax[3].axis('on')
                  ax[3].set_title('distance_from_centerline_to_boundary', fontsize=20)
                  cbar2 = fig.colorbar(distance_map2, orientation='horizontal', ax=ax[3])
                  cbar2.set_label('distance in # of cells')
                
                  fig.tight_layout()
                  f = plt.gcf()  # f = figure(n) if you know the figure number
                  #f.set_size_inches(11.69,8.27) #A4
                  f.set_size_inches(16.53,11.69) #A3
                  fname2save = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + "_centerline_plt.png")
                  plt.savefig(fname2save, papertype='a3', orientation='landscape')
                  print('Distance analysis completed, saved to file and overview figure has been created.')

                  #create shapefile from coordinates
                  UTM_X = np.arange(gt[0]+gt[1]/2,gt[0]+(skeleton.shape[1]*gt[1]),gt[1], dtype=np.float)
                  UTM_Y = np.arange(gt[3]+gt[5]/2,gt[3]+(skeleton.shape[0]*gt[5]),gt[5], dtype=np.float)
                  channel_width_m = np.zeros_like(cell_distance_cs_m, dtype=np.float)
                  UTM_X_coord = np.zeros_like(cell_distance_cs_m, dtype=np.float)
                  UTM_Y_coord = np.zeros_like(cell_distance_cs_m, dtype=np.float)
                  for i in range(len(cell_downstream_list2D)):
                      channel_width_m[i] = distance_from_centerline_to_boundary.data[cell_downstream_list2D[i,0], cell_downstream_list2D[i,1]]
                      UTM_X_coord[i] = UTM_X[cell_downstream_list2D[i,1]]
                      UTM_Y_coord[i] = UTM_Y[cell_downstream_list2D[i,0]]
                      
                  array2write = np.vstack((UTM_X_coord, UTM_Y_coord, cell_distance_cs_m, channel_width_m, sinuosity)).transpose()
                  #csv file contains UTM_X, UTM_Y, CumulativeDistance_m, ChannelWidth_m, Sinuosity
                  out_csvfile = os.path.join(OUT_CENTERLINE_SHAPE_PATH, '.'.join(file.split('.')[0:-1]) + "_centerline.csv")
                  np.savetxt(out_csvfile, array2write, fmt='%.3f', delimiter=",", header="UTM_X, UTM_Y, CumulativeDistance_m, ChannelWidth_m, Sinuosity")
                  out_shapefile = '.'.join(file.split('.')[0:-1]) + "_centerline.shp"
                  write_to_shapefile(OUT_CENTERLINE_SHAPE_PATH, out_shapefile, array2write, outSpatialRef, shp_layer='Centerline')

                  #smooth:
                  #APPROACH 1:
                  #Smooth over x-m distance and write at center of segment
                  smoothing_interval = 10
                  interp_sample_points = np.arange(0, np.max(cell_distance_cs_m), smoothing_interval)
                  cell_distance_cs_m_interp = np.interp(interp_sample_points, cell_distance_cs_m, cell_distance_cs_m)
                  sinuosity_interp = np.interp(interp_sample_points, cell_distance_cs_m, sinuosity)
                  channel_width_m_interp = np.interp(interp_sample_points, cell_distance_cs_m, channel_width_m)
                  UTM_X_coord_interp = np.interp(interp_sample_points, cell_distance_cs_m, UTM_X_coord)
                  UTM_Y_coord_interp = np.interp(interp_sample_points, cell_distance_cs_m, UTM_Y_coord)
                  array2write = np.vstack((UTM_X_coord_interp, UTM_Y_coord_interp, cell_distance_cs_m_interp, channel_width_m_interp, sinuosity_interp)).transpose()
                  out_csvfile = os.path.join(OUT_CENTERLINE_SHAPE_PATH, '.'.join(file.split('.')[0:-1]) + "_centerline_smooth_%d.csv" %(smoothing_interval))
                  np.savetxt(out_csvfile, array2write, fmt='%.3f', delimiter=",", header="UTM_X, UTM_Y, CumulativeDistance_m, ChannelWidth_m, Sinuosity")
                  out_shapefile = '.'.join(file.split('.')[0:-1]) + "_centerline_%d.shp" %(smoothing_interval)
                  write_to_shapefile(OUT_CENTERLINE_SHAPE_PATH, out_shapefile, array2write, outSpatialRef, shp_layer='Centerline')
                  
                  #channel_width_df = pd.DataFrame(data={'Cumulative_distance_m':cell_distance_cs_m, 'channel_width_m':channel_width_m, 'sinuosity': sinuosity, 'dwnstream_idx':cell_downstream_list1D})
                  #pd.rolling_mean(channel_width_df['Cumulative_distance_m'], 10)
                  

                  #finding point with lowest channel width as starting point
                  #it is better to identify direction of flow from linear regression of cumulative distance vs. channel width
#                  idx_start_x, idx_start_y = np.where(distance_from_centerline_to_boundary.data == np.min(distance_from_centerline_to_boundary))
#                  distance_from_start_to_centerline=[]
#                  for i in range(len(idx_start_x)):
#                      distance_from_start_to_centerline.append(np.sqrt( np.power(idx_start_x[i] - endpoints_idxx[0],2) +  np.power(idx_start_y[i] - endpoints_idxy[0], 2) ))
#                  idx_start_point = np.where(np.min(distance_from_start_to_centerline) == np.array(distance_from_start_to_centerline))[0]
#                  idx_start_point_x = idx_start_x[idx_start_point]
#                  idx_start_point_y = idx_start_y[idx_start_point]
#                  #now compare the values from idx_start_point_x and idx_start_point_y with endpoints
#                  if np.min(np.abs(idx_start_point_x-endpoints_idxx)) == 0 and np.min(np.abs(idx_start_point_y-endpoints_idxy)) == 0:
#                      print('Difference between endpoint and min. line distance. X: %d, Y: %d' %(np.min(np.abs(idx_start_point_x-endpoints_idxx)), np.min(np.abs(idx_start_point_y-endpoints_idxy))))
#                  else:
#                      print('Ending distances are not the same: X: %d, Y: %d' %(np.min(np.abs(idx_start_point_x-endpoints_idxx)), np.min(np.abs(idx_start_point_y-endpoints_idxy))))
#                      print('Setting to stored endpoints.')
#                      idx = np.where(np.min(np.abs(idx_start_point_x-endpoints_idxx)))[0]
#                      idx_start_point_x = endpoints_idxx[idx]
#                      idx_start_point_y = endpoints_idxy[idx]
                  


def write_to_shapefile(OUT_CENTERLINE_SHAPE_PATH, fname, array2write, outSpatialRef, shp_layer='Centerline'):
    fieldname = '.'.join(fname.split('.')[0:-1])
    shapefile_fname = os.path.join(OUT_CENTERLINE_SHAPE_PATH, fname)
    driver = ogr.GetDriverByName("ESRI Shapefile")

    # create the data source
    data_source = driver.CreateDataSource(shapefile_fname)
    
    # create the spatial reference, WGS84
    srs = outSpatialRef
    
    # create the layer
    layer = data_source.CreateLayer(shp_layer, srs, ogr.wkbPoint)
    
    # Add the fields we're interested in
    layer.CreateField(ogr.FieldDefn("UTM_X", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("UTM_Y", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("CDist_m", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("CWidth_m", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("Sinuosity", ogr.OFTReal))
    field_region = ogr.FieldDefn("Filename", ogr.OFTString)
    field_region.SetWidth(24)
    layer.CreateField(field_region)
    
    # Process the text file and add the attributes and features to the shapefile
    for row in array2write:
      # create the feature
      feature = ogr.Feature(layer.GetLayerDefn())
      # Set the attributes using the values from the delimited text file
      feature.SetField("UTM_X", row[0])
      feature.SetField("UTM_Y", row[1])
      feature.SetField("CDist_m", row[2])
      feature.SetField("CWidth_m", row[3])
      feature.SetField("Sinuosity", row[4])
      feature.SetField("Filename", fieldname)
      # create the WKT for the feature using Python string formatting
      wkt = "POINT(%f %f)" %  (float(row[0]) , float(row[1]))
    
      # Create the point from the Well Known Txt
      point = ogr.CreateGeometryFromWkt(wkt)
    
      # Set the feature geometry using the point
      feature.SetGeometry(point)
      # Create the feature in the layer (shapefile)
      layer.CreateFeature(feature)
      # Dereference the feature
      feature = None
    
    # Save and close the data source
    data_source = None

def trace_path_distance(skeleton, idx_start_point_x, idx_start_point_y, idx_end_point_x, idx_end_point_y):
    #finding path along centerline and store distance
    #walk from N, NE, E, SE, S, SE, E, NW - reversed: E
    nn_distance = np.array([1, np.sqrt(2), 1, np.sqrt(2), 1, np.sqrt(2), 1, np.sqrt(2)])
    idxx = idx_start_point_x
    idxy = idx_start_point_y
    nr_of_elements = len(np.where(skeleton)[0])
    cell_distance = np.ones(nr_of_elements, dtype=np.float)
    euclidean_distance = np.zeros(nr_of_elements, dtype=np.float)
    cell_downstream_list2D = np.zeros((nr_of_elements, 2), dtype=np.uint16) 
    cell_downstream_list1D = np.zeros((nr_of_elements), dtype=np.uint16)
    
    #last next_cell_list will be 0, no cell is 'next cell'
    for i in range(nr_of_elements-1):
        cell_neighborhood = [skeleton[idxx, idxy+1], skeleton[idxx+1, idxy+1], skeleton[idxx+1, idxy], skeleton[idxx+1, idxy-1], 
                             skeleton[idxx, idxy-1], skeleton[idxx-1, idxy-1], skeleton[idxx-1, idxy], skeleton[idxx-1, idxy+1]]
      
        nn = np.where(cell_neighborhood)[0]
        #remove previous nn from cell
        if i > 0:
            nn = nn[nn != last_nnr]
        
        if len(nn) > 1:
            nn = nn[np.where(np.min(nn_distance[nn]))[0]]
            #store distance
        cell_distance[i] = nn_distance[nn]
        euclidean_distance[i] = np.sqrt(np.sum([np.power(idx_start_point_x-idxx, 2), np.power(idx_start_point_y-idxy, 2)] ))
        #print('i: %d, nn: %d, eucl. distance: %4.2f, idxx: %d, idxy: %d' %(i+1, nn, euclidean_distance[i], idxx, idxy))
        if nn == 0:
            idxy = idxy + 1
        elif nn == 1:
            idxx = idxx + 1
            idxy = idxy + 1
        elif nn == 2:
            idxx = idxx + 1
        elif nn == 3:
            idxx = idxx + 1
            idxy = idxy - 1
        elif nn == 4:
            idxy = idxy - 1
        elif nn == 5:
            idxx = idxx - 1
            idxy = idxy - 1
        elif nn == 6:
            idxx = idxx - 1
        elif nn == 7:
            idxx = idxx - 1
            idxy = idxy + 1
        cell_downstream_list2D[i] = [idxx, idxy]
        cell_downstream_list1D[i] = np.ravel_multi_index(cell_downstream_list2D[i], skeleton.shape)
        if nn > 3:
            last_nnr = nn - 4
        elif nn < 4:
            last_nnr = nn + 4
    
    cell_downstream_list2D[-1] = [idx_end_point_x, idx_end_point_y]
    cell_downstream_list1D[-1] = np.ravel_multi_index(cell_downstream_list2D[-1], skeleton.shape)
    cell_distance_cs_m = np.cumsum(cell_distance)
    #copy into centerline
    cell_distance_grid = np.zeros_like(skeleton, dtype=float)
    cell_distance_grid.ravel()[np.where(skeleton.ravel())[0]] = cell_distance_cs_m
    euclidean_distance[-1] = np.sqrt(np.sum([np.power(idx_start_point_x-idxx, 2), np.power(idx_start_point_y-idxy, 2)] ))
    sinuosity = cell_distance_cs_m / euclidean_distance
    sinuosity[0] = 1
    return cell_distance_cs_m, cell_downstream_list2D, cell_downstream_list1D, sinuosity, cell_distance_grid
                  
                                   
def bw_remove_3nn(bw):
    nn_8 = generate_binary_structure(2, 2)
    bw_filtered = convolve(bw.astype(np.uint8), nn_8, mode='nearest')
    bw_filtered = np.ma.masked_array(bw_filtered, ~bw)
    nngt4_idxx, nngt4_idxy = np.where(bw_filtered > 4)
    nneq4_idxx, nneq4_idxy = np.where(bw_filtered == 4)
    bw_filtered[nneq4_idxx, nneq4_idxy] = 0
    bw_indices_idxx, bw_indices_idxy = np.indices(bw_filtered.shape)
    #now remove neighboring pixels with the fewest neighbor pixels
    
    for i in range(len(nngt4_idxx)):
        #bw_nn = bw[nngt4_idxx[i]-1:nngt4_idxx[i]+2, nngt4_idxy[i]-1:nngt4_idxy[i]+2].ravel()
        bw_filtered_nn = bw_filtered[nngt4_idxx[i]-1:nngt4_idxx[i]+2, nngt4_idxy[i]-1:nngt4_idxy[i]+2].ravel()
        bw_indices_idxx_nn = bw_indices_idxx[nngt4_idxx[i]-1:nngt4_idxx[i]+2, nngt4_idxy[i]-1:nngt4_idxy[i]+2].ravel()
        bw_indices_idxy_nn = bw_indices_idxy[nngt4_idxx[i]-1:nngt4_idxx[i]+2, nngt4_idxy[i]-1:nngt4_idxy[i]+2].ravel()
        #mask outside river pixels of input image
        #bw_filtered_nn[np.where(bw_nn==False)] = 0
        
        dangling_idxs = np.where(np.min(bw_filtered_nn[bw_filtered_nn > 0]) == bw_filtered_nn)[0]
        for j in range(len(dangling_idxs)):
            #print('Removing dangling arcs. %d of %d' %(j+1, len(dangling_idxs)))
            #check neighborhood of both indices
            current_idxx = bw_indices_idxx_nn[dangling_idxs[j]]
            current_idxy = bw_indices_idxy_nn[dangling_idxs[j]]
            #bw_nn2 = bw[current_idxx-1:current_idxx+2, current_idxy-1:current_idxy+2].ravel()
            bw_filtered_nn2 = bw_filtered[current_idxx-1:current_idxx+2, current_idxy-1:current_idxy+2].ravel()
            #verify numbers of links in bw_filtered_nn2: IF larger than 1, there is connected neighbor
            if len(np.ma.where(bw_filtered_nn2 > 2)[0]) > 2:
                #contains connected neighbors, no dangling arc
                #print('connected neighbours are present')
                continue
            else:
                #print('removing: %d, %d' %(current_idxx, current_idxy))
                #this is the dangling arc to be removed
                bw_filtered[current_idxx, current_idxy] = 0
                bw[current_idxx, current_idxy] = 0
        
    #bw[nngt4_idxx, nngt4_idxy] = False
    return bw

def bwmorphDiag(bw):
    # filter for 8-connectivity of the background
    f = np.array(([1, -1, 0],[-1, 1, 0],[0, 0, 0]),dtype = np.int)
    # initialize result with original image
    bw = bw.astype(np.int)
    res2 = bw.copy().astype(np.bool)
    for ii in range(4): # all orientations
        # add results where sum equals 2 -> two background pixels on the
        # diagonal with 2 foreground pixels on the crossing mini-anti-diagonal
        res2 = res2 | (convolve(np.invert(bw),f) == 2)
        f = np.rot90(f) # rotate filter to next orientation
    return res2

def find_network_endpoints(skel, endpoint_kernel='D8'):    
    #this will work only on int() images, which will be too large for this purpose
    # use thinned image to find endpoints, using scipy.ndimage.convolce. cv2.filter2D is another option
    #modified after: https://stackoverflow.com/questions/26537313/how-can-i-find-endpoints-of-binary-skeleton-image-in-opencv
    if endpoint_kernel == 'D8':
        endpoint_kernel = np.uint8([[1, 1, 1],[1,10,1],[1,1,1]])
    if endpoint_kernel == 'D4':
        endpoint_kernel = np.uint8([[0, 1, 0],[1,10,1],[0,1,0]])
    skel_filtered = convolve(skel.astype(np.uint8), endpoint_kernel, mode='nearest')

    endpoints_idxx, endpoints_idxy = np.where(skel_filtered == 11)
    endpoints = np.zeros_like(skel)
    endpoints[np.where(skel_filtered==11)] = 1
    endpoints_idx = np.where(skel_filtered.ravel() == 11)[0]
    return endpoints, endpoints_idx, endpoints_idxx, endpoints_idxy

def find_centerline_MATLAB(OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, CODE_PATH, matlab_command, OUT_CENTERLINE_SHAPE_PATH,
                    gdal_translate_command, outSpatialRef):
    ### start Matlab script for each TIF file in OUT_TIF_PATH
    cwd = os.getcwd()
    os.chdir(CODE_PATH)
    for file in os.listdir(OUT_TIF_PATH):
         if file.endswith('.tif'):
              in_tiffile = os.path.join(OUT_TIF_PATH, file)
              print("Finding Centerline %s" %in_tiffile)
    
              #verify if OUT_CENTERLINE_TIF_PATH exists
              if not os.path.exists(OUT_CENTERLINE_TIF_PATH):
                  os.makedirs(OUT_CENTERLINE_TIF_PATH)
              if os.path.exists(OUT_CENTERLINE_TIF_PATH) == False:
                  print("Can not create directory: %s" % OUT_CENTERLINE_TIF_PATH)
    
              out_tiffile1 = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + "_centerline.tif")
              out_tiffile2 = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + "_centerline_cumdist.tif")
    
              if os.path.exists(out_tiffile1) and os.path.exists(out_tiffile2):   
                  print("%s  and %s exist, skipping to next file..." %(out_tiffile1.split('/')[-1], out_tiffile2.split('/')[-1]))
              else:    
                  #call Matlab to create centerline and cumulative distance TIFs
                  subprocess.call([matlab_command+" -nosplash -nodisplay -r \"chanextract(\'%s\',\'%s\',0)\"" % (in_tiffile, out_tiffile1)],shell=True);
    
    os.chdir(cwd)

    ### post-process channel width from Matlab outputs
    for file in os.listdir(OUT_CENTERLINE_TIF_PATH):
         if file.endswith('centerline.tif'):
              in_centerline_tiffile = os.path.join(OUT_CENTERLINE_TIF_PATH, file)
              print("Vectorizing results: %s" % in_centerline_tiffile)
              cumdist_filename = file.split('.')[0] + '_cumdist.tif'
              in_cumdist_tiffile = os.path.join(OUT_CENTERLINE_TIF_PATH, cumdist_filename)
              
              #verify if OUT_CENTERLINE_SHAPE_PATH exists
              if not os.path.exists(OUT_CENTERLINE_SHAPE_PATH):
                  os.makedirs(OUT_CENTERLINE_SHAPE_PATH)
              if os.path.exists(OUT_CENTERLINE_SHAPE_PATH) == False:
                  print("Can not create directory: %s" % OUT_CENTERLINE_SHAPE_PATH)
    
              out_centerline_tiffile = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + "2.tif")
              out_cumdist_tiffile = os.path.join(OUT_CENTERLINE_TIF_PATH, cumdist_filename.split('.')[0] + "2.tif")
              
              #csv file contains RWidth_m and CumDist_m
              out_csvfile = os.path.join(OUT_CENTERLINE_TIF_PATH, file.split('.')[0] + ".csv")
              out_shapefile = os.path.join(OUT_CENTERLINE_SHAPE_PATH, file.split('.')[0] + ".shp")
    
              if os.path.exists(out_shapefile):   
                  print("%s exists, skipping to next file..." %out_shapefile)
              else:    
                  #call gdal to set background value (0 from Matlab-tif files)
                  out_tiffile_txt = out_centerline_tiffile.split('.')[0] + '.gdal_translate.out'
                  gdal_translate_subprocess_command = gdal_translate_command  + ' -a_nodata 0 ' + ' -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 ' + ' -a_nodata 0 ' + in_centerline_tiffile + ' ' + out_centerline_tiffile + '>'  + out_tiffile_txt
                  os.system(gdal_translate_subprocess_command)
                  gdal_translate_subprocess_command = gdal_translate_command + ' -a_nodata 0 ' + ' -co COMPRESS=DEFLATE -co PREDICTOR=2 -co ZLEVEL=9 ' + in_cumdist_tiffile + ' ' + out_cumdist_tiffile + '>>'  + out_tiffile_txt
                  os.system(gdal_translate_subprocess_command)
                  #delete centerline and cumdist files with NoDATA and rename *2.tif to previous filename
                  os.remove(in_centerline_tiffile)
                  os.rename(out_centerline_tiffile, in_centerline_tiffile)
                  os.remove(in_cumdist_tiffile)
                  os.rename(out_cumdist_tiffile, in_cumdist_tiffile)
                  
                  #convert tif file to CSV file (use only values > 0)
                  raster2point(in_centerline_tiffile, in_cumdist_tiffile, out_csvfile)              
    
                  #convert to point shapefile
                  csv2shapefile( out_csvfile, out_shapefile, outSpatialRef )
              #os.remove(in_centerline_tiffile)
              #os.remove(in_cumdist_tiffile)

def merge_CTRL_shapefiles(shapefile_merged_out, shapemerger_command, IN_PATH):     
    ### merge all channels into one shapefile
    #use shapemerger.py from toni @ http://gistncase.blogspot.com
    if os.path.exists(shapefile_merged_out):   
        print('Merged Shapefile %s exists.' %(shapefile_merged_out))
    else:    
        #call shapemerger to merge all centerline shapefiles
        cwd = os.getcwd()
        os.chdir(IN_PATH)
        print('Using %s' %(shapemerger_command))
        subprocess.call(['python2 %s -o %s -e -2 Centerline_Shapefiles/*.shp' %(shapemerger_command, shapefile_merged_out)], shell=True)
        os.chdir(cwd)
# Call matlab <label>__DEM_to_FAC_STR.m
#ADD projection information to output from Matlab
#ogr2ogr -s_srs epsg:32644 -t_srs epsg:32644 WHimalaya_srtm1_1km2_UTM44N_Wgs84.shp WHimalaya_srtm1_1km2.shp
              
def postprocess_shapefiles(IN_PATH, SRTM1_shp, CtrLine_shp, distance_threshold, search_radius, datapath2save, 
                           path2save, overview_fig_fn, hdf_fname, ascii_fname):
    ### Postprocess merged shapefile
    print('Postprocessing Shapefiles: %s and %s' %(os.path.basename(SRTM1_shp), os.path.basename(CtrLine_shp)) )
    cwd = os.getcwd()
    os.chdir(IN_PATH)
    #get data from SRTM file
    SRTM1_UTM_x, SRTM1_UTM_y, SRTM1_elevation, SRTM1_d_maxdd, SRTM1_DA_km2, SRTM1_gradient, SRTM1_ksn, SRTM1_distance = load_SRTM_shp_data(SRTM1_shp)
                 
    #get data from high-res Riverwidth files
    CtrLine_UTM_x, CtrLine_UTM_y, RWidth_m, CumDist_m, Shp_fname, Shp_src_ID = load_Ctrl_shp_data(CtrLine_shp)
    
    #process and sort RWidth data, store as Pandas
    RWidth_m = np.array(RWidth_m)
    CumDist_m = np.array(CumDist_m)
    Shp_fname_list = np.unique(Shp_fname)
    Shp_fname_list_nr = len(Shp_fname_list)
    Shp_fname_id = np.zeros(len(Shp_fname_list))
    Shp_fname_idx = []
    Shp_fname_txt=[]
    col_labels = ['CtrLine_UTM_x', 'CtrLine_UTM_y', 'RWidth_m', 'CumDist_m',
                  'Shp_src_ID']
    RWidth_data_sorted = []    
    for i in range(Shp_fname_list_nr):
        idx_fname = np.where(np.array(Shp_fname) == Shp_fname_list[i])[0]
        Shp_fname_txt.append(np.array(Shp_fname)[idx_fname][0])
        Shp_fname_id[i] = i
        Shp_fname_idx.append(idx_fname)
        c_river = np.array([np.array(CtrLine_UTM_x)[idx_fname],
                            np.array(CtrLine_UTM_y)[idx_fname],RWidth_m[idx_fname],
                            CumDist_m[idx_fname],
                            np.array(Shp_src_ID)[idx_fname]])
                                    
        #c_river_name=np.array(Shp_fname)[idx_fname]
        #generate array and then sort by cumulative river distance
        c_river_pd = pd.DataFrame(c_river.T, columns=col_labels)
        #you can add fitting and further filtering here
        RWidth_data_sorted.append(c_river_pd.sort_values(by='CumDist_m'))
           
    CtrlLine_pts_x = np.array(CtrLine_UTM_x)
    CtrlLine_pts_y = np.array(CtrLine_UTM_y)
    CtrlLine_pts = np.vstack((CtrlLine_pts_x, CtrlLine_pts_y)).T
    
    SRTM1_pts_x = np.array(SRTM1_UTM_x)
    SRTM1_pts_y = np.array(SRTM1_UTM_y)
    SRTM1_pts = np.vstack((SRTM1_pts_x, SRTM1_pts_y)).T
    
    #generate overview map
    if os.path.exists(overview_fig_fn) == False:
        plt.clf()
        plt.scatter(CtrlLine_pts_x, CtrlLine_pts_y, s=5, c='r', marker='+', label='Centerline')
        plt.scatter(SRTM1_pts_x, SRTM1_pts_y, s=5, c='k', marker='.', label='SRTM-FAC')
        plt.legend()
        f = plt.gcf()  # f = figure(n) if you know the figure number
        f.set_size_inches(11.69,8.27)
        plt.savefig(overview_fig_fn, papertype='a4', orientation='landscape')
    
    #Generating Kd-Tree
    print('Generating Kd-Tree from %s' %os.path.basename(CtrLine_shp))
    T = KDTree(CtrlLine_pts)
    distance, idx = T.query(SRTM1_pts,k=1)
    
    SRTM1_merged_data=[]
    Ctrline_merged_text=[]
    for i in range(len(distance)):
        if distance[i] < distance_threshold:
            #calculate mean of riverwidth for each point
            Riverwidth_c = RWidth_m[idx[i]]
            Riverwidth_c_nh = T.query_ball_point(CtrlLine_pts[idx[i]],r=search_radius)
            SRTM1_merged_data.append([SRTM1_UTM_x[i], SRTM1_UTM_y[i], SRTM1_elevation[i], 
                                      SRTM1_gradient[i], SRTM1_distance[i], SRTM1_d_maxdd[i],
                                 RWidth_m[idx[i]], np.mean(RWidth_m[Riverwidth_c_nh]), 
                                 np.std(RWidth_m[Riverwidth_c_nh]), np.median(RWidth_m[Riverwidth_c_nh]),
                                 SRTM1_DA_km2[i], CumDist_m[idx[i]], distance[i]])
            Ctrline_merged_text.append(Shp_fname[idx[i]])
    
    #SRTM1_merged_data_pd = pd.DataFrame(np.array(SRTM1_merged_data))
    SRTM1_merged_data = np.array(SRTM1_merged_data)
    
    #extract list of Centerline filenames:
    Ctrline_fname_list = np.unique(Ctrline_merged_text)
    Ctrline_fname_list_nr = len(Ctrline_fname_list)
    Ctrline_fname_id = np.zeros(len(Ctrline_merged_text))
    Ctrline_fname_idx = []
    col_labels = ['SRTM1_UTM_x', 'SRTM1_UTM_y', 'SRTM1_elevation', 'SRTM1_gradient',
                  'SRTM1_distance_from_outlet', 'SRTM1_maxdistance_from_channelhead',
                  'RWidth_m_closest', 'RWidth_m_mean', 'RWidth_m_std', 
                  'RWidth_m_median', 'SRTM1_DA_km2', 'RWidth_CumDist_m', 
                  'SRTM1_Ctrline_distance_m']
    
    SRTM1_merged_data_sorted = []    
    for i in range(Ctrline_fname_list_nr):
        idx_fname = np.where(np.array(Ctrline_merged_text) == Ctrline_fname_list[i])[0]
        #np.array(Ctrline_merged_text)[idx_fname]
        Ctrline_fname_id[idx_fname] = i
        Ctrline_fname_idx.append(idx_fname)
        c_river = np.array(SRTM1_merged_data)[idx_fname]
        #generate array and then sort by cumulative river distance
        c_river_pd = pd.DataFrame(c_river, columns=col_labels)
        #you can add fitting and further filtering here
        SRTM1_merged_data_sorted.append(c_river_pd.sort_values(by='RWidth_CumDist_m'))
        
    SRTM1_merged_data_sorted_df = pd.concat(SRTM1_merged_data_sorted)
    
    for i in range(Ctrline_fname_list_nr):
        mk_LRP_DA_subplots(SRTM1_merged_data_sorted, RWidth_data_sorted, 
                           Ctrline_fname_list, Shp_fname_txt, path2save, i)
    
    mk_single_plots(SRTM1_merged_data_sorted, SRTM1_merged_data_sorted_df, Ctrline_fname_list, Ctrline_fname_list_nr )
    
    ### Save to HDF
    SRTM1_merged_data_sorted_df.to_hdf(hdf_fname, hdf_fname.split('.')[0], mode='w')
    SRTM1_merged_data_sorted_df.to_csv(ascii_fname)
    
    if os.path.exists(datapath2save) == False:
        os.mkdir(datapath2save)
        
    #save to individual data files
    for i in range(len(SRTM1_merged_data_sorted)):
        csv2save = os.path.join(datapath2save, Ctrline_fname_list[i] + '.ascii')
        SRTM1_merged_data_sorted[i].to_csv(csv2save)

    os.chdir(cwd)
