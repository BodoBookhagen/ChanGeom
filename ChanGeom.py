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

class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self          

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
            band_data1 = Numeric.reshape( band_data1, (srcwin1[2],) )
            band_data2 = band2.ReadAsArray( srcwin1[0], y, srcwin1[2], 1 )    
            band_data2 = Numeric.reshape( band_data2, (srcwin1[2],) )
            
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
    reader = csv.DictReader(open(csvfile,"rb"),
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

def find_centerline(OUT_TIF_PATH, OUT_CENTERLINE_TIF_PATH, matlab_command, OUT_CENTERLINE_SHAPE_PATH,
                    gdal_translate_command, outSpatialRef):
    ### start Matlab script for each TIF file in OUT_TIF_PATH
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

def merge_shapefiles(shapefile_merged_out, shapemerger_command, IN_PATH):     
    ### merge all channels into one shapefile
    #use shapemerger.py from toni @ http://gistncase.blogspot.com
    if os.path.exists(shapefile_merged_out):   
        print('Merged Shapefile %s exists.' %(shapefile_merged_out))
    else:    
        #call shapemerger to merge all centerline shapefiles
        cwd = os.getcwd()
        os.chdir(IN_PATH)
        subprocess.call([shapemerger_command + ' -o %s -e -2 Centerline_Shapefiles/*.shp' %(shapefile_merged_out)],shell=True)
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
