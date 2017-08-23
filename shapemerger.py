#!/usr/bin/python
# shapemerger.py
# v.1.2 by toni @ http://gistncase.blogspot.com
#
# shapefile merger utility 

import os, sys, glob
try:
    from osgeo import ogr, gdal
except:
    try:
        import ogr, gdal
    except:
        print 'Error: shapemerge.py cannot find python OGR and GDAL modules'
        print '       for handling geometries. OsGeo4W python is recommended.'
        sys.exit(1)

def Usage():
    print
    print 'Shapefile merger'
    print 'Usage: shapemerger.py\t[-o <output.shp>] [-n] [-i] [-s <fieldname>] [-f <fieldname>]'
    print '\t\t\t[-e] [-2] [-r] [-t polygon|linestring|point]'
    print '\t\t\t<shapefiles> [-x <exclude shapefile list>]'
    print 
    print '  <shapefiles>: input shapefile list (accepts wildcards)'
    print
    print '  -o:\toutput shapefile name (default: shapemerge.shp)'
    print '  -n:\tdo not add fields for the source filenames and fids'
    print '  -i:\tdbf schemes intersection (default: union)'
    print '  -s:\tcustom name for the source filenames field (default: SOURCESHP)'
    print '  -f:\tcustom name for the source fids field (default: SOURCEFID)'
    print '  -e:\tuse existing source filenames and fids fields, if any'
    print '  -2:\tforce output geometry type to 2D'
    print '  -x:\tlist of shapefiles to exclude (accepts wildcards)'
    print '  -r:\trecursive (grabs every path of input files)'
    print '  -t:\tfilter input shapefiles by geometry type'
    sys.exit(1)

strCurrentPath = os.path.abspath('')
os.chdir(strCurrentPath)

strOutDirname = ''
strMergeName = 'shapemerge'
strOutBasename = 'shapemerge.shp'

pOgrShpDriver = ogr.GetDriverByName('Esri Shapefile')

lsFilenames=[]
lsExcludeFilenames=[]
lsVerifiedFilenames=[]
lsMergeFieldsToRetain=[]
lsMergeFieldNames=[]
lsMergeFieldOgrDefns=[]

dGeomTypeFilters={'POLYGON': ogr.wkbPolygon, 'AREAL': ogr.wkbPolygon, 'LINESTRING': ogr.wkbLineString, 'LINE': ogr.wkbLineString, 'POINT': ogr.wkbPoint, 'POINTS': ogr.wkbPoint}

bFlagExclude = 0
bFlagAddSourcefields = 1
bFlagUseExistingSfields = 0
bFlagIntersection = 0
bFlagRecursive = 0
bFlagGeomTFilter = 0
bFlag2D = 0

bGeomTypeRefIsSet = 0

strFilenameFieldname = 'SOURCESHP'
strFidFieldname = 'SOURCEFID'

i = 1

# parse commandline (first check for recursive flag)
while i < len(sys.argv): 
    arg = sys.argv[i]
    if arg.upper() == '-R':
        bFlagRecursive = 1
        del sys.argv[i]
        break
    i = i + 1

i = 1
    
while i < len(sys.argv): 
    arg = sys.argv[i]

    if arg.upper() == '-N':
        bFlagAddSourcefields = 0
        bFlagExclude = 0

    elif arg.upper() == '-I':
        try:
            strGdalVersion = int(gdal.VersionInfo()[:2])
        except:
            strGdalVersion = 0
        
        if (strGdalVersion < 19):
            print 'Warning: dbf schemes intersection supported only with gdal >= 1.9.0, ignoring.'
            bFlagIntersection = 0
        else:
            bFlagIntersection = 1
        
        bFlagExclude = 0

    elif arg.upper() == '-S':
        i = i + 1
        strFilenameFieldname = ''.join(sys.argv[i].split())[:10].upper() #upper, tronca al decimo, toglie spazi
        bFlagExclude = 0

    elif arg.upper() == '-F':
        i = i + 1
        strFidFieldname = ''.join(sys.argv[i].split())[:10].upper() #upper, tronca al decimo, toglie spazi
        bFlagExclude = 0
        
    elif arg.upper() == '-E':
        bFlagUseExistingSfields = 1
        bFlagExclude = 0

    elif arg.upper() == '-2':
        bFlag2D = 1
        bFlagExclude = 0

    elif arg.upper() == '-T':
        i = i + 1
        strGeomTFilter = (sys.argv[i]).upper()
        
        if (not (strGeomTFilter in dGeomTypeFilters)):
            print "Warning: invalid geometry type declared, ignoring." 
        else:
            intGeometryTypeRef = dGeomTypeFilters[strGeomTFilter]
            bFlagGeomTFilter = 1
        
        bFlagExclude = 0
        
    elif arg.upper() == '-X':
        bFlagExclude = 1

    elif arg.upper() == '-O':
        i = i + 1
        bFlagExclude = 0

        strOutDirname, strOutBasename = os.path.split(sys.argv[i])
        strMergeName, strExtension = os.path.splitext(strOutBasename)
            
    else:
        # expands wildcards: exclude, input list/recursive input list
        if (bFlagRecursive):
            strInputShapePath, strInputShapeBasename = os.path.split(arg)
            strInputShapePath = strInputShapePath or '.'
            
            lsInputShapeRoots=[]
            lsInputShapeRoots.extend(glob.glob(strInputShapePath)) #expands paths with wildcards

        if (bFlagExclude and bFlagRecursive):
            for strRootPath in lsInputShapeRoots:
                for strRoot,strDirs,strFiles in os.walk(strRootPath):
                    lsExcludeFilenames.extend ( glob.glob(strRoot + os.sep + strInputShapeBasename) )
        
        elif (bFlagExclude and (not bFlagRecursive)):
            lsExcludeFilenames.extend ( glob.glob(arg) )
        
        elif (bFlagRecursive):
            for strRootPath in lsInputShapeRoots:
                for strRoot,strDirs,strFiles in os.walk(strRootPath):
                    lsFilenames.extend ( glob.glob(strRoot + os.sep + strInputShapeBasename) )
        else:
            lsFilenames.extend ( glob.glob(arg) )
                
    i = i + 1

for strItem in lsExcludeFilenames:
    try:
        lsFilenames.remove(strItem)
    except:
        pass

if len(lsFilenames) == 0:
    Usage()
    sys.exit(1)

if not os.path.isdir(strOutDirname or '.'):
    print strOutDirname + ': invalid directory, outputting here'
    strOutDirname = ''

if os.path.isfile( os.path.join(strOutDirname, strOutBasename) ):
    strSuffix = ''
    j = 0
    while os.path.isfile( os.path.join(strOutDirname,strMergeName + strSuffix + '.shp') ):
        j = j + 1
        strSuffix = '_' + str(j)
  
    if (strMergeName != 'shapemerge'):
        print strOutBasename + ': file exists, outputting in ' + strMergeName + strSuffix + '.shp'
  
    if (j > 0):
        strMergeName = strMergeName + strSuffix
        strOutBasename = strMergeName + '.shp'  
  
# geometry type check / "source attributes" renaming
strSuffix = ''
j = 0 
for i,strFile in enumerate(lsFilenames):
    try:
        pOgrDatasource = pOgrShpDriver.Open(strFile)
        pOgrLayer = pOgrDatasource.GetLayer()
        pOgrLayerSchema = pOgrLayer.GetLayerDefn()
        
        if ( (not bGeomTypeRefIsSet) and (not bFlagGeomTFilter) ):

            # picks up the first geometry type as reference if geometry type filter is not set
            intGeometryTypeRef = pOgrLayer.GetGeomType()
            
            if (bFlag2D):
                # fetch first feature to flatten geometry and grab TypeRef
                pOgrFeature = pOgrLayer.GetNextFeature()
                pOgrFeatureGeometry = pOgrFeature.GetGeometryRef()
                pOgrClonedGeometry = pOgrFeatureGeometry.Clone()
                pOgrClonedGeometry.FlattenTo2D()
                intGeometryTypeRef = pOgrClonedGeometry.GetGeometryType() 
                #pOgrFeature.Destroy()
            
            bGeomTypeRefIsSet = i # register index of the reference shapefile
                
        #check if layers the picked geometry type (ignore 3D flag)
        if pOgrLayer.GetGeomType() != intGeometryTypeRef:
            if ((pOgrLayer.GetGeomType() == ogr.wkbPolygon) or (pOgrLayer.GetGeomType() == ogr.wkbPolygon25D)) and ((intGeometryTypeRef == ogr.wkbPolygon) or (intGeometryTypeRef == ogr.wkbPolygon25D)) \
                or ((pOgrLayer.GetGeomType() == ogr.wkbLineString) or (pOgrLayer.GetGeomType() == ogr.wkbLineString25D)) and ((intGeometryTypeRef == ogr.wkbLineString) or (intGeometryTypeRef == ogr.wkbLineString25D)) \
                or ((pOgrLayer.GetGeomType() == ogr.wkbPoint) or (pOgrLayer.GetGeomType() == ogr.wkbPoint25D)) and ((intGeometryTypeRef == ogr.wkbPoint) or (intGeometryTypeRef == ogr.wkbPoint25D)) :
                    pass
            else:
                if (not bFlagGeomTFilter):
                    #print 'Error: shapefiles must have the same geometry type.'
                    #print '\t' + os.path.basename(lsFilenames[0]) + ': ' + ogr.GeometryTypeToName(intGeometryTypeRef).upper()
                    #print '\t' + os.path.basename(lsFilenames[i]) + ': ' + ogr.GeometryTypeToName(pOgrLayer.GetGeomType()).upper() 
                    #sys.exit(1)
                    print 'Warning: ' + os.path.basename( strFile ) + ' type ' + ogr.GeometryTypeToName(pOgrLayer.GetGeomType()).upper() + ' unfits first file type ' + ogr.GeometryTypeToName(intGeometryTypeRef).upper() + ' (' + os.path.basename(lsFilenames[bGeomTypeRefIsSet]) +'), skipping.'
                continue
                
        if ( (not bFlagUseExistingSfields) and bFlagAddSourcefields ):
            # check for "source fields" in current shapefile to eventually rename them
            while ( (pOgrLayerSchema.GetFieldIndex(strFilenameFieldname) != -1) or (pOgrLayerSchema.GetFieldIndex(strFidFieldname) != -1) ):
                j = j + 1
                strSuffix = '_' + str(j)
                strFilenameFieldname = strFilenameFieldname[:(10 - len(strSuffix))] + strSuffix
                strFidFieldname      = strFidFieldname[:(10 - len(strSuffix))] + strSuffix

        lsVerifiedFilenames.append( strFile )
        pOgrDatasource.Destroy()
        
    except:
        print 'Warning: ' + strFile + ' invalid shapefile, skipping.'
        #continue

if (len(lsVerifiedFilenames) == 0):
    print 'Error: all the input shapefiles are invalid.'
    sys.exit(1)
else:
    lsFilenames = lsVerifiedFilenames

if (j >0):
    print "Warning: \"source fieldnames\" renamed to: " +  strFilenameFieldname + ", " + strFidFieldname

# building attributes
if (bFlagAddSourcefields):
    # "source fields" attributes
    lsMergeFieldOgrDefns.append( ogr.FieldDefn(strFilenameFieldname, ogr.OFTString) )
    lsMergeFieldOgrDefns[0].SetWidth(255)
    lsMergeFieldOgrDefns.append( ogr.FieldDefn(strFidFieldname, ogr.OFTInteger) )
    lsMergeFieldOgrDefns[1].SetWidth(9) #max one billion features! to be fixed!
    
    lsMergeFieldNames.append(strFilenameFieldname)
    lsMergeFieldNames.append(strFidFieldname)


for i,strFile in enumerate(lsFilenames):
    # shapefiles attributes
    
    if ( i>0 and bFlagIntersection and ( (bFlagAddSourcefields and (len(lsMergeFieldNames)==2)) or (len(lsMergeFieldNames)==0) ) ):
        # no more scheme intersection to do, exit loop 
        break
    
    strFullname = strFile
    strDirname, strBasename = os.path.split(strFullname)
    strShapename, strExtension = os.path.splitext(strBasename)
    
    pOgrShpDatasource = pOgrShpDriver.Open (strFullname)
    pOgrShpLayer = pOgrShpDatasource.GetLayer()
    pOgrShpLayerSchema = pOgrShpLayer.GetLayerDefn()
     
    for j in range( pOgrShpLayerSchema.GetFieldCount() ):
        # parse shapefile attributes
        pOgrShpFielddefn = pOgrShpLayerSchema.GetFieldDefn(j)
        strShpFieldname = pOgrShpFielddefn.GetNameRef().upper()
        
        try:
            iMergeListsFieldIdx = lsMergeFieldNames.index(strShpFieldname)
        except:
            iMergeListsFieldIdx = -1
 
        if (iMergeListsFieldIdx == -1):
            # field doesn't exist: create it. If flag intersection is on, parse only the first file
            
            if ( (not bFlagIntersection) or (bFlagIntersection and (i == 0)) ):
                
                pOgrFielddefn = ogr.FieldDefn( strShpFieldname, pOgrShpFielddefn.GetType() )
                pOgrFielddefn.SetWidth( pOgrShpFielddefn.GetWidth() )
                pOgrFielddefn.SetPrecision( pOgrShpFielddefn.GetPrecision() )
                
                #ONLY TEXT:
                #pOgrFielddefn = ogr.FieldDefn( strShpFieldname, ogr.OFTString )
                #pOgrFielddefn.SetWidth( pOgrShpFielddefn.GetWidth() + pOgrShpFielddefn.GetPrecision() )
                
                lsMergeFieldOgrDefns.append( pOgrFielddefn )
                lsMergeFieldNames.append( strShpFieldname )
                
                # add field to the retain list for intersection
                lsMergeFieldsToRetain.append( len(lsMergeFieldNames) - 1 )
        
        else:
            # field exists, check properties
            pOgrFielddefn = lsMergeFieldOgrDefns[iMergeListsFieldIdx]
            
            #same data type
            if ( pOgrFielddefn.GetType() == pOgrShpFielddefn.GetType() ):
            
                if pOgrFielddefn.GetWidth() < pOgrShpFielddefn.GetWidth():
                    pOgrFielddefn.SetWidth( pOgrShpFielddefn.GetWidth() )
            
                if pOgrFielddefn.GetPrecision() < pOgrShpFielddefn.GetPrecision():
                    pOgrFielddefn.SetPrecision( pOgrShpFielddefn.GetPrecision() )
                    
            else:
                # different data type: "translate" it into string
                pOgrFielddefn.SetType( ogr.OFTString )
                pOgrFielddefn.SetWidth( pOgrFielddefn.GetWidth() + pOgrFielddefn.GetPrecision() )
                
                if ( pOgrFielddefn.GetWidth() < pOgrShpFielddefn.GetWidth() + pOgrShpFielddefn.GetPrecision() ): 
                    pOgrFielddefn.SetWidth( pOgrShpFielddefn.GetWidth() + pOgrShpFielddefn.GetPrecision() )
                                        
            if (bFlagIntersection):
                # add field to the retain list for intersection
                lsMergeFieldsToRetain.append(iMergeListsFieldIdx)
            
    if ( bFlagIntersection ):
        # performing dbf intersection
        
        if (bFlagAddSourcefields):
            lsMergeFieldsToRetain.append( lsMergeFieldNames.index(strFilenameFieldname) )
            lsMergeFieldsToRetain.append( lsMergeFieldNames.index(strFidFieldname) )
            
        lsMergeFieldsToDelete = range( len(lsMergeFieldNames)-1, -1, -1) # reverse: from last index to 0
        
        for iIdx in lsMergeFieldsToRetain:
            try:
                lsMergeFieldsToDelete.remove(iIdx)
            except:
                pass
       
        for iIdx in lsMergeFieldsToDelete:
            # add or remove fields items
            del lsMergeFieldOgrDefns[iIdx]
            del lsMergeFieldNames[iIdx]            

    pOgrShpDatasource.Destroy()
    lsMergeFieldsToRetain = []

# generate the merge layer
strMergeShpFullname = os.path.join(strOutDirname, strOutBasename)

pOgrMergeShpDatasource = pOgrShpDriver.CreateDataSource( strMergeShpFullname )
pOgrMergeShpLayer = pOgrMergeShpDatasource.CreateLayer( strMergeName, geom_type = intGeometryTypeRef )

for i in lsMergeFieldOgrDefns:
    # create the fields
    pOgrMergeShpLayer.CreateField(i,1)

pOgrMergeShpLayerSchema = pOgrMergeShpLayer.GetLayerDefn()

# Loop for populating the merge shapefile
for i,strFile in enumerate(lsFilenames):
    strFullname = strFile
    strDirname, strBasename = os.path.split(strFullname)
    strShapename, strExtension = os.path.splitext(strBasename)
    
    pOgrShpDatasource = pOgrShpDriver.Open (strFullname)
    pOgrShpLayer = pOgrShpDatasource.GetLayer()
    pOgrShpLayerSchema = pOgrShpLayer.GetLayerDefn()
            
    pOgrShpLayer.ResetReading()
    
    pOgrShpFeature = pOgrShpLayer.GetNextFeature()    
    
    bFlagExistsFilenameFld = pOgrShpLayerSchema.GetFieldIndex(strFilenameFieldname)
    bFlagExistsFidFld = pOgrShpLayerSchema.GetFieldIndex(strFidFieldname)
    
    while pOgrShpFeature:
        pOgrMergeShpfeature = ogr.Feature(pOgrMergeShpLayerSchema)
        pOgrMergeShpfeature.SetFrom(pOgrShpFeature,1)
        
        if (bFlag2D):
            # performing 2d conversion
            pOgrMergeGeometry2d = pOgrMergeShpfeature.GetGeometryRef().Clone()
            pOgrMergeGeometry2d.FlattenTo2D()
            pOgrMergeShpfeature.SetGeometryDirectly(pOgrMergeGeometry2d)
        
        if (bFlagAddSourcefields):
            # "source fields" values, if any and if they should be populated
            if ( (bFlagUseExistingSfields and bFlagExistsFilenameFld == -1) or (not bFlagUseExistingSfields) ):
                #pOgrMergeShpfeature.SetField( strFilenameFieldname, strBasename )
                pOgrMergeShpfeature.SetField( strFilenameFieldname, strShapename)
            if ( (bFlagUseExistingSfields and bFlagExistsFidFld == -1) or (not bFlagUseExistingSfields) ):
                pOgrMergeShpfeature.SetField( strFidFieldname, pOgrShpFeature.GetFID() )
            
        pOgrMergeShpLayer.CreateFeature(pOgrMergeShpfeature)
        
        pOgrMergeShpfeature.Destroy()
        pOgrShpFeature.Destroy()
        pOgrShpFeature = pOgrShpLayer.GetNextFeature()
            
    pOgrShpDatasource.Destroy()

# flush all to disk before quitting
pOgrMergeShpDatasource.SyncToDisk()
pOgrMergeShpDatasource.Destroy()

print
print strMergeShpFullname
sys.exit(0)
