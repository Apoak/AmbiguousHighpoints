import os
os.environ['PROJ_LIB'] = 'C:\\Users\\andre\OneDrive_CalPoly\\Documents\\SeniorProject\\Code\\pyGeo\\lib\\site-packages\\osgeo\\data\\proj'
# os.environ['GDAL_DATA'] = 'C:\\Users\\Sai kiran\\anaconda3\\envs\sai\\Library\\share'
import numpy as np

from multiprocessing import Pool
from osgeo import gdal, ogr, osr
from pathlib import Path

SINGLE_FILE = "LiDAR/Storey/USGS_1M_11_x27y435_NV_WestCentral_EarthMRI_2020_D20.tif"
TEST_DIR = Path("LiDAR/Storey")
OUT_DIR = "reprojected/storey/"
shp_path = "ShapeFiles/tl_2024_us_county.shp"
SHP_OUT = "ShapeFiles/out/"

def create_shapefile(boundary, path):
    """Create a shapefile from a boundary"""

    # path = SHP_OUT + "boundary_buff.shp"
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(path)

    if data_source is None:
        print("Failed to create file")
    else:
        print("Shapefile created successfully!")

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    
    layer = data_source.CreateLayer("boundary", srs, ogr.wkbPolygon)

    # Add an ID field
    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    layer.CreateField(id_field)

    # Create a new feature
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetField("id", 1)
    feature.SetGeometry(boundary)

    layer.CreateFeature(feature)

    # Save and close everything
    data_source = None
    return path

def reproject_shapefile(target, shapefile, flag, out_path=SHP_OUT + "boundary_reprojected2.shp"):
    # CONTAINS COUNTY FILTERING
    """"Changes the shapefile projection to match the tif file and then creates a new shapefile"""
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile, 0)  # 0 means read-only mode
    layer = data_source.GetLayer()

    # Print geometries and attributes of each feature
    if shapefile == shp_path:
        layer.SetAttributeFilter("NAMELSAD = 'Storey County' AND STATEFP = '32'")

    #set spatial reference and transformation
    sourceprj = layer.GetSpatialRef()

    # Target is a tif or shapefile
    if flag == "tif":
        targetprj = osr.SpatialReference(wkt = target.GetProjection())
    else:
        targetprj = target

    transform = osr.CoordinateTransformation(sourceprj, targetprj)

    to_fill = ogr.GetDriverByName("Esri Shapefile")
    # out_path = SHP_OUT + "boundary_reprojected.shp"
    ds = to_fill.CreateDataSource(out_path)
    outlayer = ds.CreateLayer('', targetprj, ogr.wkbPolygon)
    outlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))

    #apply transformation
    i = 0

    for feature in layer:
        transformed = feature.GetGeometryRef()
        transformed.Transform(transform)

        geom = ogr.CreateGeometryFromWkb(transformed.ExportToWkb())
        defn = outlayer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', i)
        feat.SetGeometry(geom)
        outlayer.CreateFeature(feat)
        i += 1
        feat = None

    ds = None
    return out_path

def get_shapefile_layer(shapefile_path):
    # CONTAINS COUNTY FILTERING
    """Returns the layer of a shapefile"""
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile_path, 0)  # 0 means read-only mode

    if data_source is None:
        print("Failed to open file")
        return

    layer = data_source.GetLayer()
    layer.SetAttributeFilter("NAMELSAD = 'Storey County' AND STATEFP = '32'")
    return layer.GetSpatialRef()

def get_shapefile_metadata(shapefile_path):
    """Prints metadata of a shapefile"""
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile_path, 0)  # 0 means read-only mode

    if data_source is None:
        print("Failed to open file")
        return

    layer = data_source.GetLayer()
    layer_definition = layer.GetLayerDefn()

    print("Layer name:", layer.GetName())
    print("Number of features:", layer.GetFeatureCount())
    print("Spatial Reference:", layer.GetSpatialRef().ExportToWkt())

    print("\nField definitions:")
    for i in range(layer_definition.GetFieldCount()):
        field_defn = layer_definition.GetFieldDefn(i)
        # print(f"  {field_defn.GetName()} ({field_defn.GetTypeName()})")

    print("\nFeatures:")
    for feature in layer:
        # print("Feature ID:", feature.GetFID())
        for i in range(layer_definition.GetFieldCount()):
            field_defn = layer_definition.GetFieldDefn(i)
            field_name = field_defn.GetName()
            field_value = feature.GetField(i)
            # print(f"  {field_name}: {field_value}")

        geometry = feature.GetGeometryRef()
        # print("  Geometry:", geometry.ExportToWkt())

    data_source = None