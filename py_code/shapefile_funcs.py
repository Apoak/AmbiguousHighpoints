import os
os.environ['PROJ_LIB'] = 'C:\\Users\\andre\\OneDrive_CalPoly\\Documents\\SeniorProject\\Code\\pyGeo\\lib\\site-packages\\osgeo\\data\\proj'

import numpy as np

from multiprocessing import Pool
from osgeo import gdal, ogr, osr
from pathlib import Path

def create_shapefile(object, path, type = None):
    """Create a shapefile from a boundary"""

    # path = SHP_OUT + "boundary_buff.shp"
    # Set up shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")
    # Create data source (shapefile)
    data_source = driver.CreateDataSource(path)

    if data_source is None:
        print("Failed to create file")
    else:
        print("Shapefile created successfully!")

    # Create spatial reference, WGS84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    
    if type == "polygon" or type == None:
        layer = data_source.CreateLayer("boundary", srs, ogr.wkbPolygon)
    if type == "point":
        layer = data_source.CreateLayer("point", srs, ogr.wkbPoint)
    if type == "line":
        layer = data_source.CreateLayer("line", srs, ogr.wkbLineString)
    if type == "points":
        layer = data_source.CreateLayer("points", srs, ogr.wkbMultiPoint)
    # Add an ID field
    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    layer.CreateField(id_field)

    # Create a new feature
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetField("id", 1)
    feature.SetGeometry(object)

    layer.CreateFeature(feature)

    # Save and close everything
    data_source = None
    return path


def get_shapefile_layer(shapefile_path):
    # CONTAINS COUNTY FILTERING
    """Returns the layer of a shapefile"""
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile_path, 0)  # 0 means read-only mode

    if data_source is None:
        print("Failed to open file")
        return

    layer = data_source.GetLayer()
    # layer.SetAttributeFilter("NAMELSAD = 'Storey County' AND STATEFP = '32'")
    for feature in layer:
        polygon = feature.GetGeometryRef()
        data_source = None
        return polygon.Clone()
    # return layer.GetSpatialRef()


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


def shp_to_geojson(shapefile_path, out_path):
    """Convert a shapefile to GeoJSON"""
    driver = ogr.GetDriverByName("GeoJSON")
    data_source = driver.CreateDataSource(out_path)

    if data_source is None:
        print("Failed to create GeoJSON file")
        return

    shp_driver = ogr.GetDriverByName("ESRI Shapefile")
    shp_data_source = shp_driver.Open(shapefile_path, 0)  # 0 means read-only mode

    if shp_data_source is None:
        print("Failed to open shapefile")
        return

    layer = shp_data_source.GetLayer()
    out_layer = data_source.CopyLayer(layer, "layer_name")

    if out_layer is None:
        print("Failed to copy layer")
        return

    data_source = None
    shp_data_source = None
    print("Shapefile converted to GeoJSON successfully!")


def geojson_to_shp(geojson_path, out_path):
    """Convert a GeoJSON file to a shapefile"""
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(out_path)

    if data_source is None:
        print("Failed to create shapefile")
        return

    geojson_driver = ogr.GetDriverByName("GeoJSON")
    geojson_data_source = geojson_driver.Open(geojson_path, 0)  # 0 means read-only mode

    if geojson_data_source is None:
        print("Failed to open GeoJSON file")
        return

    layer = geojson_data_source.GetLayer()
    out_layer = data_source.CopyLayer(layer, "layer_name")

    if out_layer is None:
        print("Failed to copy layer")
        return

    data_source = None
    geojson_data_source = None
    print("GeoJSON converted to shapefile successfully!")


def create_gpkg(object, path):
    # Create the geometry: a linear ring + polygon
    
    object = object.GetGeometryRef(0) # Get the first geometry from the shapefile
    # Set up GPKG driver
    driver = ogr.GetDriverByName("GPKG")
    gpkg = driver.CreateDataSource(path)

    # Spatial reference (WGS84)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    # Create layer with point geometry
    layer = gpkg.CreateLayer('Boundary points', srs, ogr.wkbPoint)

    # Add an ID field
    field_id = ogr.FieldDefn("ID", ogr.OFTInteger)
    layer.CreateField(field_id)

    # Add points to the layer
    points = np.array([object.GetPoint(i)[:2] for i in range(object.GetPointCount())])
    print(points)
    for i, (x, y) in enumerate(points):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x, y)

        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(point)
        feature.SetField("ID", i + 1)
        layer.CreateFeature(feature)
        feature = None  # clean up

    gpkg = None  # save and close
    print(f"Saved {len(points)} points to {path}")
