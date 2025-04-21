# I DELETED THIS FROM PATH \msys64\ucrt64\binpi
import argparse
import glob
import math
import os
from pyproj import Geod
os.environ['PROJ_LIB'] = 'C:\\Users\\andre\OneDrive_CalPoly\\Documents\\SeniorProject\\Code\\pyGeo\\lib\\site-packages\\osgeo\\data\\proj'
# os.environ['GDAL_DATA'] = 'C:\\Users\\Sai kiran\\anaconda3\\envs\sai\\Library\\share'

import signal
import subprocess
import numpy as np

from multiprocessing import Pool
from functools import partial
from osgeo import gdal, ogr, osr
from pathlib import Path
from shapefile_funcs import *

"""
1. Open a .tif file using gdal
2. Look at the metadata to determine what the data is
3. Convert the file from utm to lat long, play around with this after the conversion
"""
SINGLE_FILE = "LiDAR/Storey/USGS_1M_11_x27y435_NV_WestCentral_EarthMRI_2020_D20.tif"
TEST_DIR = Path("LiDAR/Storey")
OUT_DIR = "reprojected/storey/"
shp_path = "ShapeFiles/tl_2024_us_county.shp"
SHP_OUT = "ShapeFiles/out/"


def get_corners(ds):
    """ Return list of corner coordinates from a gdal Dataset """
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    return xmin, xmax, ymin, ymax


def create_point(x, y):
    """Creates a point at the given coordinates"""
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(x, y)
    point.AssignSpatialReference(None)
    print("Point created!")
    return point.Clone()


def create_circle(x, y, radius):
    """Creates a circle around a point with a given radius"""
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(x, y)
    point.AssignSpatialReference(None)
    circle = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)

    for angle in range(0, 360):
        x = point.GetX() + radius * math.cos(math.radians(angle))
        y = point.GetY() + radius * math.sin(math.radians(angle))
        ring.AddPoint(x, y)
    ring.CloseRings()
    circle.AddGeometry(ring)
    circle.CloseRings()  # Close the ring to form a polygon
    print("Circle created!")
    return circle.Clone(), point.Clone()


def scale_boundary(dx, dy, county=None):
    """Creates a secondary boundary given a specific distance value (d) to exapnd by
            Find centroid of the county
            Scale each point in the boundary by a scaler vector"""
    if county is None:
        county = get_boundary(shp_path)

    # Get the boundary of the county
    boundary = county.GetGeometryRef(0)
    # Get the centroid of the county
    # centroid = boundary.Centroid()
    coords = np.array([boundary.GetPoint(i)[:2] for i in range(boundary.GetPointCount())])
   
    # Compute Centroid
    centroid = coords.mean(axis=0)  # Mean of x and y coordinates
    coords -= centroid # Translate to Origin

    # Apply Scaling
    scaling_matrix = np.array([[dx, 0], [0, dy]])  # 2x2 scaling matrix
    coords = coords @ scaling_matrix.T  # Apply transformation
    coords += centroid # Translate Back to Original Position

    # Convert back to OGR Polygon
    new_ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in coords:
        new_ring.AddPoint(x, y)
    
    new_polygon = ogr.Geometry(ogr.wkbPolygon)
    new_polygon.AddGeometry(new_ring)

    return new_polygon.Clone()


def get_boundary(path):
    # CONTAINS COUNTY FILTERING
    """Opens a shapefile, filters it to a specific county, and then returns the polygon
    """
    
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(path, 0)  # 0 means read-only mode

    if data_source is None:
        print("Failed to open file")
    else:
        print("Shapefile opened successfully!")
    
    layer = data_source.GetLayer()
    # print(layer.GetSpatialRef().ExportToWkt())
    # Print geometries and attributes of each feature
    layer.SetAttributeFilter("NAMELSAD = 'Storey County' AND STATEFP = '32'")
    
    # print(layer.GetFeatureCount())
    # print(layer.GetSpatialRef().ExportToWkt())
    # print(layer.GetExtent())
    for feature in layer:
        polygon = feature.GetGeometryRef()
        data_source = None
        return polygon.Clone()


def buffer_boundary(path, buffer_distance):
    # CONTAINS COUNTY FILTERING
    # PATH IS OF A SPECIFIC COUNTY
    """Converts polygon from lat/lon to cartesian coordinates in order to make distance in meters.
    Then Creates a buffer in meters around a boundary"""

    # driver = ogr.GetDriverByName("ESRI Shapefile")
    # data_source = driver.Open(path, 0)  # 0 means read-only mode
    # layer = data_source.GetLayer()
    # # layer.SetAttributeFilter("NAMELSAD = 'Storey County' AND STATEFP = '32'")
    # feature = layer.GetNextFeature()
    # boundary = feature.GetGeometryRef()

    # source_srs = osr.SpatialReference()
    # # target_srs = osr.SpatialReference()
    # # target_srs.ImportFromEPSG(32611)  # UTM Zone 11N

    # centroid = boundary.Centroid()
    # lon, lat, _ = centroid.GetPoint()

    # # Determine the UTM zone based on the longitude
    # utm_zone = int((lon + 180) / 6) + 1
    # is_northern = lat >= 0

    # target_srs = osr.SpatialReference()
    # target_srs.SetUTM(utm_zone, is_northern)

    # # Create coordinate transformation
    # transform = osr.CoordinateTransformation(source_srs, target_srs)

    # # Directly transform the polygon (no manual iteration needed!)
    # boundary.Transform(transform)
    # buffer = boundary.Buffer(buffer_distance)
    buffer = path.Buffer(buffer_distance)
    print("Buffer created!")
    return buffer.Clone()


def take_difference(boundary, buffer):
    difference = buffer.Difference(boundary)
    return difference.Clone()

def distance_point_point(point1, point2):
    """Calculates the distance between two points"""
    # Get lat and lon from points, calculate distance
    geod = Geod(ellps="WGS84")
    lon1, lat1 = point1.GetX(), point1.GetY()
    lon2, lat2 = point2.GetX(), point2.GetY()
    az12, az21, distance = geod.inv(lon1, lat1, lon2, lat2)
    # Create line
    line = ogr.Geometry(ogr.wkbLineString)
    line.AddPoint(lon1, lat1)  # Point 1: Longitude, Latitude
    line.AddPoint(lon2, lat2)  # Point 2: Longitude, Latitude
    return distance, line.Clone()  # in meters


def distance_point_boundary(center, boundary):
    """Calculates the distance between a point and a boundary"""
    # point = ogr.Geometry(ogr.wkbPoint)
    # point.AddPoint(x, y)
    # point.AssignSpatialReference(None)
    center_lon = center.GetX()
    center_lat = center.GetY()
    geod = Geod(ellps="WGS84")
    min_distance = float("inf")

    line = ogr.Geometry(ogr.wkbLineString)

    # Step 2: Add points to the line (each point is a pair of coordinates)
    
    # Iterate through the points of the boundary polygon
    for point in boundary.GetGeometryRef(0).GetPoints():
        boundary_lon = point[0]
        boundary_lat = point[1]
        az12, az21, distance = geod.inv(center_lon, center_lat, boundary_lon, boundary_lat)
        if distance < min_distance:
            min_boundary = point
            min_distance = distance
    line.AddPoint(center.GetX(), center.GetY())  # Point 1: Longitude, Latitude
    line.AddPoint(min_boundary[0], min_boundary[1])  # Point 2: Longitude, Latitude
    return min_distance, line.Clone()  # in meters


def reproject_to_web_mercator(geom):
    
    """
    Reprojects a geometry (point or polygon) from WGS84 (EPSG:4326)
    to Web Mercator (EPSG:3857).

    Parameters:
        geom (ogr.Geometry): An OGR geometry in EPSG:4326

    Returns:
        ogr.Geometry: A new reprojected geometry in EPSG:3857
    """
    if not isinstance(geom, ogr.Geometry):
        raise ValueError("Input must be an ogr.Geometry object")

    # Source: WGS84
    source_srs = osr.SpatialReference()
    source_srs.ImportFromEPSG(4326)

    # Target: Web Mercator
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(3857)

    # Coordinate transformation
    transform = osr.CoordinateTransformation(source_srs, target_srs)

    # Clone the geometry to avoid modifying original
    reprojected_geom = geom.Clone()
    reprojected_geom.Transform(transform)

    return reprojected_geom


