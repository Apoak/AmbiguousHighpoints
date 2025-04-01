# I DELETED THIS FROM PATH \msys64\ucrt64\binpi
import argparse
import glob
import math
import os
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
# "C:\Users\andre\OneDrive_CalPoly\Documents\SeniorProject\Code\LiDAR\USGS_1M_11_x54y441_NV_EastCentral_2021_D21.tif"
# This is kinda helpful to see how to open stuff with gdal and
# reproject the data

def run_command(command_string):
    print("> " + command_string)
    retval = subprocess.call(command_string, shell=True)
    return


    
# def get_extent(ds):
#     """ Return list of corner coordinates from a gdal Dataset """
#     xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
#     width, height = ds.RasterXSize, ds.RasterYSize
#     xmax = xmin + width * xpixel
#     ymin = ymax + height * ypixel

#     return xmin, xmax, ymin, ymax

# def create_vrts(tile_dir, input_files, skip_boundary):
#     """
#     Creates virtual rasters (VRTs) for the input files on disk.
#     Returns:
#     1) The filename for a VRT of all of the input files in EPSG:4326
#     2) The bounding polygon for all of the input files in EPSG:4326
#     """
#     print("Creating VRTs and computing boundary")

#     # Some input datasets have more than one projection.
#     #
#     # All files in a VRT must have the same projection.  Thus, we must first
#     # create an intermediate VRT that warps the inputs.
#     #
#     # However, computing the boundary can be much faster on the raw
#     # data, because we can use overview (pyramid) images at lower
#     # resolution, and these would not be present in the warped
#     # version.  So, for each projection, we create a VRT of all the
#     # inputs with that projection, compute the boundary of each one,
#     # and then we take the union of all of these boundaries to get the
#     # boundary of all of the inputs.
    
#     # Map of projection name to list of files with that projection
#     projection_map = {}  
#     for input_file in input_files:
#         ds = gdal.Open(input_file)
#         prj = ds.GetProjection()
#         srs = osr.SpatialReference(wkt=prj)
#         # gdalbuildvrt needs the projections to match pretty much exactly within a VRT file.
#         # It's not enough to use the projection name; sometimes two projections have the
#         # same name but slightly different datums, and gdalbuildvrt will fail.
#         # So instead we use the whole WKT string.
#         projection_name = srs.ExportToWkt()
#         projection_map.setdefault(projection_name, []).append(input_file)

#     return projection_map

# Maybe I can put all of the file names in a list and then do this with threads

def process_dir(dir):
    """Given a directory containing LiDAR files. For all files:
    reproject all files to lat/lon, find the max value and cooresponding coordinates.
    then compare results to find and return true max.
    """
    suffix = "_reprojected"
    file_type = ".tif"
    boundary = get_boundary(shp_path)
    county_max = county_x = county_y = 0.0
    
    for raster in dir.iterdir():
        if not os.path.isdir(raster):
            # rename outfile to xyz_reprojected.tif
            # ds = gdal.Open(raster)
            # name = str(raster)[:-4].split("\\")[-1]
            # outfile = os.path.join(OUT_DIR, f"{name}_reprojected.tif")
            
            outfile = reproject_raster(raster)
            repro_raster = gdal.Open(outfile)   # Open the file of the new preojection

            local_max, local_x, local_y = max_point(repro_raster)   # Get max of the current raster

            if check_point(boundary, local_x, local_y):     # Check if point is in the boundary
                if local_max > county_max:      # Update the county max
                    county_max = local_max
                    county_x = local_x
                    county_y = local_y
            ds = None
    
    return county_max, county_x, county_y

def get_repro_raster_list(dir):
    """Given a directory containing LiDAR files. 
    For all files:
      reproject all files to lat/lon"""
    
    raster_list = [raster for raster in dir.iterdir() if not os.path.isdir(raster)]
    with Pool() as pool:
        reprojected_list = pool.map(reproject_raster, raster_list)
    # print(reprojected_list)
    return reprojected_list

def process_dir_parallel(dir):
    """
    Given a directory containing LiDAR files. 
    For all files:
        Reproject all files to lat/lon
        Filter rasters:
            * keep rasters that have at least one corner in the boundary
        Find the max value and cooresponding coordinates
        Then compare results to find and return true max
    """
    suffix = "_reprojected"
    file_type = ".tif"

    boundary = get_boundary(shp_path)
    buffer = buffer_boundary(boundary, 0.01)

    county_max = county_x = county_y = 0.0
    
    # Get list of raster files (.tif files)
    # file_list = [raster for raster in dir.iterdir() if not os.path.isdir(raster)]
    reprojected_list = get_repro_raster_list(dir) # Get list of reprojected rasters
    filtered_rasters = filter_raster_from_list(reprojected_list, boundary, buffer) # Filter rasters that are within the boundary

    with Pool() as pool:
        # process_func = partial(process_file, boundary=boundary, buffer=buffer)
        # results = pool.map(process_func, file_list)
        results = pool.map(process_file, filtered_rasters)  # Results is a list of tuples containing the max value and coordinates

    for local_max, local_x, local_y in results:
        # if check_point(boundary, local_x, local_y):     # Check if point is in the boundary
        if check_point_buffer(boundary, local_x, local_y, 0.01):     # Check if point is in the buffer
            if local_max > county_max:      # Update the county max
                county_max = local_max
                county_x = local_x
                county_y = local_y
    
    return county_max, county_x, county_y
    
def process_file(raster):
    """Helper func for proccess parallel: Given a raster file, reproject it to lat/lon, find the max value and cooresponding coordinates."""    
    from osgeo import gdal  # Ensure it's imported inside the function for safety
    gdal.UseExceptions()
    ds = gdal.Open(raster)
    # outfile = reproject_raster(raster)
    # repro_raster = gdal.Open(outfile)   # Open the file of the new preojection
    local_max, local_x, local_y = max_point(ds)   # Get max of the current raster
    ds = None
    return local_max, local_x, local_y

def filter_raster(raster, boundary, buffer):
    """Given a raster in the same projection as the boundary(buffer), check if any of the corners of the raster are within the boundary"""
    # xmin, xmax, ymin, ymax = get_corners(raster)
    xmin, xmax, ymin, ymax = get_corners(raster)
    # Check if corners are in buffer but not in boundary
    if check_point_buffer(boundary, xmin, ymin, 0.01, buffer) or check_point_buffer(boundary, xmin, ymax, 0.01, buffer) or check_point_buffer(boundary, xmax, ymin, 0.01, buffer) or check_point_buffer(boundary, xmax, ymax, 0.01, buffer):
        return True
    return False

def filter_raster_from_list(reprojected_list, boundary, buffer):
    """Given a raster in the same projection as the boundary(buffer), check if any of the corners of the raster are within the boundary"""
    print("Original number of rasters: ", len(reprojected_list))
    # xmin, xmax, ymin, ymax = get_corners(raster)
    
    filtered_rasters = []
    for raster in reprojected_list:
        # Get corners of the raster
        ds = gdal.Open(raster)
        xmin, xmax, ymin, ymax = get_corners(ds)
        # Check if corners are in buffer but not in boundary
        if check_point_buffer(boundary, xmin, ymin, 0.01, buffer) or check_point_buffer(boundary, xmin, ymax, 0.01, buffer) or check_point_buffer(boundary, xmax, ymin, 0.01, buffer) or check_point_buffer(boundary, xmax, ymax, 0.01, buffer):
        # if check_point(boundary, xmin, ymin) or check_point(boundary, xmin, ymax) or check_point(boundary, xmax, ymin) or check_point(boundary, xmax, ymax):
            filtered_rasters.append(raster)
    print("Number of rasters within boundary: ", len(filtered_rasters))
    print(filtered_rasters)
    return filtered_rasters

    # print("Number of rasters within boundary: ", len(filtered_rasters))
    # return filtered_rasters

def check_point(boundary, x, y):
    """Given a boundary and x,y coordinates create a point and determine if it's in the boundary"""
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(x, y)
    point.AssignSpatialReference(None)
    
    if boundary.Contains(point):
        return True
        # print("The point is inside the polygon.")
    else:
        return False
        # print("The point is outside the polygon.")
    
def check_point_buffer(boundary, x, y, buffer_distance, buffer=None):
    """Given a boundary and x,y coordinates create a point and determine 
    if it's in the buffer BUT NOT IN THE boundary"""
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(x, y)
    point.AssignSpatialReference(None)

    if buffer is None:
        buffer = boundary.Buffer(buffer_distance)
    
    if buffer.Contains(point) and not boundary.Contains(point):
        return True
        # print("The point is inside the polygon.")
    else:
        return False
        # print("The point is outside the polygon.")


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
    return buffer.Clone()


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


# def reproject_raster(output_file, ds):
def reproject_raster(raster):
    """Reproject a raster dataset to lat/lon (EPSG:4326)"""
    from osgeo import gdal  
    gdal.UseExceptions()

    # Rename the output file to xyz_reprojected.tif
    ds = gdal.Open(raster)
    name = str(raster)[:-4].split("\\")[-1]        
    outfile = os.path.join(OUT_DIR, f"{name}_reprojected.tif")

    # Below this comment works 

    # Reproject the raster to lat/lon
    ds_reprojected = gdal.Warp(outfile, ds, dstSRS="EPSG:4326")
    ds = None
    ds_reprojected = None
    return outfile


def get_boundary(path):
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


def get_corners(ds):
    """ Return list of corner coordinates from a gdal Dataset """
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    return xmin, xmax, ymin, ymax


def get_metadata(dataset):
    """Prints metadata of a dataset"""
    metadata = dataset.GetMetadata()
    print("metadata: ", metadata)
    print("\nProjection ", dataset.GetProjection())  # Check coordinate reference system
    print("\nRaster X, Raster Y ", dataset.RasterXSize, dataset.RasterYSize)  # Check image dimensions
    print("Raster Count", dataset.RasterCount) 
    print("\nGeoTransform: ", dataset.GetGeoTransform())
    # Try conerting the data type here:


def max_point(dataset):# Put the elevations into a numpy array
    """Find the max elevation value and the associated coordinates in lat Lon"""
    data_array = dataset.GetRasterBand(1).ReadAsArray()
    max_index = np.unravel_index(np.argmax(data_array), data_array.shape)  # (row, col)

    # Max height
    max_value = data_array[max_index]

    # Get the raster geotransform
    geotransform = dataset.GetGeoTransform()

    # Convert index (row, col) to real-world coordinates
    x_origin, pixel_width, _, y_origin, _, pixel_height = geotransform
    max_x = x_origin + max_index[1] * pixel_width
    max_y = y_origin + max_index[0] * pixel_height  # Note: y decreases as row index increases

    # Print results
    # print(f"Max Value: {max_value}") # Previously 2329.8767
    # print(f"Max Value Coordinates: ({max_x}, {max_y})")

    return max_value, max_x, max_y
