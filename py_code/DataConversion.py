# DataConversion.py

import os
import numpy as np
from multiprocessing import Pool, Manager
from functools import partial
from osgeo import gdal, ogr, osr
from pathlib import Path
from shapefile_funcs import *
from geometry import *

# Set the PROJ_LIB environment variable to the path of the proj data directory 
os.environ['PROJ_LIB'] = 'C:\\Users\\andre\\OneDrive_CalPoly\\Documents\\SeniorProject\\Code\\pyGeo\\lib\\site-packages\\osgeo\\data\\proj'

"""
1. Open a .tif file using gdal
2. Look at the metadata to determine what the data is
3. Convert the file from utm to lat long, play around with this after the conversion
"""


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


def get_repro_raster_list(dir, out_dir):
    raster_list = [raster for raster in dir.iterdir() if not os.path.isdir(raster)]
    args = [(raster, out_dir) for raster in raster_list]
   
    with Pool() as pool:
        reprojected_list = pool.starmap(reproject_raster, args)

    print("Rasters reprojected!")
    return reprojected_list


def reproject_raster(raster, out_dir):
    from osgeo import gdal
    gdal.UseExceptions()

    ds = gdal.Open(raster)
    repro_name = str(raster)[:-4].split("\\")[-1]
    repro_outfile = os.path.join(out_dir, f"{repro_name}_reprojected.tif")

    ds_reprojected = gdal.Warp(repro_outfile, ds, dstSRS="EPSG:4326")

    ds = None
    ds_reprojected = None
    return repro_outfile


def process_dir_parallel(dir, out_dir, shp_path, boundary = None, buffer = None):
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

    if boundary is None and buffer is None:
        boundary = get_boundary(shp_path)
        buffer = buffer_boundary(boundary, 0.01)

    county_max = county_x = county_y = 0.0
    
    # Get list of raster files (.tif files)
    # file_list = [raster for raster in dir.iterdir() if not os.path.isdir(raster)]
    reprojected_list = get_repro_raster_list(dir, out_dir) # Get list of reprojected rasters
    # filtered_rasters = filter_raster_from_list(reprojected_list, boundary, buffer) # Filter rasters that are within the boundary

    with Pool() as pool:
        # process_func = partial(process_file, boundary=boundary, buffer=buffer)
        # results = pool.map(process_func, file_list)
        # results = pool.map(process_file, filtered_rasters)  # Results is a list of tuples containing the max value and coordinates
        results = pool.map(process_file, reprojected_list)  # Results is a list of tuples containing the max value and coordinates
    # for local_max, local_x, local_y in results:
    #     # if check_point(boundary, local_x, local_y):     # Check if point is in the boundary
    #     if check_point_buffer(boundary, local_x, local_y, 0.01, buffer):     # Check if point is in the buffer
    #         if local_max > county_max:      # Update the county max
    #             county_max = local_max
    #             county_x = local_x
    #             county_y = local_y
    
    return results


def find_highest_point(reproj_raster_list, boundary, buffer, search_flag, num_points, altitude_range=None):
    """Given a list of reproj_raster_list, find the max value and cooresponding coordinates"""
    county_max = county_x = county_y = 0.0
    max_list = []
    if altitude_range is None:
        altitude_range = [0, float("inf")]

    for local_max, local_x, local_y in reproj_raster_list:
        # Check if the point is within the altitude range
        if local_max < altitude_range[0] or local_max > altitude_range[1]:
            continue  # Skip this point if it's outside the range
        
        if search_flag == "buffer":
            if check_point_buffer(boundary, local_x, local_y, 0.01, buffer):     # Check if point is in the buffer
                if local_max > county_max:      # Update the county max
                    county_max = local_max
                    county_x = local_x
                    county_y = local_y
                    entry = local_max, local_x, local_y
                    max_list.append(entry)
        elif search_flag == "boundary":
            if check_point(boundary, local_x, local_y):
                if local_max > county_max:      # Update the county max
                    county_max = local_max
                    county_x = local_x
                    county_y = local_y
                    entry = local_max, local_x, local_y
                    max_list.append(entry)
    sorted_max_list = sorted(max_list, key=lambda x: x[0], reverse=True)  # Sort the list by max value in descending order
    length = len(sorted_max_list)
    if length < num_points:
        return sorted_max_list[:length]  # Return the max value and coordinates from the sorted list
    return sorted_max_list[:num_points]  # Return the max value and coordinates from the sorted list
    

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
    
    # print("Original number of rasters: ", len(reprojected_list))

    difference = take_difference(boundary, buffer)
    filtered_rasters = []
    # difference = os.path.join(SHP_OUT, "difference.shp")

    for raster in reprojected_list:
        # Get corners of the raster
        ds = gdal.Open(raster)
        geo_transform = ds.GetGeoTransform()
        # Create a bounding box from the corners of the raster file
        bounding_box = ogr.Geometry(ogr.wkbPolygon)
        # Get the corners of the raster
        minX, maxX, minY, maxY = get_corners(ds)
        # Create a linear ring to define the bounding box
        ring = ogr.Geometry(ogr.wkbLinearRing)
        # Add points to the ring to create a rectangle
        ring.AddPoint(minX, minY)
        ring.AddPoint(minX, maxY)
        ring.AddPoint(maxX, maxY)
        ring.AddPoint(maxX, minY)
        ring.AddPoint(minX, minY)  # Close ring
        bounding_box.AddGeometry(ring)

        
        # layer = difference.GetLayer()
        # polygon = layer.GetNextFeature().GetGeometryRef()
        
        if difference.Intersects(bounding_box):
            filtered_rasters.append(raster)  # Intersection exists
    ds = None
    # print("Number of rasters within boundary: ", len(filtered_rasters))
    # print(filtered_rasters)
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


