
from DataConversion import *
from shapefile_funcs import *
from geometry import *
from InOut import *

SINGLE_FILE = "reprojected/storey/USGS_1M_11_x28y436_NV_WestCentral_EarthMRI_2020_D20_reprojected.tif"
TEST_DIR = Path("LiDAR/Storey")
OUT_DIR = "reprojected/storey/"
shp_path = "ShapeFiles/tl_2024_us_county.shp"
SHP_OUT = "ShapeFiles/out/"
def test_process_dir_parallel():
    print(process_dir_parallel(TEST_DIR))

def test_buffer_reprojection():
    '''FLOW SHOULD BE:
    1. get boundary from shapefile
    2. Convert the boundary to a web mercator projection
    3. put the boundary in a shapefile
    4. buffer the boundary'''
    
    boundary = get_boundary(shp_path)
    boundary_shapefile = create_shapefile(boundary, path = SHP_OUT + "boundary.shp") # Create a shapefile from the boundary
    reprojected_list = get_repro_raster_list(TEST_DIR) # Get list of reprojected rasters
    sample_raster = reprojected_list[0] # Get a sample raster to reproject the shapefile to
    reproj_shp = reproject_shapefile(gdal.Open(sample_raster), boundary_shapefile, flag="tif") # Reproject boundary shapefile to the same as QGIS
    reproj_boundary = get_shapefile_layer(reproj_shp) # Get the reprojected boundary
    buffer = buffer_boundary(reproj_boundary, 1000) # Calculate the buffer of the polygon
    buffer_shp = create_shapefile(buffer, path = SHP_OUT + "Buffer.shp") # Create a shapefile from the buffer
    reproj_buffer = reproject_shapefile(gdal.Open(sample_raster), buffer_shp, flag="tif", out_path=SHP_OUT + "buffer_webMercator.shp") # Reproject buffer shapefile to the same as QGIS
    get_shapefile_metadata(reproj_buffer) # Get the metadata of the shapefile
    get_shapefile_metadata(reproj_shp) # Get the metadata of the shapefile

    # layer = get_shapefile_layer(shp_path) # Get the layer of the shapefile containing all counties
    # reprojected_buffer = reproject_shapefile(layer, buffer_shapefile, flag="shape", out_path = SHP_OUT + "buffer_reprojectednew.shp") # Reproject buffer shapefile to the same as QGIS
    

def test_nearest_point():
    '''This works so far, I think that distance between the point and the boundary will always be less than the radius,
    so I can filter distances that are less than the radius'''

    point, radius = print_welcome_message()
     # Get altitude range from the user
    altitude_range = get_altitude_range() # Get the altitude range from the user
    # Prompt user for a point and radius
    if point != None:
        search_radius = create_circle(point[0], point[1], radius) # Create a circle around the point
        search_in_radius(search_radius, point, altitude_range) # Search for the highest point in the radius
    else:
        search_boundary(altitude_range) # Search for the highest point in the boundary

    # CREATE SHAPEFILES
    # create_shapefile(buff_bound_line, path = SHP_OUT + "distance_line.shp", type = 'line') # Create a shapefile from the line
    # create_shapefile(HighPoint_line, path = SHP_OUT + "highPoint_distance_line.shp", type = 'line') # Create a shapefile from the line
    # create_shapefile(intersection, path = SHP_OUT + "boundary-circle.shp", type = 'polygon') # Create a shapefile from the circle
    # create_shapefile(circle, path = SHP_OUT + "circle.shp", type = 'polygon') # Create a shapefile from the circle
    # create_shapefile(center, path = SHP_OUT + "center.shp", type = 'point') # Create a shapefile from the center
    # create_shapefile(boundary_HighPoint, path = SHP_OUT + "boundary_highPoint.shp", type = 'point') # Create a shapefile from the center
    # create_gpkg(boundary, SHP_OUT + "boundary_points.gpkg") # Create a geopackage from the boundary

def search_boundary():
    # GET POLYGONS
    boundary = get_boundary(shp_path)
    buffer = buffer_boundary(boundary, .01) # Still in degrees
    difference = take_difference(boundary, buffer)

    # Get the list of reproj rasters
    reproj_raster_list = process_dir_parallel(TEST_DIR)

    # CALCULATE HIGH POINTS
    boundary_max_list = find_highest_point(reproj_raster_list, boundary, buffer, "boundary", altitude_range) # Get the max raster from the list
    buffer_max_list = find_highest_point(reproj_raster_list, boundary, buffer, "buffer", altitude_range) # Get the max raster from the list
    buffer_max = buffer_max_list[0] # Get the max value from the list
    boundary_max = boundary_max_list[0] # Get the max value from the list

    print("\nFive highest points within boundary:")
    print_high_points(boundary_max_list) # Print the high points
    print("\nFive highest points within buffer:")
    print_high_points(buffer_max_list) # Print the high points
    print(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}")
    print(f"buffer Max: {buffer_max[0]}, buffer X: {buffer_max[1]}, buffer Y: {buffer_max[2]}")

    # CREATE CIRCLE around the point in the BUFFER
    circle, buffer_HighPoint = create_circle(buffer_max[1], buffer_max[2], .01) 
    intersection = boundary.Intersection(circle) 

    # CALCULATE DISTANCE
    # Get the distance between the point and the boundary
    buf_bound_distance, buff_bound_line = distance_point_boundary(buffer_HighPoint, intersection)
    # Get the distance between the buffer high point and the boundary highpoint
    boundary_HighPoint = create_point(boundary_max[1], boundary_max[2]) # Create a point from the boundary high point
    highPoint_distance, HighPoint_line = distance_point_point(buffer_HighPoint, boundary_HighPoint) # Get the distance between the two points

    print(f"Buf boundary distance: {buf_bound_distance} meters") # Print the distance
    print(f"High Point to High Point distance: {highPoint_distance} meters") # Print the distance
    print(f"Difference in height: {abs(buffer_max[0] - boundary_max[0])} meters") # Print the difference in height
    # print(f"circle:  {circle.ExportToWkt()}") # Print the circle

    # CREATE SHAPEFILES
    # create_shapefile(buff_bound_line, path = SHP_OUT + "distance_line.shp", type = 'line') # Create a shapefile from the line
    # create_shapefile(HighPoint_line, path = SHP_OUT + "highPoint_distance_line.shp", type = 'line') # Create a shapefile from the line
    # create_shapefile(intersection, path = SHP_OUT + "boundary-circle.shp", type = 'polygon') # Create a shapefile from the circle
    # create_shapefile(circle, path = SHP_OUT + "circle.shp", type = 'polygon') # Create a shapefile from the circle
    # create_shapefile(center, path = SHP_OUT + "center.shp", type = 'point') # Create a shapefile from the center
    # create_shapefile(boundary_HighPoint, path = SHP_OUT + "boundary_highPoint.shp", type = 'point') # Create a shapefile from the center
    # create_gpkg(boundary, SHP_OUT + "boundary_points.gpkg") # Create a geopackage from the boundary


def search_in_radius(search_radius, search_point):
    """This function will search for the highest point in a radius around a point"""

    boundary = search_radius
    # Get altitude range from the user
    altitude_range = get_altitude_range() # Get the altitude range from the user
    # Get the list of reproj rasters
    reproj_raster_list = process_dir_parallel(TEST_DIR)
    # CALCULATE HIGH POINTS
    boundary_max_list = find_highest_point(reproj_raster_list, boundary, None, "boundary", altitude_range) # Get the max raster from the list
    boundary_max = boundary_max_list[0] # Get the max value from the list

    print("\nFive highest points within search radius:")
    print_high_points(boundary_max_list) # Print the high points
    print(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}")

    # create_shapefile(circle, path = SHP_OUT + "searchRadius.shp", type = 'polygon') # Create a shapefile from the circle
    # create_shapefile(center, path = SHP_OUT + "searchPoint.shp", type = 'point') # Create a shapefile from the center

    
def main():
    gdal.UseExceptions() 
    # print(process_dir(TEST_DIR))
    
    """TESTING process_dir_parallel"""
    # test_process_dir_parallel()
    

    """BUFFER TESTING"""
    # test_buffer_reprojection()

    '''TESTING NEAREST POINT'''
    test_nearest_point()
    

    
if __name__ == '__main__':
    main()