from DataConversion import *
from shapefile_funcs import *
from geometry import *
from InOut import *

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
    

def test_io(): 
    print_welcome_message()
    # Get input path from the user
    input_path = get_Lidar_path()
    # Get output path from the user
    output_path = get_output_path()
    # Check if the input and output paths are None
    
    if input_path is None:
        input_path = TEST_DIR
        
    if output_path is None:
        output_path = OUT_DIR
    # Get county filter from the user
    county_filter = get_county_filter()
    # Prompt user for a point and radius
    point, radius = circle_question() 
    # Get number of points to search for
    num_points = get_number_of_points() 
    # Get altitude range from the user
    altitude_range = get_altitude_range() 


def test_nearest_point(input_text_file):
    '''This works so far, I think that distance between the point and the boundary will always be less than the radius,
    so I can filter distances that are less than the radius'''
    input_path, output_path, county_name, state_filter, altitude_range, num_points, radius, point, create_shape = read_text_file(input_text_file) # Read the text file
    print_welcome_message()
    
    if input_path is None:
        input_path = TEST_DIR
        
    if output_path is None:
        output_path = OUT_DIR
    
    county_filter = "NAMELSAD = " + "'" + county_name + "'" + " AND STATEFP = " + "'" + state_filter + "'"# Filter for Salt Lake County

    print(f"County and state filter: {county_filter}") # Print the filtered county
    if point != None and radius != None: # If the point and radius are not None, search in a circle
        search_radius, center = create_circle(point[0], point[1], radius) # Create a circle around the point
        search_in_radius(search_radius, center, altitude_range, input_path, output_path, num_points) # Search for the highest point in the radius
    else:
        search_boundary(altitude_range, input_path, output_path, num_points, county_filter, create_shape, county_name) # Search for the highest point in the boundary

    
   
def search_boundary(altitude_range, input_path, output_path, num_points, county_filter, create_shape, county_name):
    """This function will search for the highest point in the boundary and buffer, writes results to output file"""
    output_file = "output.txt"

    # GET POLYGONS
    boundary = get_boundary(shp_path, county_filter)
    buffer = buffer_boundary(boundary, .01) # Still in degrees
    difference = take_difference(boundary, buffer)

    # Get the list of reproj rasters
    reproj_raster_list = process_dir_parallel(input_path, output_path, boundary, buffer)
    # print(f"Reproj raster list: {reproj_raster_list}") # Print the reproj raster list

    # CALCULATE HIGH POINTS
    boundary_max_list = find_highest_point(reproj_raster_list, boundary, buffer, "boundary", num_points, altitude_range) # Get the max raster from the list
    buffer_max_list = find_highest_point(reproj_raster_list, boundary, buffer, "buffer", num_points, altitude_range) # Get the max raster from the list
    buffer_max = buffer_max_list[0] # Get the max value from the list
    boundary_max = boundary_max_list[0] # Get the max value from the list

    # CREATE CIRCLE around the point in the BUFFER
    circle, buffer_HighPoint = create_circle(buffer_max[1], buffer_max[2], .01) 
    intersection = boundary.Intersection(circle) 

    # CALCULATE DISTANCE
    # Get the distance between the point and the boundary
    buff_bound_distance, buff_bound_line = distance_point_boundary(buffer_HighPoint, intersection)
    # Get the distance between the buffer high point and the boundary highpoint
    boundary_HighPoint = create_point(boundary_max[1], boundary_max[2]) # Create a point from the boundary high point
    highPoint_distance, HighPoint_line = distance_point_point(buffer_HighPoint, boundary_HighPoint) # Get the distance between the two points
    
    print_output(boundary_max_list, buffer_max_list, boundary_max, buffer_max, buff_bound_distance, highPoint_distance)
    county_name = county_name.replace(" ", "_")
    write_text_file(boundary_max_list, buffer_max_list, boundary_max, buffer_max, buff_bound_distance, highPoint_distance, county_name)
    
    # CREATE SHAPEFILES
    if create_shape:
        create_shapefile(buff_bound_line, path = SHP_OUT + county_name + "distance_line.shp", type = 'line') # Create a shapefile from the line
        create_shapefile(HighPoint_line, path = SHP_OUT + county_name + "highPoint_distance_line.shp", type = 'line') # Create a shapefile from the line
        create_shapefile(boundary, path = SHP_OUT + county_name + "County.shp", type = 'polygon') # Create a shapefile from the circle
        create_shapefile(circle, path = SHP_OUT + county_name + "circle.shp", type = 'polygon') # Create a shapefile from the circle
        create_shapefile(buffer_HighPoint, path = SHP_OUT + county_name +  "center.shp", type = 'point') # Create a shapefile from the center
        create_shapefile(boundary_HighPoint, path = SHP_OUT + county_name + "boundary_highPoint.shp", type = 'point') # Create a shapefile from the center
        create_gpkg(boundary, SHP_OUT + county_name + "boundary_points.gpkg") # Create a geopackage from the boundary


def search_in_radius(search_radius, search_point, altitude_range, input_path, output_path, num_points):
    """This function will search for the highest point in a radius around a point"""

    boundary = search_radius
    # Get the list of reproj rasters
    reproj_raster_list = process_dir_parallel(input_path, output_path)
    # CALCULATE HIGH POINTS
    boundary_max_list = find_highest_point(reproj_raster_list, boundary, None, "boundary", num_points, altitude_range) # Get the max raster from the list
    boundary_max = boundary_max_list[0] # Get the max value from the list

    print("\nFive highest points within search radius:")
    print_high_points(boundary_max_list) # Print the high points
    print(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}")

    # create_shapefile(search_radius, path = SHP_OUT + "searchRadius.shp", type = 'polygon') # Create a shapefile from the circle
    # create_shapefile(search_point, path = SHP_OUT + "searchPoint.shp", type = 'point') # Create a shapefile from the center


def test_boundary():
    county_filter = "NAMELSAD = 'Salt Lake County' AND STATEFP = '49'"
    boundary = get_boundary(shp_path, county_filter)
    buffer = buffer_boundary(boundary, .01) # Still in degrees
    # create_shapefile(boundary, path = SHP_OUT + "SaltLakeCounty.shp", type = 'polygon') # Create a shapefile from the circle
    # create_shapefile(buffer, path = SHP_OUT + "SaltLakeCounty_buffer.shp", type = 'polygon') # Create a shapefile from the circle


# CREATE SHAPEFILES
# create_shapefile(buff_bound_line, path = SHP_OUT + "distance_line.shp", type = 'line') # Create a shapefile from the line
# create_shapefile(HighPoint_line, path = SHP_OUT + "highPoint_distance_line.shp", type = 'line') # Create a shapefile from the line
# create_shapefile(intersection, path = SHP_OUT + "boundary_circle.shp", type = 'polygon') # Create a shapefile from the circle
# create_shapefile(circle, path = SHP_OUT + "circle.shp", type = 'polygon') # Create a shapefile from the circle
# create_shapefile(center, path = SHP_OUT + "center.shp", type = 'point') # Create a shapefile from the center
# create_shapefile(boundary_HighPoint, path = SHP_OUT + "boundary_highPoint.shp", type = 'point') # Create a shapefile from the center
# create_gpkg(boundary, SHP_OUT + "boundary_points.gpkg") # Create a geopackage from the boundary

