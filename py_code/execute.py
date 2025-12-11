from DataConversion import *
from shapefile_funcs import *
from geometry import *
from InOut import *


def nearest_point(input_text_file, SHP_OUT):
    '''This works so far, I think that distance between the point and the boundary will always be less than the radius,
    so I can filter distances that are less than the radius'''
    input_path, output_path, county_name, state_filter, altitude_range, num_points, radius, point, create_shape, shp_path = read_text_file(input_text_file) # Read the text file
    print_welcome_message()
    
    county_filter = "NAMELSAD = " + "'" + county_name + "'" + " AND STATEFP = " + "'" + state_filter + "'"# Filter for Salt Lake County

    print(f"County and state filter: {county_filter}") # Print the filtered county
    if point != None and radius != None: # If the point and radius are not None, search in a circle
        search_radius, center = create_circle(point[0], point[1], radius) # Create a circle around the point
        search_in_radius(search_radius, center, altitude_range, input_path, output_path, num_points, shp_path, SHP_OUT, county_name) # Search for the highest point in the radius
    else:
        search_boundary(altitude_range, input_path, output_path, num_points, county_filter, create_shape, county_name, shp_path, SHP_OUT) # Search for the highest point in the boundary

    
   
def search_boundary(altitude_range, input_path, output_path, num_points, county_filter, create_shape, county_name, shp_path, SHP_OUT):
    """This function will search for the highest point in the boundary and buffer, writes results to output file"""
    output_file = "output.txt"

    # GET POLYGONS
    boundary = get_boundary(shp_path, county_filter)
    buffer = buffer_boundary(boundary, .01) # Still in degrees
    difference = take_difference(boundary, buffer)

    # Get the list of reproj rasters
    reproj_raster_list = process_dir_parallel(input_path, output_path, shp_path, boundary, buffer)
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
        create_shapefile(buff_bound_line, path = SHP_OUT + county_name + "_distance_line.shp", type = 'line') # Create a shapefile from the line
        create_shapefile(HighPoint_line, path = SHP_OUT + county_name + "_highPoint_distance_line.shp", type = 'line') # Create a shapefile from the line
        create_shapefile(boundary, path = SHP_OUT + county_name + "_County.shp", type = 'polygon') # Create a shapefile from the circle
        create_shapefile(circle, path = SHP_OUT + county_name + "_circle.shp", type = 'polygon') # Create a shapefile from the circle
        create_shapefile(buffer_HighPoint, path = SHP_OUT + county_name +  "_center.shp", type = 'point') # Create a shapefile from the center
        create_shapefile(boundary_HighPoint, path = SHP_OUT + county_name + "_boundary_highPoint.shp", type = 'point') # Create a shapefile from the center
        create_gpkg(boundary, SHP_OUT + county_name + "_boundary_points.gpkg") # Create a geopackage from the boundary
        create_shapefile(buffer, path = SHP_OUT + county_name + "_Buffer.shp") # Create a shapefile from the buffer

def search_in_radius(search_radius, search_point, altitude_range, input_path, output_path, num_points, shp_path, SHP_OUT, county_name):
    """This function will search for the highest point in a radius around a point"""

    boundary = search_radius
    # Get the list of reproj rasters
    reproj_raster_list = process_dir_parallel(input_path, output_path, shp_path)
    # CALCULATE HIGH POINTS
    boundary_max_list = find_highest_point(reproj_raster_list, boundary, None, "boundary", num_points, altitude_range) # Get the max raster from the list
    boundary_max = boundary_max_list[0] # Get the max value from the list

    print(f"\n{num_points} highest points within search radius:")
    print_high_points(boundary_max_list) # Print the high points
    print(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}")
    county_name = county_name.replace(" ", "_")
    write_text_file_no_buffer(boundary_max_list, boundary_max, county_name) # Write the results to a text file

    create_shapefile(search_radius, path = SHP_OUT + "searchRadius.shp", type = 'polygon') # Create a shapefile from the circle
    create_shapefile(search_point, path = SHP_OUT + "searchPoint.shp", type = 'point') # Create a shapefile from the center


# # CREATE SHAPEFILES
# create_shapefile(buff_bound_line, path = SHP_OUT + "distance_line.shp", type = 'line') # Create a shapefile from the line
# create_shapefile(HighPoint_line, path = SHP_OUT + "highPoint_distance_line.shp", type = 'line') # Create a shapefile from the line
# create_shapefile(intersection, path = SHP_OUT + "boundary_circle.shp", type = 'polygon') # Create a shapefile from the circle
# create_shapefile(circle, path = SHP_OUT + "circle.shp", type = 'polygon') # Create a shapefile from the circle
# create_shapefile(center, path = SHP_OUT + "center.shp", type = 'point') # Create a shapefile from the center
# create_shapefile(boundary_HighPoint, path = SHP_OUT + "boundary_highPoint.shp", type = 'point') # Create a shapefile from the center
# create_gpkg(boundary, SHP_OUT + "boundary_points.gpkg") # Create a geopackage from the boundary

