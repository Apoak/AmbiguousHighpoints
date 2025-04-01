from DataConversion import *
from shapefile_funcs import *

SINGLE_FILE = "LiDAR/Storey/USGS_1M_11_x27y435_NV_WestCentral_EarthMRI_2020_D20.tif"
TEST_DIR = Path("LiDAR/Storey")
OUT_DIR = "reprojected/storey/"
shp_path = "ShapeFiles/tl_2024_us_county.shp"
SHP_OUT = "ShapeFiles/out/"
# def main():
#     """
#     Works on singular raster tile
#     Changes raster projection from UTM to Lat Lon, then finds the max elevation value and the associated coordinates in lat Lon.
#     """
#     # Turn on gdal exceptions
#     gdal.UseExceptions() 

#     # Open dataset and get METADATA
#     dataset = gdal.Open(TEST_DIR)
#     get_metadata(dataset)

#     # Name an output file, call GDAL warp, change from UTM to LAT/LONG
#     output_file = "reprojected.tif"
#     reproject(output_file, dataset)
    
#     # Open the new file and get METADATA
#     out_data = gdal.Open(output_file)
#     get_metadata(out_data)

#     # Calculate the high point and associated coordinates
#     max_elevation, max_x, max_y = max_point(out_data)
#     point = ogr.Geometry(ogr.wkbPoint)
#     point.AddPoint(max_x, max_y)

#     boundary = get_boundary(shp_path)
#     point.AssignSpatialReference(None)
#     print(point.ExportToWkt())

#     # print(boundary.ExportToWkt())
#     print(boundary.GetGeometryType())
    
#     if boundary.Contains(point):
#         print("The point is inside the polygon.")
#     else:
#         print("The point is outside the polygon.")
    # Close dataset
   
def main():
    """
    HAVING TROUBLE WITH THE BUFFER_BOUNDARY FUNCTION 
    I NEED TO CHANGE THE PROJECTION OF THE BOUNDARY TO UTM
    IN ORDER TO MAKE THE BUFFER DISTANCE IN METERS
    TRY OUT THE SHAPEFILE METADATA FUNCTION.
    """
    gdal.UseExceptions() 
    # print(process_dir(TEST_DIR))
    
    """
    TESTING process_dir_parallel:
    Creates a list of reprojected files, filters them, 
    calculates max elevation and coordinates in the filtered raster tiles
    """
    print(process_dir_parallel(TEST_DIR))
    # get_shapefile_metadata(shp_path)
    # boundary = get_boundary(shp_path)
    # print(boundary)
    # boundary = feature.GetGeometryRef()
    # print(boundary.GetGeometryType())
    # get_shapefile_metadata(shp_path)
    # county = get_boundary(shp_path)
    
    # reprojected_shapefile_path = reproject_shapefile(gdal.Open(SINGLE_FILE), shp_path, flag="tif", out_path = SHP_OUT + "buffer_reprojected3.shp")
    #buffer = buffer_boundary(reprojected_shapefile_path, 10000) # Calculate the buffer of the polygon

    # THIS WORKS: CREATES A SHAPEFILE OF BUFFERED COUNTY
    
    """CONSIDER REPROJECTING THE BOUNDARIES TO UTM?"""
    """Create list of raster files, reproject them, then filter them"""
    # reprojected_list = get_repro_raster_list(TEST_DIR)
    # boundary = get_boundary(shp_path)
    # buffer = (buffer_boundary(boundary, .01))
    # filter_raster_from_list(reprojected_list, boundary, buffer)

    """BUFFER TESTING"""
    # buffer = (buffer_boundary(boundary, .01))
    # print(buffer)
    # buffer_shapefile = create_shapefile(buffer, path = SHP_OUT + "buffer3.shp") # Create a shapefile from the buffer
    # layer = get_shapefile_layer(shp_path) # Get the layer of the shapefile containing all counties
    # reprojected_buffer = reproject_shapefile(layer, buffer_shapefile, flag="shape", out_path = SHP_OUT + "buffer_reprojectednew.shp") # Reproject buffer shapefile to the same as QGIS
    

    # create_shp_file(boundary)
    # print(boundary.ExportToWkt()) 
    # #print(polygon.GetGeometryType())
    # # print(TEST_DIR)

    
if __name__ == '__main__':
    main()