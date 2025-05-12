# main.py
import sys
from execute import *

SINGLE_FILE = "reprojected/storey/USGS_1M_11_x28y436_NV_WestCentral_EarthMRI_2020_D20_reprojected.tif"
TEST_DIR = Path("LiDAR/Storey")
# OUT_DIR = "reprojected/storey/"
shp_path = "ShapeFiles/tl_2024_us_county.shp"
SHP_OUT = "ShapeFiles/out/"


# "NAMELSAD = 'Storey County' AND STATEFP = '32'"
# NAMELSAD = Salt Lake County AND STATEFP = 49
# LiDAR/SaltLake
# reprojected/slc/


def main():

    gdal.UseExceptions() 
    if len(sys.argv) != 2:
        print("Usage: python main.py <input_text_file>")
        sys.exit(1)
    
    input_text_file = sys.argv[1]
    

    # print(process_dir(TEST_DIR))
    
    """TESTING process_dir_parallel"""
    # test_process_dir_parallel()
    
    """TESTING IO"""
    # test_io()

    """BUFFER TESTING"""
    # test_buffer_reprojection()
    
    '''TESTING NEAREST POINT'''
    test_nearest_point(input_text_file)

    '''TESTING COUNTY INPUT'''
    # county_filter = get_county_filter()
    # print(f"County filter: {county_filter}")
    # boundary = get_boundary(shp_path, county_filter) # Get the boundary from the shapefile
    # # TETON COUNTY = 22 Wyoming = 56
    # create_shapefile(boundary, path = SHP_OUT + "TetonCounty.shp", type='polygon') # Create a shapefile from the boundary

    
if __name__ == '__main__':
    main()
