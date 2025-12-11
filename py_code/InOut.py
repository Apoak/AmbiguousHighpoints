from pathlib import Path

def print_welcome_message():
    """Prints a welcome message"""
    print("\nWelcome to the High Point Holmes, LiDAR Analysis Program!\n")
    print("This program will help you find the highest points in your LiDAR data.")
    # print("Please follow the prompts to enter your data.")
    # print("\n") # Print a new line for spacing


def circle_question():
    answer = input("\nWould you like to specify a point of interest and a search radius? (y/n): ") # Ask the user if they want to see the instructions
    
    if answer == "":
        return None, None # Return None if the user does not want to specify a point or radius
    
    elif answer.lower() == "y":
        point = get_point() # Get the point from the user
        radius = get_radius() # Get the radius from the user
        return point, radius
    
    elif answer.lower() == "n":
        return None, None # Return None if the user does not want to specify a point or radius
    
    else:
        print("Please enter 'y' or 'n'")
        return circle_question()


def get_Lidar_path():
    """Prompts the user for the path to the LiDAR data"""
    print("\nEnter the path to the LiDAR data.")
    print(" * To use the default path, press enter.")
    print(" * To enter a path, type it in the form of /path/to/data")
    lidar_path = input("\nEnter the path to the LiDAR data: ") # Get the path from the user
    
    if lidar_path == "":
        return None # Return None if the user does not want to specify a path
    
    return Path(lidar_path) # Return the path


def get_output_path():
    """Prompts the user for the path to the output directory"""
    print("\nEnter the path to the output directory.")
    print(" * To use the default path, press enter.")
    print(" * To enter a path, type it in the form of /path/to/data")
    output_path = input("\nEnter the path to the output directory: ") # Get the path from the user
    
    if output_path == "":
        return None # Return None if the user does not want to specify a path
    
    return output_path # Return the path


def get_altitude_range():
    print("\nEnter the altitude range in meters, the default range is 0 to infinity.")
    print(" * To use the default range, press enter.")
    print(" * To enter a range, type it in the form of 0,1000 ...")
    print(" * Type 'inf' for infinity")
    print(" * Example: 0,inf or 0,1000")
    altitude_range = input("\nEnter the altitude range in meters: ") # Get the altitude range from the user
    
    if altitude_range == "":
        print("Using default range of 0 to infinity")
        return [0, float("inf")] # Return the default range
    
    altitude_range = altitude_range.split(",") # Split the altitude range into a list
    altitude_range = [float(i) for i in altitude_range] # Convert the list to floats
   
    print(f"Using altitude range of {altitude_range[0]} to {altitude_range[1]}\n\n")
    return altitude_range # Return the altitude range


def get_point():
    print("\nEnter the point of interest in the form of latitude and longitude.")
    print(" * To enter a point, type it in the form of 39.261840852281, -119.7051741128619")
    print(" * Note: the point must not be inside of provided LiDAR data.")
    point = input("\nEnter the point of interest: ") # Get the point from the user
    
    if point == "":
        print("Please enter a point")
        return get_point()
    
    point = point.split(",") # Split the point into a list

    try:
        point = [float(i) for i in point] # Convert the list to floats
    except ValueError:
        print("\nPoint must be a valid decimal number")
        return get_point()
    
    if len(point) != 2:
        print("Point must be in the form of x,y")
        return get_point() # Get the point again if it is not in the correct format
    
    print(f"Using point of {point[0]}, {point[1]}\n")
    return point # Return the point


def get_radius():
    print("\nEnter the radius in meters, the default radius is 1000 meters.")
    print(" * To use the default radius, press enter.")
    print(" * To enter a radius, type it in the form of 1000")
    print(" * Note: the radius must not be inside of provided LiDAR data.")
    radius = input("\nEnter the radius in meters: ") # Get the radius from the user
    
    if radius == "":
        return get_radius() # Return the default radius
    
    radius = float(radius) # Convert the radius to a float
    if radius <= 0:
        print("Radius must be greater than 0")
        return get_radius() # Get the radius again if it is less than or equal to 0
    
    print(f"Using radius of {radius} meters\n")
    return radius # Return the radius


def get_county_filter():
    """Prompts the user for a county filter"""
    # NAMELSAD = 'Storey County' AND STATEFP = '32'
    county = "NAMELSAD = " # Initialize the county filter
    state = " AND STATEFP = " # Initialize the state filter
    nameLsad = input("\nEnter the county NAMELSAD: ") # Get the county NAMELSAD from the user
    stateEFP = input("Enter the state FIPS code: ") # Get the state FIPS code from the user
    
    if nameLsad == "" or stateEFP == "":
        print("No county filter provided")
        return None
    
    county_filter = county + "'" + nameLsad + "'" + state + "'" + stateEFP + "'" # Create the county filter
    return county_filter # Return the county filter


def get_number_of_points():
    """Prompts the user for the number of points to search for"""
    print("\nEnter the number of points to search for, the default is 5.")
    print(" * To use the default number of points, press enter.")
    print(" * To enter a number, type it in the form of 5")
    num_points = input("\nEnter the number of points to search for: ") # Get the number of points from the user
    
    if num_points == "":
        return 5 # Return the default number of points
    
    try:
        num_points = int(num_points) # Convert the number of points to an int
    except ValueError:
        print("\nNumber of points must be a valid integer")
        return get_number_of_points()
    
    if num_points <= 0:
        print("Number of points must be greater than 0")
        return get_number_of_points() # Get the number of points again if it is less than or equal to 0
    
    print(f"Using {num_points} points\n")
    return num_points # Return the number of points


def print_high_points(HighPoint_list):
    """Prints the high point of a geometry"""
    for count, point in enumerate(HighPoint_list):
        print(f"{count + 1}. Elevation: {point[0]} Lat: {point[1]} Lon: {point[2]}")
    print("\n") # Print a new line for spacing


def read_text_file(file_path):
    """Reads a text file containing input directory, 
    output directory, county filter, state filter, altitude range, 
    and number of points, radius, and point of interest"""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Parse the lines into variables
    try:
        lidar_path = Path(lines[0].strip())
    except ValueError:
        print("Invalid path provided.")
        print("Please provide a valid path to the LiDAR data.")
        lidar_path = get_Lidar_path()

    output_path = lines[1].strip()
    county_filter = lines[2].strip()
    state_filter = lines[3].strip()
    try:
        altitude_range = list(map(float, lines[4].strip().split(',')))
    except ValueError:
        print("Using default range of 0 to infinity.")
        altitude_range = [0, float("inf")]

    try:
        num_points = int(lines[5].strip())
    except ValueError:
        print("Using default of 5.")
        num_points = 5
    try:
        radius = float(lines[6].strip())
    except ValueError:
        print("No radius provided.")
        radius = None
    try:
        point_of_interest = list(map(float, lines[7].strip().split(',')))
    except ValueError:
        print("No point of interest provided.")
        point_of_interest = None
    
    bool = lines[8].strip().lower() 
    if bool == "true":
        create_shape = True
    elif bool == "false":
        create_shape = False
    else:
        print("Invalid input for create shape. Defaulting to False.")
        create_shape = False
    shp_path = Path(lines[9].strip()) if len(lines) > 9 else None
    if not lidar_path.exists():
        print(f"LiDAR path {lidar_path} does not exist. Please provide a valid path.")
        exit(1)  # Exit if the LiDAR path does not exist
    if not Path(output_path).exists():
        print(f"Output path {output_path} does not exist. Please provide a valid path.")
        exit(1)  # Exit if the output path does not exist
    if not shp_path or not shp_path.exists():
        print(f"Shapefile path {shp_path} does not exist. Please provide a valid path.")
        exit(1)  # Exit if the shapefile path does not exist

    return lidar_path, output_path, county_filter, state_filter, altitude_range, num_points, radius, point_of_interest, create_shape, shp_path


def write_text_file(boundary_max_list, buffer_max_list, boundary_max, buffer_max, buff_bound_distance, highPoint_distance, county_name):
    """Writes the high points to a text file"""
    with open(county_name + "_high_points.txt", 'w') as file:
        file.write(county_name + "\n")
        file.write("High Points in boundary:\n")
        # count = 1
        for count, point in enumerate(boundary_max_list):
            file.write(f"{count + 1}. Elevation: {point[0]} Lat: {point[1]} Lon: {point[2]}\n")
            # file.write(f"{count + 1}. Elevation: {point}\n")
            # count += 1
        file.write("\n")
        file.write("High Points in buffer:\n")
        count = 1
        for count, point in enumerate(buffer_max_list):
            file.write(f"{count + 1}. Elevation: {point[0]} Lat: {point[1]} Lon: {point[2]}\n")
            # file.write(f"{count + 1}. Elevation: {point}\n")
            count += 1
        file.write("\n")
        file.write(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}\n")
        file.write(f"Buffer Max: {buffer_max[0]}, buffer X: {buffer_max[1]}, buffer Y: {buffer_max[2]}\n")
        
        file.write(f"Distance from buffer high point to boundary: {buff_bound_distance}\n")
        file.write(f"Distance from buffer high point to boundary high point: {highPoint_distance}\n")

def write_text_file_no_buffer(boundary_max_list, boundary_max, county_name):
    """Writes the high points to a text file without buffer"""
    with open(county_name + "_high_points.txt", 'w') as file:
        file.write(county_name + "\n")
        file.write("High Points in boundary:\n")
        for count, point in enumerate(boundary_max_list):
            file.write(f"{count + 1}. Elevation: {point[0]} Lat: {point[1]} Lon: {point[2]}\n")
        file.write("\n")
        file.write(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}\n")


def print_no_buffer_output(boundary_max_list, boundary_max, county_name):
    """Prints the high points without buffer"""
    print("\nFive highest points within boundary:")
    print_high_points(boundary_max_list) # Print the high points
    print(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}")
    county_name = county_name.replace(" ", "_")
    write_text_file_no_buffer(boundary_max_list, boundary_max, county_name)

    
def print_output(boundary_max_list, buffer_max_list, boundary_max, buffer_max, buff_bound_distance, highPoint_distance):
    print("\nFive highest points within boundary:")
    print_high_points(boundary_max_list) # Print the high points
    print("\nFive highest points within buffer:")
    print_high_points(buffer_max_list) # Print the high points
    print(f"Boundary Max: {boundary_max[0]}, Boundary X: {boundary_max[1]}, Boundary Y: {boundary_max[2]}")
    print(f"buffer Max: {buffer_max[0]}, buffer X: {buffer_max[1]}, buffer Y: {buffer_max[2]}")
    print(f"Buf boundary distance: {buff_bound_distance} meters") # Print the distance
    print(f"High Point to High Point distance: {highPoint_distance} meters") # Print the distance
    print(f"Difference in height: {abs(buffer_max[0] - boundary_max[0])} meters") # Print the difference in height
