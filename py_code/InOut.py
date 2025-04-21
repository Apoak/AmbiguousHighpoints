def print_welcome_message():
    """Prints a welcome message"""
    print("\nWelcome to the High Point Holmes, LiDAR Analysis Program!\n")
    print("This program will help you find the highest points in your LiDAR data.")
    print("Please follow the prompts to enter your data.")
    print("\n") # Print a new line for spacing

    answer = input("Would you like to specify a point of interest and a search radius? (y/n): ") # Ask the user if they want to see the instructions
   
    if answer.lower() == "y":
        point = get_point() # Get the point from the user
        radius = get_radius() # Get the radius from the user
        return point, radius
    
    elif answer.lower() == "n":
        return None, None # Return None if the user does not want to specify a point or radius


def prompt_input():
    """Prompts the user for input and returns the file name"""
    print("Please enter the path to the LiDAR data file:")
    print("The file must be in .las or .laz format.")
    print("Example: /path/to/file.las")
    file_name = input("\nEnter the path to the LiDAR data file: ") # Get the file name from the user
    
    if file_name == "":
        print("No file name provided, exiting...")
        exit() # Exit if no file name is provided
    
    return file_name # Return the file name


def get_altitude_range():
    print("\nEnter the altitude range in meters, the default range is 0 to infinity.")
    print("To use the default range, press enter.")
    print("To enter a range, type it in the form of 0,1000 ...")
    print("Type 'inf' for infinity")
    print("Example: 0,inf or 0,1000")
    altitude_range = input("\nEnter the altitude range in meters: ") # Get the altitude range from the user
    
    if altitude_range == "":
        print("Using default range of 0 to infinity")
        return [0, float("inf")] # Return the default range
    
    altitude_range = altitude_range.split(",").strip() # Split the altitude range into a list
    altitude_range = [float(i) for i in altitude_range] # Convert the list to floats
   
    print(f"Using altitude range of {altitude_range[0]} to {altitude_range[1]}\n")
    return altitude_range # Return the altitude range


def get_point():
    print("\nEnter the point of interest in the form of latitude and longitude.")
    print("To enter a point, type it in the form of 39.261840852281, -119.7051741128619")
    print("Note: the point must not be inside of provided LiDAR data.")
    point = input("\nEnter the point of interest: ") # Get the point from the user
    
    if point == "":
        print("Please enter a point")
        return get_point()
    
    point = point.split(",").strip() # Split the point into a list
    point = [float(i) for i in point] # Convert the list to floats
    
    if len(point) != 2:
        print("Point must be in the form of x,y")
        return get_point() # Get the point again if it is not in the correct format
    
    print(f"Using point of {point[0]}, {point[1]}\n")
    return point # Return the point


def get_radius():
    print("\nEnter the radius in meters, the default radius is 1000 meters.")
    print("To use the default radius, press enter.")
    print("To enter a radius, type it in the form of 1000")
    print("Note: the radius must not be inside of provided LiDAR data.")
    radius = input("\nEnter the radius in meters: ") # Get the radius from the user
    
    if radius == "":
        print("Using default radius of 1000 meters")
        return 1000 # Return the default radius
    
    radius = float(radius) # Convert the radius to a float
    if radius <= 0:
        print("Radius must be greater than 0")
        return get_radius() # Get the radius again if it is less than or equal to 0
    
    print(f"Using radius of {radius} meters\n")
    return radius # Return the radius


def print_high_points(HighPoint_list):
    """Prints the high point of a geometry"""
    for count, tuple in enumerate(HighPoint_list):
        print(f"{count + 1}. Elevation: {tuple[0]} Lat: {tuple[1]} Lon: {tuple[2]}")

