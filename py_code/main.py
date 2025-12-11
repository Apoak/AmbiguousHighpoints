# main.py
import sys
from execute import *




def main():
    SHP_OUT = "ShapeFiles/out/"
    # Check if SHP_OUT directory exists, if not create it
    if not Path(SHP_OUT).exists():
        Path(SHP_OUT).mkdir(parents=True, exist_ok=True)
    
    gdal.UseExceptions() 
    if len(sys.argv) != 2:
        print("Usage: python main.py <input_text_file>")
        sys.exit(1)
    
    input_text_file = sys.argv[1]
    
    '''NEAREST POINT'''
    nearest_point(input_text_file, SHP_OUT)
    
if __name__ == '__main__':
    main()
