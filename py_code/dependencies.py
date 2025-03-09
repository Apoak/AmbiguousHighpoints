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
from osgeo import gdal, ogr, osr
from pathlib import Path