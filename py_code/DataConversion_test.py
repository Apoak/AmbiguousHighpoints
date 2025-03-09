import unittest
import DataConversion
from osgeo import gdal, ogr, osr
import os

class TestDataConversion(unittest.TestCase):

    def test_make_boundary(self):
        dx = dy = 2
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(0, 0)
        ring.AddPoint(4, 0)
        ring.AddPoint(4, 4)
        ring.AddPoint(0, 4)
        ring.AddPoint(0, 0)  # Close the ring
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)

        test = ogr.Geometry(ogr.wkbLinearRing)
        test.AddPoint(0, 0)
        test.AddPoint(4, 0)
        test.AddPoint(4, 4)
        test.AddPoint(0, 4)
        test.AddPoint(0, 0)  # Close the ring
        test_polygon = ogr.Geometry(ogr.wkbPolygon)
        test_polygon.AddGeometry(test)
        self.assertEqual(DataConversion.make_boundary(dx, dy, polygon), test_polygon)

if __name__ == '__main__':
    unittest.main()

import unittest
from osgeo import ogr
from DataConversion import get_boundary, buffer_boundary, scale_boundary, create_shapefile, check_point

class TestDataConversion(unittest.TestCase):

    def setUp(self):
        self.shp_path = "ShapeFiles/tl_2024_us_county.shp"
        self.boundary = get_boundary(self.shp_path)

    def test_get_boundary(self):
        self.assertIsNotNone(self.boundary, "Boundary should not be None")

    def test_buffer_boundary(self):
        buffer_distance = 1000
        buffered_boundary = buffer_boundary(self.boundary, buffer_distance)
        self.assertIsNotNone(buffered_boundary, "Buffered boundary should not be None")

    def test_scale_boundary(self):
        dx, dy = 1.1, 1.1
        scaled_boundary = scale_boundary(dx, dy, self.boundary)
        self.assertIsNotNone(scaled_boundary, "Scaled boundary should not be None")

    def test_create_shp_file(self):
        create_shapefile(self.boundary)
        # Check if the shapefile was created (you can add more checks if needed)
        self.assertTrue(os.path.exists("ShapeFiles/out/boundary_buff.shp"), "Shapefile should be created")

    def test_check_point(self):
        point_x, point_y = -119.75, 39.5  # Example coordinates
        is_inside = check_point(self.boundary, point_x, point_y)
        self.assertIsInstance(is_inside, bool, "check_point should return a boolean")

if __name__ == '__main__':
    unittest.main()