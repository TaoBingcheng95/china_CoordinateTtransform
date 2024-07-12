import os
from pathlib import Path
from glob import glob
from osgeo import gdal, osr, ogr

from china_transform.RasterTransform import transform_raster

if __name__ == "__main__":
    src_data = "test_data\do_2207.tif"
    dst_data = "test_data\do_2207_gcj02_new.tif"

    transform_raster(src_data, dst_data, "wgs84", "gcj02")
