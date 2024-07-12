import os
# from glob import glob
# from shutil import copyfile
import numpy as np
import math
from osgeo import gdal, osr, gdal_array, gdalconst
from .CoordTransform import gcj02tobd09, bd09togcj02, wgs84togcj02, gcj02towgs84, wgs84tomercator, mercatortowgs84, \
    bd09tomercator, mercatortobd09, wkt_wgs84
from .CoordTransform import CoordTrans


class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y


def affine_abc(pixel_points, geo_points):
    n = len(pixel_points)
    pixelX_square = 0.0
    pixelX_pixelY = 0.0
    pixelY_square = 0.0
    pixelX = 0.0
    pixelY = 0.0
    pixelX_geoX = 0.0
    pixelY_geoX = 0.0
    geoX = 0.0
    for i in range(0, n):
        pixelX_square += math.pow(pixel_points[i].x, 2)
        pixelX_pixelY += pixel_points[i].x * pixel_points[i].y
        pixelY_square += math.pow(pixel_points[i].y, 2)
        pixelX += pixel_points[i].x
        pixelY += pixel_points[i].y
        pixelX_geoX += pixel_points[i].x * geo_points[i].x
        pixelY_geoX += pixel_points[i].y * geo_points[i].x
        geoX += geo_points[i].x
    a = np.array([[pixelX_square, pixelX_pixelY, pixelX], [pixelX_pixelY, pixelY_square, pixelY], [pixelX, pixelY, n]])
    b = np.array([[pixelX_geoX], [pixelY_geoX], [geoX]])
    at = np.linalg.inv(a)
    result = at.dot(b)
    return result[0, 0], result[1, 0], result[2, 0]


def affine_def(pixel_points, geo_points):
    n = len(pixel_points)
    pixelX_square = 0.0
    pixelX_pixelY = 0.0
    pixelY_square = 0.0
    pixelX = 0.0
    pixelY = 0.0
    pixelX_geoY = 0.0
    pixelY_geoY = 0.0
    geoY = 0.0
    for i in range(0, n):
        pixelX_square += math.pow(pixel_points[i].x, 2)
        pixelX_pixelY += pixel_points[i].x * pixel_points[i].y
        pixelY_square += math.pow(pixel_points[i].y, 2)
        pixelX += pixel_points[i].x
        pixelY += pixel_points[i].y
        pixelX_geoY += pixel_points[i].x * geo_points[i].y
        pixelY_geoY += pixel_points[i].y * geo_points[i].y
        geoY += geo_points[i].y
    a = np.array([[pixelX_square, pixelX_pixelY, pixelX], [pixelX_pixelY, pixelY_square, pixelY], [pixelX, pixelY, n]])
    b = np.array([[pixelX_geoY], [pixelY_geoY], [geoY]])
    at = np.linalg.inv(a)
    result = at.dot(b)
    return result[0, 0], result[1, 0], result[2, 0]


# def abc_def(path_str):
#     pfw = open(path_str, 'r')
#     affineOption = pfw.readlines()
#     pfw.close()
#     a = float(affineOption[0].strip('\n'))
#     d = float(affineOption[1].strip('\n'))
#     b = float(affineOption[2].strip('\n'))
#     e = float(affineOption[3].strip('\n'))
#     c = float(affineOption[4].strip('\n'))
#     f = float(affineOption[5].strip('\n'))
#     return a, b, c, d, e, f


def generate_worldfile(world_file, geo_transform):
    ul_x = geo_transform[0] + geo_transform[1] / 2 + geo_transform[2] / 2
    ul_y = geo_transform[3] + geo_transform[4] / 2 + geo_transform[5] / 2
    with open(world_file, 'wt') as tfw_f:  # os.path.join(psz_dir, f"{psz_stem}.tfw")
        tfw_f.write("%0.8f\n" % geo_transform[1])
        tfw_f.write("%0.8f\n" % geo_transform[2])
        tfw_f.write("%0.8f\n" % geo_transform[4])
        tfw_f.write("%0.8f\n" % geo_transform[5])
        tfw_f.write("%0.8f\n" % ul_x)
        tfw_f.write("%0.8f\n" % ul_y)
    return None


def generate_prj(prj_file, projection):
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(projection)
    src_srs.MorphToESRI()
    src_wkt = src_srs.ExportToWkt()
    with open(prj_file, 'wt') as prj_f:  # os.path.join(psz_dir, f"{psz_stem}.prj")
        prj_f.write(src_wkt)
    del src_srs
    return None


def transform_raster(psz_file, dst_file, coord_original, coord_target, associate_file=False):
    coord_original = coord_original.lower()
    coord_target = coord_target.lower()
    if coord_original not in ['gcj02', 'wgs84', 'bd09']:
        print("coord_original value must be one of gcj02, wgs84 or bd09")
        return False
    if coord_target not in ['gcj02', 'wgs84', 'bd09']:
        print("coord_target value must be one of gcj02, wgs84 or bd09")
        return False

    try:
        ds = gdal.Open(psz_file, gdalconst.GA_ReadOnly)
    except RuntimeError as e:
        print(e)
        return False
    if ds is None:
        print(f"{psz_file} : Can not open such file or directory")
        return False

    psz_dir = os.path.dirname(psz_file)
    psz_stem = os.path.splitext(os.path.basename(psz_file))[0]

    im_driver = ds.GetDriver()
    im_width = ds.RasterXSize
    im_height = ds.RasterYSize
    # im_band = ds.RasterCount
    im_GeoTransform = ds.GetGeoTransform()
    im_proj = ds.GetProjection()
    if associate_file:
        if im_proj != "":
            prj_file = os.path.join(psz_dir, f"{psz_stem}.prj")
            generate_prj(prj_file, im_proj)
        world_file = os.path.join(psz_dir, f"{psz_stem}.tfw")
        generate_worldfile(world_file, im_GeoTransform)

    # GeoTranform -> World File
    # 0->4, 1->0, 2->1, 3->5, 4->2, 5->3
    # World File -> GeoTranform
    # 0->1, 1->2, 2->4, 3->5, 4->0, 5->3
    a = im_GeoTransform[1]
    b = im_GeoTransform[2]
    d = im_GeoTransform[4]
    e = im_GeoTransform[5]
    c = im_GeoTransform[0] + a / 2.0 + b / 2.0
    f = im_GeoTransform[3] + d / 2.0 + e / 2.0

    p0 = [im_width / 4, im_height / 4]
    p1 = [im_width / 4 + im_width / 2, im_height / 4]
    p2 = [im_width / 4 + im_width / 2, im_height / 4 + im_height / 2]
    p3 = [im_width / 4, im_height / 4 + im_height / 2]
    go0 = [a * p0[0] + b * p0[1] + c, d * p0[0] + e * p0[1] + f]
    go1 = [a * p1[0] + b * p1[1] + c, d * p1[0] + e * p1[1] + f]
    go2 = [a * p2[0] + b * p2[1] + c, d * p2[0] + e * p2[1] + f]
    go3 = [a * p3[0] + b * p3[1] + c, d * p3[0] + e * p3[1] + f]
    # 坐标系转换
    if coord_original == 'gcj02' and coord_target == 'wgs84':
        gt0 = CoordTrans.gcj02_to_wgs84(go0[0], go0[1])
        gt1 = CoordTrans.gcj02_to_wgs84(go1[0], go1[1])
        gt2 = CoordTrans.gcj02_to_wgs84(go2[0], go2[1])
        gt3 = CoordTrans.gcj02_to_wgs84(go3[0], go3[1])
    elif coord_original == 'wgs84' and coord_target == 'gcj02':
        gt0 = CoordTrans.wgs84_to_gcj02(go0[0], go0[1])
        gt1 = CoordTrans.wgs84_to_gcj02(go1[0], go1[1])
        gt2 = CoordTrans.wgs84_to_gcj02(go2[0], go2[1])
        gt3 = CoordTrans.wgs84_to_gcj02(go3[0], go3[1])
    elif coord_original == 'gcj02' and coord_target == 'bd09':
        gt0 = CoordTrans.gcj02_to_bd09(go0[0], go0[1])
        gt1 = CoordTrans.gcj02_to_bd09(go1[0], go1[1])
        gt2 = CoordTrans.gcj02_to_bd09(go2[0], go2[1])
        gt3 = CoordTrans.gcj02_to_bd09(go3[0], go3[1])
    elif coord_original == 'bd09' and coord_target == 'gcj02':
        gt0 = CoordTrans.bd09_to_gcj02(go0[0], go0[1])
        gt1 = CoordTrans.bd09_to_gcj02(go1[0], go1[1])
        gt2 = CoordTrans.bd09_to_gcj02(go2[0], go2[1])
        gt3 = CoordTrans.bd09_to_gcj02(go3[0], go3[1])
    elif coord_original == 'wgs84' and coord_target == 'bd09':
        gtm0 = CoordTrans.wgs84_to_gcj02(go0[0], go0[1])
        gtm1 = CoordTrans.wgs84_to_gcj02(go1[0], go1[1])
        gtm2 = CoordTrans.wgs84_to_gcj02(go2[0], go2[1])
        gtm3 = CoordTrans.wgs84_to_gcj02(go3[0], go3[1])
        gt0 = CoordTrans.gcj02_to_bd09(gtm0[0], gtm0[1])
        gt1 = CoordTrans.gcj02_to_bd09(gtm1[0], gtm1[1])
        gt2 = CoordTrans.gcj02_to_bd09(gtm2[0], gtm2[1])
        gt3 = CoordTrans.gcj02_to_bd09(gtm3[0], gtm3[1])
    elif coord_original == 'bd09' and coord_target == 'wgs84':
        gtm0 = CoordTrans.bd09_to_gcj02(go0[0], go0[1])
        gtm1 = CoordTrans.bd09_to_gcj02(go1[0], go1[1])
        gtm2 = CoordTrans.bd09_to_gcj02(go2[0], go2[1])
        gtm3 = CoordTrans.bd09_to_gcj02(go3[0], go3[1])
        gt0 = CoordTrans.gcj02_to_wgs84(gtm0[0], gtm0[1])
        gt1 = CoordTrans.gcj02_to_wgs84(gtm1[0], gtm1[1])
        gt2 = CoordTrans.gcj02_to_wgs84(gtm2[0], gtm2[1])
        gt3 = CoordTrans.gcj02_to_wgs84(gtm3[0], gtm3[1])
    else:
        gt0 = [go0[0], go0[1]]
        gt1 = [go1[0], go1[1]]
        gt2 = [go2[0], go2[1]]
        gt3 = [go3[0], go3[1]]
    pl = [Point(p0[0], p0[1]), Point(p1[0], p1[1]), Point(p2[0], p2[1]), Point(p3[0], p3[1])]
    gl = [Point(gt0[0], gt0[1]), Point(gt1[0], gt1[1]), Point(gt2[0], gt2[1]), Point(gt3[0], gt3[1])]
    ar, br, cr = affine_abc(pl, gl)
    dr, er, fr = affine_def(pl, gl)

    new_geotransform = [0, 1, 0, 0, 0, 1]  # default
    new_geotransform[0] = cr - ar / 2.0 - br / 2.0
    new_geotransform[1] = ar
    new_geotransform[2] = br
    new_geotransform[3] = fr - dr / 2.0 - er / 2.0
    new_geotransform[4] = dr
    new_geotransform[5] = er

    dst_dir = os.path.dirname(dst_file)
    dst_stem = os.path.splitext(os.path.basename(dst_file))[0]
    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    options = ["TFW=YES"]  #
    dst_ds = im_driver.CreateCopy(str(dst_file), ds, 0, options=options)
    if dst_ds.GetDriver().ShortName == 'GTiff':
        dst_ds.SetGeoTransform([0, 0, 0, 0, 0, 0])
    else:
        dst_ds.SetGeoTransform([0, 1, 0, 0, 0, 1])
    dst_ds.SetGeoTransform(new_geotransform)
    if coord_target in ['gcj02', 'bd09']:
        dst_ds.SetProjection(wkt_wgs84)  # for gcj02 or bd09
        dst_proj = wkt_wgs84
    else:
        dst_proj = dst_ds.GetGeotransform()

    dst_prj_file = os.path.join(dst_dir, f"{dst_stem}.prj")
    generate_prj(dst_prj_file, dst_proj)

    del ds
    del dst_ds

    return None


if __name__ == "__main__":
    src_data = "chla_2207.tif"
    dst_data = ".chla_2207_gcj02.tif"
    transform_raster(src_data, dst_data, "wgs84", "gcj02", associate_file=False)
