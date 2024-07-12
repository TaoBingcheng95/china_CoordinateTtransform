"""
GPS(wgs84)、高德(gcj02)、百度(BD09)三种最常用的地图坐标系之间的相互转换。

未来计划实现
wgs84 Pseudo-Mercator投影转换
bd09 bd09mc投影转换
"""
import math

wkt_wgs84 = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],' \
            'AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",' \
            '0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]] '

wkt_pm = """PROJCS["WGS 84 / Pseudo-Mercator",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,
298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],
UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Mercator_1SP"],
PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",
0],UNIT["metre",1],AXIS["Easting",EAST],AXIS["Northing",NORTH],EXTENSION["PROJ4","+proj=merc +a=6378137 +b=6378137 
+lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs"],AUTHORITY["EPSG","3857"]] """

x_pi = 3.14159265358979324 * 3000.0 / 180.0
pi = 3.1415926535897932384626  # π
a = 6378245.0  # 长半轴
ee = 0.00669342162296594323  # 扁率

# 百度墨卡托投影纠正矩阵
LLBAND = [75, 60, 45, 30, 15, 0]
LL2MC = [
    [-0.0015702102444, 111320.7020616939, 1704480524535203, -10338987376042340, 26112667856603880, -35149669176653700,
     26595700718403920, -10725012454188240, 1800819912950474, 82.5],
    [0.0008277824516172526, 111320.7020463578, 647795574.6671607, -4082003173.641316, 10774905663.51142,
     -15171875531.51559, 12053065338.62167, -5124939663.577472, 913311935.9512032, 67.5],
    [0.00337398766765, 111320.7020202162, 4481351.045890365, -23393751.19931662, 79682215.47186455, -115964993.2797253,
     97236711.15602145, -43661946.33752821, 8477230.501135234, 52.5],
    [0.00220636496208, 111320.7020209128, 51751.86112841131, 3796837.749470245, 992013.7397791013, -1221952.21711287,
     1340652.697009075, -620943.6990984312, 144416.9293806241, 37.5],
    [-0.0003441963504368392, 111320.7020576856, 278.2353980772752, 2485758.690035394, 6070.750963243378,
     54821.18345352118, 9540.606633304236, -2710.55326746645, 1405.483844121726, 22.5],
    [-0.0003218135878613132, 111320.7020701615, 0.00369383431289, 823725.6402795718, 0.46104986909093,
     2351.343141331292, 1.58060784298199, 8.77738589078284, 0.37238884252424, 7.45]]

# 百度墨卡托转回到百度经纬度纠正矩阵
MCBAND = [12890594.86, 8362377.87, 5591021, 3481989.83, 1678043.12, 0]
MC2LL = [[1.410526172116255e-8, 0.00000898305509648872, -1.9939833816331, 200.9824383106796, -187.2403703815547,
          91.6087516669843, -23.38765649603339, 2.57121317296198, -0.03801003308653, 17337981.2],
         [-7.435856389565537e-9, 0.000008983055097726239, -0.78625201886289, 96.32687599759846, -1.85204757529826,
          -59.36935905485877, 47.40033549296737, -16.50741931063887, 2.28786674699375, 10260144.86],
         [-3.030883460898826e-8, 0.00000898305509983578, 0.30071316287616, 59.74293618442277, 7.357984074871,
          -25.38371002664745, 13.45380521110908, -3.29883767235584, 0.32710905363475, 6856817.37],
         [-1.981981304930552e-8, 0.000008983055099779535, 0.03278182852591, 40.31678527705744, 0.65659298677277,
          -4.44255534477492, 0.85341911805263, 0.12923347998204, -0.04625736007561, 4482777.06],
         [3.09191371068437e-9, 0.000008983055096812155, 0.00006995724062, 23.10934304144901, -0.00023663490511,
          -0.6321817810242, -0.00663494467273, 0.03430082397953, -0.00466043876332, 2555164.4],
         [2.890871144776878e-9, 0.000008983055095805407, -3.068298e-8, 7.47137025468032, -0.00000353937994,
          -0.02145144861037, -0.00001234426596, 0.00010322952773, -0.00000323890364, 826088.5]]


class LLT:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class CoordTrans(object):
    x_pi = 3.14159265358979324 * 3000.0 / 180.0
    pi = 3.1415926535897932384626  # π
    a = 6378245.0  # 长半轴
    ee = 0.00669342162296594323  # 偏心率平方
    rf = 1 / (1 - (1 - ee) ** 0.5)  #
    b = a - a / rf

    # ee = 1 - (1 - 1 / rf) ** 2       rf=298.3

    @classmethod
    def gcj02_to_bd09(cls, lng, lat):
        """
        火星坐标系(GCJ-02)转百度坐标系(BD-09)
        谷歌、高德——>百度
        :param lng:火星坐标经度
        :param lat:火星坐标纬度
        :return:
        """
        z = math.sqrt(lng * lng + lat * lat) + 0.00002 * math.sin(lat * cls.x_pi)
        theta = math.atan2(lat, lng) + 0.000003 * math.cos(lng * cls.x_pi)
        bd_lng = z * math.cos(theta) + 0.0065
        bd_lat = z * math.sin(theta) + 0.006
        return [bd_lng, bd_lat]

    @classmethod
    def bd09_to_gcj02(cls, bd_lon, bd_lat):
        """
        百度坐标系(BD-09)转火星坐标系(GCJ-02)
        百度——>谷歌、高德
        :param bd_lat:百度坐标纬度
        :param bd_lon:百度坐标经度
        :return:转换后的坐标列表形式
        """
        x = bd_lon - 0.0065
        y = bd_lat - 0.006
        z = math.sqrt(x * x + y * y) - 0.00002 * math.sin(y * cls.x_pi)
        theta = math.atan2(y, x) - 0.000003 * math.cos(x * cls.x_pi)
        gg_lng = z * math.cos(theta)
        gg_lat = z * math.sin(theta)
        return [gg_lng, gg_lat]

    @classmethod
    def wgs84_to_gcj02(cls, lng, lat):
        """
        WGS84转GCJ02(火星坐标系)
        :param lng:WGS84坐标系的经度
        :param lat:WGS84坐标系的纬度
        :return:
        """
        if cls.out_of_china(lng, lat):  # 判断是否在国内
            return [lng, lat]
        dlat = cls._transformlat(lng - 105.0, lat - 35.0)
        dlng = cls._transformlng(lng - 105.0, lat - 35.0)
        radlat = lat / 180.0 * cls.pi
        magic = math.sin(radlat)
        magic = 1 - cls.ee * magic * magic
        sqrtmagic = math.sqrt(magic)
        dlat = (dlat * 180.0) / ((cls.a * (1 - cls.ee)) / (magic * sqrtmagic) * cls.pi)
        dlng = (dlng * 180.0) / (cls.a / sqrtmagic * math.cos(radlat) * cls.pi)
        mglat = lat + dlat
        mglng = lng + dlng
        return [mglng, mglat]

    @classmethod
    def gcj02_to_wgs84(cls, lng, lat):
        """
        GCJ02(火星坐标系)转GPS84
        :param lng:火星坐标系的经度
        :param lat:火星坐标系纬度
        :return:
        """
        if cls.out_of_china(lng, lat):
            return [lng, lat]
        dlat = cls._transformlat(lng - 105.0, lat - 35.0)
        dlng = cls._transformlng(lng - 105.0, lat - 35.0)
        radlat = lat / 180.0 * cls.pi
        magic = math.sin(radlat)
        magic = 1 - cls.ee * magic * magic
        sqrtmagic = math.sqrt(magic)
        dlat = (dlat * 180.0) / ((cls.a * (1 - cls.ee)) / (magic * sqrtmagic) * cls.pi)
        dlng = (dlng * 180.0) / (cls.a / sqrtmagic * math.cos(radlat) * cls.pi)
        mglat = lat + dlat
        mglng = lng + dlng
        return [lng * 2 - mglng, lat * 2 - mglat]

    @classmethod
    def bd09_to_wgs84(cls, bd_lon, bd_lat):
        lon, lat = cls.bd09_to_gcj02(bd_lon, bd_lat)
        return cls.gcj02_to_wgs84(lon, lat)

    @classmethod
    def wgs84_to_bd09(cls, lon, lat):
        lon, lat = cls.wgs84_to_gcj02(lon, lat)
        return cls.gcj02_to_bd09(lon, lat)

    @classmethod
    def _transformlat(cls, lng, lat):
        pi = cls.pi
        ret = -100.0 + 2.0 * lng + 3.0 * lat + 0.2 * lat * lat + 0.1 * lng * lat + 0.2 * math.sqrt(math.fabs(lng))
        ret += (20.0 * math.sin(6.0 * lng * pi) + 20.0 * math.sin(2.0 * lng * pi)) * 2.0 / 3.0
        ret += (20.0 * math.sin(lat * pi) + 40.0 * math.sin(lat / 3.0 * pi)) * 2.0 / 3.0
        ret += (160.0 * math.sin(lat / 12.0 * pi) + 320 * math.sin(lat * pi / 30.0)) * 2.0 / 3.0
        return ret

    @classmethod
    def _transformlng(cls, lng, lat):
        pi = cls.pi
        ret = 300.0 + lng + 2.0 * lat + 0.1 * lng * lng + 0.1 * lng * lat + 0.1 * math.sqrt(math.fabs(lng))
        ret += (20.0 * math.sin(6.0 * lng * pi) + 20.0 * math.sin(2.0 * lng * pi)) * 2.0 / 3.0
        ret += (20.0 * math.sin(lng * pi) + 40.0 * math.sin(lng / 3.0 * pi)) * 2.0 / 3.0
        ret += (150.0 * math.sin(lng / 12.0 * pi) + 300.0 * math.sin(lng / 30.0 * pi)) * 2.0 / 3.0
        return ret

    @staticmethod
    def out_of_china(lng, lat):
        """
        判断是否在国内，不在国内不做偏移
        :param lng:
        :param lat:
        :return:
        """
        return not (73.66 < lng < 135.05 and 3.86 < lat < 53.55)

    @staticmethod
    def wgs84tomercator(lng, lat):
        """
        wgs84投影到墨卡托
        :param lng:
        :param lat:
        :return:
        """
        x = lng * 20037508.34 / 180
        y = math.log(math.tan((90 + lat) * math.pi / 360)) / (math.pi / 180) * 20037508.34 / 180
        return x, y

    @staticmethod
    def mercatortowgs84(x, y):
        """
        墨卡托投影坐标转回wgs84
        :param x:
        :param y:
        :return:
        """
        lng = x / 20037508.34 * 180
        lat = 180 / math.pi * (2 * math.atan(math.exp(y / 20037508.34 * 180 * math.pi / 180)) - math.pi / 2)
        return lng, lat

    @staticmethod
    def getRange(cC, cB, T):
        if cB is not None:
            cC = max(cC, cB)
        if T is not None:
            cC = min(cC, T)
        return cC

    @staticmethod
    def getLoop(cC, cB, T):
        while cC > T:
            cC -= T - cB
        while cC < cB:
            cC += T - cB
        return cC

    @staticmethod
    def convertor(cC, cD):
        if cC is None or cD is None:
            print('null')
            return None
        T = cD[0] + cD[1] * abs(cC.x)
        cB = abs(cC.y) / cD[9]
        cE = cD[2] + cD[3] * cB + cD[4] * cB * cB + cD[5] * cB * cB * cB + cD[6] * cB * cB * cB * cB + cD[
            7] * cB * cB * cB * cB * cB + cD[8] * cB * cB * cB * cB * cB * cB
        if cC.x < 0:
            T = T * -1
        else:
            T = T
        if cC.y < 0:
            cE = cE * -1
        else:
            cE = cE
        return [T, cE]

    @staticmethod
    def convertLL2MC(T):
        cD = None
        T.x = getLoop(T.x, -180, 180)
        T.y = getRange(T.y, -74, 74)
        cB = T
        for cC in range(0, len(LLBAND), 1):
            if cB.y >= LLBAND[cC]:
                cD = LL2MC[cC]
                break
        if cD is not None:
            for cC in range(len(LLBAND) - 1, -1, -1):
                if cB.y <= -LLBAND[cC]:
                    cD = LL2MC[cC]
                    break
        cE = convertor(T, cD)
        return cE

    @staticmethod
    def convertMC2LL(cB):
        cC = LLT(abs(cB.x), abs(cB.y))
        cE = None
        for cD in range(0, len(MCBAND), 1):
            if cC.y >= MCBAND[cD]:
                cE = MC2LL[cD]
                break
        T = convertor(cB, cE)
        return T

    @staticmethod
    def bd09tomercator(lng, lat):
        """
        bd09投影到百度墨卡托
        :param lng:
        :param lat:
        :return:
        """
        baidut = LLT(lng, lat)
        return convertLL2MC(baidut)

    @staticmethod
    def mercatortobd09(x, y):
        """
        墨卡托投影坐标转回bd09
        :param x:
        :param y:
        :return:
        """
        baidut = LLT(x, y)
        return convertMC2LL(baidut)


def gcj02tobd09(lng, lat):
    """
    火星坐标系(GCJ02)转百度坐标系(BD09)
    :param lng:火星坐标经度
    :param lat:火星坐标纬度
    :return:
    """
    z = math.sqrt(lng * lng + lat * lat) + 0.00002 * math.sin(lat * x_pi)
    theta = math.atan2(lat, lng) + 0.000003 * math.cos(lng * x_pi)
    bd_lng = z * math.cos(theta) + 0.0065
    bd_lat = z * math.sin(theta) + 0.006
    return [bd_lng, bd_lat]


def bd09togcj02(bd_lon, bd_lat):
    """
    百度坐标系(BD09)转火星坐标系(GCJ02)
    :param bd_lat:百度坐标纬度
    :param bd_lon:百度坐标经度
    :return:转换后的坐标列表形式
    """
    x = bd_lon - 0.0065
    y = bd_lat - 0.006
    z = math.sqrt(x * x + y * y) - 0.00002 * math.sin(y * x_pi)
    theta = math.atan2(y, x) - 0.000003 * math.cos(x * x_pi)
    gg_lng = z * math.cos(theta)
    gg_lat = z * math.sin(theta)
    return [gg_lng, gg_lat]


def wgs84togcj02(lng, lat):
    """
    WGS84转GCJ02(火星坐标系)
    :param lng:WGS84坐标系的经度
    :param lat:WGS84坐标系的纬度
    :return:
    """
    if out_of_china(lng, lat):  # 判断是否在国内
        return lng, lat
    dlat = transformlat(lng - 105.0, lat - 35.0)
    dlng = transformlng(lng - 105.0, lat - 35.0)
    radlat = lat / 180.0 * pi
    magic = math.sin(radlat)
    magic = 1 - ee * magic * magic
    sqrtmagic = math.sqrt(magic)
    dlat = (dlat * 180.0) / ((a * (1 - ee)) / (magic * sqrtmagic) * pi)
    dlng = (dlng * 180.0) / (a / sqrtmagic * math.cos(radlat) * pi)
    mglat = lat + dlat
    mglng = lng + dlng
    return [mglng, mglat]


def gcj02towgs84(lng, lat):
    """
    GCJ02(火星坐标系)转GPS84
    :param lng:火星坐标系的经度
    :param lat:火星坐标系纬度
    :return:
    """
    if out_of_china(lng, lat):
        return lng, lat
    dlat = transformlat(lng - 105.0, lat - 35.0)
    dlng = transformlng(lng - 105.0, lat - 35.0)
    radlat = lat / 180.0 * pi
    magic = math.sin(radlat)
    magic = 1 - ee * magic * magic
    sqrtmagic = math.sqrt(magic)
    dlat = (dlat * 180.0) / ((a * (1 - ee)) / (magic * sqrtmagic) * pi)
    dlng = (dlng * 180.0) / (a / sqrtmagic * math.cos(radlat) * pi)
    mglat = lat + dlat
    mglng = lng + dlng
    return [lng * 2 - mglng, lat * 2 - mglat]


def transformlat(lng, lat):
    ret = -100.0 + 2.0 * lng + 3.0 * lat + 0.2 * lat * lat + 0.1 * lng * lat + 0.2 * math.sqrt(math.fabs(lng))
    ret += (20.0 * math.sin(6.0 * lng * pi) + 20.0 * math.sin(2.0 * lng * pi)) * 2.0 / 3.0
    ret += (20.0 * math.sin(lat * pi) + 40.0 *
            math.sin(lat / 3.0 * pi)) * 2.0 / 3.0
    ret += (160.0 * math.sin(lat / 12.0 * pi) + 320 *
            math.sin(lat * pi / 30.0)) * 2.0 / 3.0
    return ret


def transformlng(lng, lat):
    ret = 300.0 + lng + 2.0 * lat + 0.1 * lng * lng + 0.1 * lng * lat + 0.1 * math.sqrt(math.fabs(lng))
    ret += (20.0 * math.sin(6.0 * lng * pi) + 20.0 * math.sin(2.0 * lng * pi)) * 2.0 / 3.0
    ret += (20.0 * math.sin(lng * pi) + 40.0 * math.sin(lng / 3.0 * pi)) * 2.0 / 3.0
    ret += (150.0 * math.sin(lng / 12.0 * pi) + 300.0 * math.sin(lng / 30.0 * pi)) * 2.0 / 3.0
    return ret


def out_of_china(lng, lat):
    """
    判断是否在国内，不在国内不做偏移
    :param lng:
    :param lat:
    :return:
    """
    if lng < 72.004 or lng > 137.8347:
        return True
    if lat < 0.8293 or lat > 55.8271:
        return True
    return False


def wgs84tomercator(lng, lat):
    """
    wgs84投影到墨卡托
    :param lng:
    :param lat:
    :return:
    """
    x = lng * 20037508.34 / 180
    y = math.log(math.tan((90 + lat) * math.pi / 360)) / (math.pi / 180) * 20037508.34 / 180
    return x, y


def mercatortowgs84(x, y):
    """
    墨卡托投影坐标转回wgs84
    :param x:
    :param y:
    :return:
    """
    lng = x / 20037508.34 * 180
    lat = 180 / math.pi * (2 * math.atan(math.exp(y / 20037508.34 * 180 * math.pi / 180)) - math.pi / 2)
    return lng, lat


def getRange(cC, cB, T):
    if cB is not None:
        cC = max(cC, cB)
    if T is not None:
        cC = min(cC, T)
    return cC


def getLoop(cC, cB, T):
    while cC > T:
        cC -= T - cB
    while cC < cB:
        cC += T - cB
    return cC


def convertor(cC, cD):
    if cC is None or cD is None:
        print('null')
        return None
    T = cD[0] + cD[1] * abs(cC.x)
    cB = abs(cC.y) / cD[9]
    cE = cD[2] + cD[3] * cB + cD[4] * cB * cB + cD[5] * cB * cB * cB + cD[6] * cB * cB * cB * cB + cD[
        7] * cB * cB * cB * cB * cB + cD[8] * cB * cB * cB * cB * cB * cB
    if cC.x < 0:
        T = T * -1
    else:
        T = T
    if cC.y < 0:
        cE = cE * -1
    else:
        cE = cE
    return [T, cE]


def convertLL2MC(T):
    cD = None
    T.x = getLoop(T.x, -180, 180)
    T.y = getRange(T.y, -74, 74)
    cB = T
    for cC in range(0, len(LLBAND), 1):
        if cB.y >= LLBAND[cC]:
            cD = LL2MC[cC]
            break
    if cD is not None:
        for cC in range(len(LLBAND) - 1, -1, -1):
            if cB.y <= -LLBAND[cC]:
                cD = LL2MC[cC]
                break
    cE = convertor(T, cD)
    return cE


def convertMC2LL(cB):
    cC = LLT(abs(cB.x), abs(cB.y))
    cE = None
    for cD in range(0, len(MCBAND), 1):
        if cC.y >= MCBAND[cD]:
            cE = MC2LL[cD]
            break
    T = convertor(cB, cE)
    return T


def bd09tomercator(lng, lat):
    """
    bd09投影到百度墨卡托
    :param lng:
    :param lat:
    :return:
    """
    baidut = LLT(lng, lat)
    return convertLL2MC(baidut)


def mercatortobd09(x, y):
    """
    墨卡托投影坐标转回bd09
    :param x:
    :param y:
    :return:
    """
    baidut = LLT(x, y)
    return convertMC2LL(baidut)


if __name__ == '__main__':
    lng = 128.543
    lat = 37.065
    result1 = CoordTrans.gcj02_to_bd09(lng, lat)
    result2 = CoordTrans.bd09_to_gcj02(lng, lat)
    result3 = CoordTrans.wgs84_to_gcj02(lng, lat)
    result4 = CoordTrans.gcj02_to_wgs84(lng, lat)
    result5 = CoordTrans.bd09_to_wgs84(lng, lat)
    result6 = CoordTrans.wgs84_to_bd09(lng, lat)

    print("gcj02_to_bd09 : ", result1)
    print("bd09_to_gcj02 : ", result2)
    print("wgs84_to_gcj02 : ", result3)
    print("gcj02_to_wgs84 : ", result4)
    print("bd09_to_wgs84 : ", result5)
    print("wgs84_to_bd09 : ", result6)
