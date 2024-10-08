# china_coordinate_transform

## 基本介绍

### 一、谷歌、高德、百度经纬坐标相互转换（CoordTransform模块）

+ 无需第三方包

+ GPS（wgs84）、高德（gcj02）、百度（BD09）三种最常用的地图坐标系之间的相互转换。

1. WGS84原始坐标系，一般用国际GPS纪录仪记录下来的经纬度，通过GPS定位拿到的原始经纬度，Google和高德地图定位的的经纬度（国外）都是基于WGS84坐标系的；但是在国内是不允许直接用WGS84坐标系标注的，必须经过加密后才能使用；

2. GCJ02坐标系，又名“火星坐标系”，是我国国测局独创的坐标体系，由WGS84加密而成，在国内，必须至少使用GCJ－02坐标系，或者使用在GCJ02加密后再进行加密的坐标系，如百度坐标系。高德和Google在国内都是使用GCJ02坐标系，可以说GCJ02是国内最广泛使用的坐标系；

3. 百度坐标系:bd-09，百度坐标系是在GCJ－02坐标系的基础上再次加密偏移后形成的坐标系，只适用于百度地图。(目前百度API提供了从其它坐标系转换为百度坐标系的API，但却没有从百度坐标系转为其他坐标系的API);

**注意！转换公式是近似计算的，实际公式因为安全等原因未公开，wgs84 和 gcj02 不是直接互相转换的。模块中公式在不同经纬度地区误差不同，但不影响一般实际的手机导航使用。**

以此模块为核心，衍生出针对矢量文件转换的VectorTransform模块和针对栅格文件转换的RasterTransform模块，可用于对项目中的地理数据进行处理。

### 二、中国常用大地测量投影坐标系转换（ChinaCoordProjection 模块）

* 需要安装 pyproj 包。 **pip install pyproj** 

+ 中国常用大地测量（投影）坐标系转换（标准转换，不涉及使用7参数平移、旋转、缩放）。

+ 主要有wgs84、西安80、北京54、新北京、cgcs2000等我国常用坐标系之间的相互转换。

+ 本模块是对 pyproj 中的相关方法进行了二次封装，以方便使用！采用epsg中记录的各个坐标系参数

+ pyproj 官方文档：http://pyproj4.github.io/pyproj/stable/

**note**：

+ pyproj 通过 crs 生成transformer 是非常耗时的操作，经测试AMD Ryzen 5 2500U with Radeon Vega Mobile 2.00GHz平台，创建一个需要0.2~0.3秒。在批量转换的任务中，如果对每个点生成一个transformer，效率是不能接受的！考虑到，一个确定的转换任务中，需要重新生成transformer的地方在：输入点位于不同投影分度带。

+ 这里采用的策略是参考单例设计模式的思想，在批量转换的时候，先获取某个点转换的特征：

        1. 输入坐标projection
        
        2. 输入坐标中央经线
        
        3. 输入坐标是否带有分度带带号
        
        4. 输出坐标的projection
        
        5. 输出坐标中央经线
        
        6. 输出坐标是否带有分度带带号

+ 根据这些特征可以唯一的描述一个transformer，可以根据这些特征，生成一个transformer_key，用transformer_key作为键将transformer对象保存至字典中。在生成transformer之前检查拥有此特征的transformer是否已经生成，如果已经生成则直接返回，如果未生成，则生成。

+ 因为每次批量转换的任务所涉及的范围一般不会跨几个分度带，因此不会生成特别多的transformer，所以空间占用和整体的时间效率都能接受。

## 三、高德地图地理编码和逆地理编码 web api (amap 模块)

* 需要安装 requests 包。 **pip install requests**

1. 高德地图地理编码和逆地理编码 web api 的 python封装，以方便使用。

2. 从高德地图 js 库中查出的两点经纬度坐标距离计算（简单的计算方法，计算两经纬度坐标之间的大圆弧长，主流的在线地图测距工具算法，精确度不高，常用于生活场景）。

高德地图地理编码和逆地理编码 web api文档：[https://lbs.amap.com/api/webservice/guide/api/georegeo](https://lbs.amap.com/api/webservice/guide/api/georegeo)