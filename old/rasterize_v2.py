# https://stackoverflow.com/questions/34480328/numpy-summing-up-a-list-of-vectors
# https://stackoverflow.com/questions/53093073/why-is-numpy-sum-10-times-slower-than-the-operator
# https://stackoverflow.com/questions/34618301/python-multiprocessing-using-pickle-kwargs-and-function-references
# geopandas tries to load entire geometry into a col. I think... very slow and memory intensive

from osgeo import gdal, ogr
import pandas as pd
import fiona
import shapely
import rasterio
import pprint
from pyproj import Proj, transform
from rasterio import features
from rasterio.features import bounds as calculate_bounds
import geojson
from fiona import transform
import math
import affine
import multiprocess 
import glob
import functools
import operator
from dbfread import DBF
import os
import random
import numpy as np
import cupy as cp


all_features = fiona.open("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/range_maps/MAMMALS.shp")
shapes = [feature["geometry"] for feature in all_features[10:20]]

bounds = (-180, -90, 180, 90)
res = (0.01, 0.01)
width = max(int(math.ceil((bounds[2] - bounds[0]) / float(res[0]))), 1)
height = max(int(math.ceil((bounds[3] - bounds[1]) / float(res[1]))), 1)
trans = affine.Affine(res[0], 0, bounds[0], 0, -res[1], bounds[3])





def get_rast(i):
    print(i)        
    shapes = [feature["geometry"] for feature in all_features[10:20]]
    return rasterio.features.rasterize(shapes,
                                        out_shape = (height, width),
                                        transform = trans,
                                        default_value=1,
                                        all_touched=False)

test = list(map(get_rast, range(30)))
hmm = np.add.reduce(test)




stack = np.stack(test)
cp_stack = cp.array(stack)

test = cp_stack.sum(axis=0)

out = test[0]
for i in range(1,30):
        print(i)
        out = out + test[i]

a = cp.array(test[3])
b = cp.array(test[4])
c = a + b

c = test[3] + test[4]


rast = rasterio.open(
    "/home/chrisgraywolf/analysis_desktop/PA_matching/hmm.tif",
    'w',
    driver='GTiff',
    height=height,
    width=width,
    count=1,
    dtype=rasterio.uint8,
    crs='EPSG:4326',
    compress='lzw',
    transform=trans,
)

rast.write(test[0], 1)
rast.close()



as dst:
    dst.write(data




all_features = fiona.open("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp")

bounds = (-180, -90, 180, 90)
res = (0.01, 0.01)
width = max(int(math.ceil((bounds[2] - bounds[0]) / float(res[0]))), 1)
height = max(int(math.ceil((bounds[3] - bounds[1]) / float(res[1]))), 1)
trans = affine.Affine(res[0], 0, bounds[0], 0, -res[1], bounds[3])

# https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html
shapes = [feature["geometry"] for feature in all_features[0:100000]]

rast = rasterio.features.rasterize(shapes,
                                        out_shape = (height, width),
                                        transform = trans,
                                        default_value=1,
                                        all_touched=False) # too slow.. sticking w/ fasterize






select_features = list(map(lambda x: all_features[x], range(100)))

select_features_geom = shapely.geometry.shape(geojson.GeometryCollection(select_features))





[shapely.geometry.shape(select_features[1]['geometry'])]







output_name = "/home/stats/wolfch/Tom2/temp/"+str(id_no)
    if os.path.isfile(output_name):
        return pd.read_csv(output_name).to_dict('index')[0]
    print file_name, id_no, select_indices

    select_features = map(lambda(x): all_features[x], select_indices)
    Behr_features = map(lambda(x): transform.transform_geom('EPSG:4326', Behr, x['geometry']), select_features)
    Behr_features_gc = shapely.geometry.shape(geojson.GeometryCollection(Behr_features))
    # https://github.com/mapbox/rasterio/blob/master/rasterio/rio/rasterize.py
    bounds = Behr_features_gc.bounds
    bounds = (bounds[0] - 100e3, bounds[1] - 100e3, bounds[2] + 100e3, bounds[3] + 100e3) # weird bug? https://trac.osgeo.org/gdal/ticket/7176
    # ^-- note wgs84_bounding_box rasterizes fine with this line (could clip extents o.w.).  I get 513 million km2 for world area (Google says 510)
    # perhaps bug is related to shapefile edge touching outer edge or raster border or something?
    res = (100e3, 100e3)
    width = max(int(math.ceil((bounds[2] - bounds[0]) / float(res[0]))), 1)
    height = max(int(math.ceil((bounds[3] - bounds[1]) / float(res[1]))), 1)
    trans = affine.Affine(res[0], 0, bounds[0], 0, -res[1], bounds[3])
    #
    num_cells = rasterio.features.rasterize(Behr_features_gc,
                                        out_shape = (height, width),
                                        transform = trans,
                                        default_value=1,
                                        all_touched=True).sum()
    area_km2 = num_cells * 1e4
    out = {'area_km2':area_km2, 'id_no':id_no}
    pd.DataFrame(out,index=[0]).to_csv(output_name,index=False)





output_name = "/home/stats/wolfch/Tom2/temp/"+str(id_no)
    if os.path.isfile(output_name):
        return pd.read_csv(output_name).to_dict('index')[0]
    print file_name, id_no, select_indices
    all_features = fiona.open(file_name) # can't pickle even w/ multiprocess instead of multiprocessing
    select_features = map(lambda(x): all_features[x], select_indices)
    Behr_features = map(lambda(x): transform.transform_geom('EPSG:4326', Behr, x['geometry']), select_features)
    Behr_features_gc = shapely.geometry.shape(geojson.GeometryCollection(Behr_features))
    # https://github.com/mapbox/rasterio/blob/master/rasterio/rio/rasterize.py
    bounds = Behr_features_gc.bounds
    bounds = (bounds[0] - 100e3, bounds[1] - 100e3, bounds[2] + 100e3, bounds[3] + 100e3) # weird bug? https://trac.osgeo.org/gdal/ticket/7176
    # ^-- note wgs84_bounding_box rasterizes fine with this line (could clip extents o.w.).  I get 513 million km2 for world area (Google says 510)
    # perhaps bug is related to shapefile edge touching outer edge or raster border or something?
    res = (100e3, 100e3)
    width = max(int(math.ceil((bounds[2] - bounds[0]) / float(res[0]))), 1)
    height = max(int(math.ceil((bounds[3] - bounds[1]) / float(res[1]))), 1)
    trans = affine.Affine(res[0], 0, bounds[0], 0, -res[1], bounds[3])
    #
    num_cells = rasterio.features.rasterize(Behr_features_gc,
                                        out_shape = (height, width),
                                        transform = trans,
                                        default_value=1,
                                        all_touched=True).sum()
    area_km2 = num_cells * 1e4
    out = {'area_km2':area_km2, 'id_no':id_no}
    pd.DataFrame(out,index=[0]).to_csv(output_name,index=False)








require(sf)
require(fasterize)
require(raster)

setwd("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/PAs/")

PAs = read_sf("WDPA_Jan2020-shapefile-polygons.shp")

r = raster(resolution=c(0.01,0.01))

start = Sys.time()
test = fasterize(PAs, r)
Sys.time() - start # 4.7 sec :O



