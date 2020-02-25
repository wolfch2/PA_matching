# https://stackoverflow.com/questions/34480328/numpy-summing-up-a-list-of-vectors
# https://stackoverflow.com/questions/53093073/why-is-numpy-sum-10-times-slower-than-the-operator
# https://stackoverflow.com/questions/34618301/python-multiprocessing-using-pickle-kwargs-and-function-references
# geopandas tries to load entire geometry into a col. I think... very slow and memory intensive
# https://stackoverflow.com/questions/10741346/numpy-most-efficient-frequency-counts-for-unique-values-in-an-array
# https://medium.com/distributed-computing-with-ray/ray-for-the-curious-fa0e019e17d3

from osgeo import gdal, ogr
import pandas as pd
import fiona
import shapely
import rasterio
from pyproj import Proj, transform
from rasterio import features
from rasterio.features import bounds as calculate_bounds
import geojson
from fiona import transform
import math
import affine
import multiprocessing
import glob
import functools
import operator
from dbfread import DBF
import os
import random
import numpy as np
import matplotlib.pyplot as plt
import ctypes
from collections import defaultdict
import psutil
import ray
import time

@ray.remote
class MyActor(object):
    def __init__(self, elev_mat_id):
        self.all_features = fiona.open("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/range_maps/MAMMALS.shp")
        self.elev = rasterio.open('/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/elev.tif')
        self.elev_mat = elev_mat_id
        self.richness = (0 * self.elev_mat).astype('int16')
    def add_range(self, i):
        print(i)
        shapes = [self.all_features[index]["geometry"] for index in [i]]
        range_rast = rasterio.features.rasterize(shapes,
                                        out_shape = self.elev.shape,
                                        transform = self.elev.transform,
                                        default_value= 1,
                                        all_touched=False)
        range_rast[np.logical_or(self.elev_mat < 3, self.elev_mat > 10)] = 0
        self.richness += range_rast
    def get_richness(self):
        return self.richness

elev = rasterio.open('/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/elev.tif')
elev_mat = np.round(elev.read(1) / 100).astype('int8')

ray.init(object_store_memory=20*1024*1024*1024)
elev_mat_id = ray.put(elev_mat) # https://stackoverflow.com/questions/54582073/sharing-objects-across-workers-using-pyarrow
actors = [MyActor.remote(elev_mat_id) for _ in range(12)]
pool = ray.experimental.ActorPool(actors)
list(pool.map_unordered(lambda actor, v: actor.add_range.remote(v), list(range(2000))))

results = ray.get([actor.get_richness.remote() for actor in actors])

richness = results[0]
for i in range(1,12):
    richness += results[i]

plt.matshow(richness)
plt.show()


