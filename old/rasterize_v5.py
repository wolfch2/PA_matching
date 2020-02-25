# https://stackoverflow.com/questions/34480328/numpy-summing-up-a-list-of-vectors
# https://stackoverflow.com/questions/53093073/why-is-numpy-sum-10-times-slower-than-the-operator
# https://stackoverflow.com/questions/34618301/python-multiprocessing-using-pickle-kwargs-and-function-references
# geopandas tries to load entire geometry into a col. I think... very slow and memory intensive
# https://stackoverflow.com/questions/10741346/numpy-most-efficient-frequency-counts-for-unique-values-in-an-array

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

all_features = fiona.open("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/range_maps/MAMMALS.shp")
elev = rasterio.open('/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/elev.tif')
elev_mat = elev.read(1)

richness = multiprocessing 0 * elev_mat

shared_arr = multiprocessing.Array(ctypes.c_short, elev.shape[0] * elev.shape[1])
shared_arr[:] = np.zeros(shape=elev.shape[0] * elev.shape[1]).astype('int16')


def get_rast(i):
    print(i)
    i = [i]
    shapes = [all_features[index]["geometry"] for index in i]
    range_rast = rasterio.features.rasterize(shapes,
                                        out_shape = elev.shape,
                                        transform = elev.transform,
                                        default_value= 1,
                                        all_touched=False)
    range_rast[np.logical_or(elev_mat < 0, elev_mat > 1000)] = 0
    return 1 # range_rast

def main():
    print("starting set up")
    pool = multiprocessing.Pool(2)
    out = list(pool.map(get_rast, list(range(20))))
    pool.close() 
    pool.join()

if __name__ == '__main__':
    main()

