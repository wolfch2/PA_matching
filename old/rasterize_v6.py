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

class ActorPool(object):
    def __init__(self, actors):
        self.actors = actors
    #
    def map(self, fn, values):
        values = list(values)
        idle = list(self.actors)
        running = {}
        results = {}
        for i, v in enumerate(values):
            if not idle:
                [r], _ = ray.wait(list(running), num_returns=1)
                j, a = running.pop(r)
                results[j] = ray.get(r)
                idle.append(a)
            a = idle.pop()
            running[fn(a, v)] = (i, a)
            i += 1
        while running:
            [r], _ = ray.wait(list(running), num_returns=1)
            j, _ = running.pop(r)
            results[j] = ray.get(r)
        return [results[i] for i in range(len(results))]


@ray.remote
class MyActor(object):
    def __init__(self):
        self.value = 0
    def f(self, v):
        self.value = self.value + v
    def get_value(self):
        return self.value

ray.init()
actors = [MyActor.remote() for _ in range(4)]
pool = ActorPool(actors)
ffs = pool.map(lambda actor, v: actor.f.remote(v), list(range(10,20)))

results = ray.get([actor.get_value.remote() for actor in actors])



