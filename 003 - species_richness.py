# https://stackoverflow.com/questions/34480328/numpy-summing-up-a-list-of-vectors
# https://stackoverflow.com/questions/53093073/why-is-numpy-sum-10-times-slower-than-the-operator
# https://stackoverflow.com/questions/34618301/python-multiprocessing-using-pickle-kwargs-and-function-references
# geopandas tries to load entire geometry into a col. I think... very slow and memory intensive
# https://stackoverflow.com/questions/10741346/numpy-most-efficient-frequency-counts-for-unique-values-in-an-array
# https://medium.com/distributed-computing-with-ray/ray-for-the-curious-fa0e019e17d3
# https://stackoverflow.com/questions/54582073/sharing-objects-across-workers-using-pyarrow

# import affine
# import ctypes
import fiona
# import functools
# import geojson
# import geopandas as gpd
# import glob
# import math
# import matplotlib.pyplot as plt
# import multiprocessing
import numpy as np
# import operator
# import os
import pandas as pd
# import psutil
# import random
import rasterio
# problem with ray https://github.com/ray-project/ray/issues/14535
import os
os.environ["RAY_OBJECT_STORE_ALLOW_SLOW_STORAGE"] = "1"
import ray
# import shapely
# from collections import defaultdict
from dbfread import DBF, FieldParser, InvalidValue
# from fiona import transform
# from osgeo import gdal, ogr
# from pyproj import Proj, transform
from rasterio import features
from rasterio.features import bounds as calculate_bounds

class MyFieldParser(FieldParser):
    def parse(self, field, data):
        try:
            return FieldParser.parse(self, field, data)
        except ValueError:
            return InvalidValue(data)

def get_indices(project_dir, keep_ids):
    species_data = pd.read_csv(project_dir + "data_processed/species_data.csv")
    species_data.loc[pd.isna(species_data['lower']), 'lower'] = -1e7
    species_data.loc[pd.isna(species_data['upper']), 'upper'] = 1e7
    ids = []
    for gp in ['MAMMALS','AMPHIBIANS','REPTILES']: # + ['BIRDS_' + str(x) for x in range(1,9)]:
        print(gp)
        shp_path = project_dir + "data_input/range_maps/" + gp + ".shp"
        dbf_path = project_dir + "data_input/range_maps/" + gp + ".dbf"
        gp_data = pd.DataFrame(iter(DBF(dbf_path, encoding='utf-8', char_decode_errors='ignore', parserclass=MyFieldParser)))
        if 'BIRDS_' in gp:
            gp_data['presence'] = gp_data['presenc']
        gp_ids = gp_data.loc[(gp_data['presence'] <= 2) & (gp_data['origin'] <= 2)]['id_no'].unique().tolist()
        gp_ids = list(set(keep_ids) & set(gp_ids))
        for gp_id in gp_ids:
            ids += [{'path':shp_path,
                'id':int(gp_id),
                'rows':gp_data.index[gp_data['id_no'] == gp_id].tolist(),
                'lower':float(species_data.loc[species_data['taxonid'] == gp_id]['lower']),
                'upper':float(species_data.loc[species_data['taxonid'] == gp_id]['upper'])}]
    return ids

@ray.remote
class MyActor(object):
    def __init__(self, elev_mat_id, project_dir):
        self.elev = rasterio.open(project_dir + 'data_input/elev.tif')
        self.elev_mat = elev_mat_id
        self.richness = (0 * self.elev_mat).astype('int16')
    def add_range(self, id_info):
        print(id_info['id'])
        all_features = fiona.open(id_info['path'])
        shapes = [all_features[index]["geometry"] for index in id_info['rows']]
        range_rast = rasterio.features.rasterize(shapes,
                                        out_shape = self.elev.shape,
                                        transform = self.elev.transform,
                                        default_value = 1,
                                        all_touched = False)
        range_rast[np.logical_or(self.elev_mat < id_info['lower'], self.elev_mat > id_info['upper'])] = 0
        self.richness += range_rast
    def get_richness(self):
        return self.richness

def calculate_richness(project_dir, keep_ids, rast_name, num_cpu=12):
    id_info_list = get_indices(project_dir, keep_ids)
    elev = rasterio.open(project_dir + 'data_input/elev.tif')
    elev_mat = elev.read(1)
    #
    ray.init(object_store_memory=45*1024*1024*1024)
    elev_mat_id = ray.put(elev_mat)
    actors = [MyActor.remote(elev_mat_id, project_dir) for _ in range(num_cpu)]
    pool = ray.util.ActorPool(actors)
    list(pool.map_unordered(lambda actor, v: actor.add_range.remote(v), id_info_list))
    #
    results = ray.get([actor.get_richness.remote() for actor in actors])
    ray.shutdown()
    richness = results[0].copy() # EDIT: adding .copy() as needed with up to date numpy
    for i in range(1,num_cpu):
        richness += results[i]
    #
    rast = rasterio.open(
        project_dir + 'data_processed/' + rast_name + '.tif',
        'w',
        driver='GTiff',
        height=elev.shape[0],
        width=elev.shape[1],
        count=1,
        dtype=rasterio.int16,
        crs='EPSG:4326',
        compress='lzw',
        transform=elev.transform,
    )
    rast.write(richness, 1)
    rast.close()

def main():
    project_dir = "/home/onyxia/work/PA_matching/"
    species_data = pd.read_csv(project_dir + "data_processed/species_data.csv")
    #
    for gp in list(species_data['group'].unique()):
        calculate_richness(project_dir, species_data[species_data['group'] == gp]['taxonid'].tolist(), gp)

if __name__== "__main__":
  main()

