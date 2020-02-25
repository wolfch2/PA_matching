import ee
import rasterio
import pandas as pd
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive

ee.Initialize()

################################################## export rasters for analysis

scale = 0.008983153 # checked w/ small region export at 1000 m

GEE_rasts = {
    #'land' = ee.Image("MODIS/MOD44W/MOD44W_005_2000_02_24").select('water_mask').eq(0).rename(['land']).unmask(0),
    #'gtopo30' = ee.Image('USGS/GTOPO30'),
    'loss':ee.Image("UMD/hansen/global_forest_change_2018_v1_6").select('loss').reduceResolution(reducer=ee.Reducer.mean(), bestEffort=True),
    'lossyear':ee.Image("UMD/hansen/global_forest_change_2018_v1_6").select('lossyear').reduceResolution(reducer=ee.Reducer.mode(), bestEffort=True),
    'gain':ee.Image("UMD/hansen/global_forest_change_2018_v1_6").select('gain').reduceResolution(reducer=ee.Reducer.mean(), bestEffort=True),
    'cover':ee.Image("UMD/hansen/global_forest_change_2018_v1_6").select('treecover2000').reduceResolution(reducer=ee.Reducer.mean(), bestEffort=True),
    'mask':ee.Image("UMD/hansen/global_forest_change_2018_v1_6").select('datamask').reduceResolution(reducer=ee.Reducer.mode(), bestEffort=True),
    'travel_time':ee.Image("Oxford/MAP/accessibility_to_cities_2015_v1_0")
}

for rast_name in GEE_rasts.keys():
    ee.batch.Export.image.toDrive(GEE_rasts[rast_name],
        description=rast_name,
        crsTransform=[scale, 0, -180, 0, -scale, 90],
        crs='EPSG:4326',
        dimensions=str(int(360/scale)) + 'x' + str(int(180/scale)),
        folder="GEE_rasts",maxPixels=1e11).start()

cover_loss = ee.Image("UMD/hansen/global_forest_change_2018_v1_6").select('treecover2000'). \
    multiply(ee.Image("UMD/hansen/global_forest_change_2018_v1_6").select('loss')).reduceResolution(reducer=ee.Reducer.mean(), bestEffort=True)

ee.batch.Export.image.toDrive(cover_loss,
        description='cover_loss_v2',
        crsTransform=[scale, 0, -180, 0, -scale, 90],
        crs='EPSG:4326',
        dimensions=str(int(360/scale)) + 'x' + str(int(180/scale)),
        folder="GEE_rasts",maxPixels=1e11).start()


while True:
    states = pd.DataFrame(map(lambda x: x.status(), ee.batch.Task.list()))['state']
    if 'RUNNING' in states.unique() or 'READY' in states.unique():
        print 'Waiting on GEE.  Tasks: ' + str(((states == 'READY') | (states == 'RUNNING')).sum())
    time.sleep(60)
    else:
        break

################################################## download from Google Drive

gauth = GoogleAuth()
gauth.LoadCredentialsFile("/home/chrisgraywolf/analysis_desktop/mycreds.txt")
if gauth.access_token_expired:
    gauth.Refresh()

drive = GoogleDrive(gauth)

root_files = drive.ListFile({'q': "'root' in parents and trashed=false"}).GetList()
root_files_df = pd.DataFrame(root_files)
folder_id = root_files_df['id'][root_files_df['title'] == "GEE_rasts"].tolist()[0]
files = drive.ListFile({'q': "'{}' in parents and trashed=false".format(folder_id)}).GetList()
for file in files:
    file.GetContentFile("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/" + file['title'])
    # file.Delete()

