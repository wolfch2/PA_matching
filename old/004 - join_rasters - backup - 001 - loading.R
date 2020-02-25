require(raster)
require(sf)
require(fasterize)
require(rgdal)
require(gdalUtils)

rasterOptions(maxmemory = 1e+09) # want to handle all cells of 1-km global rasters in memory
rasterOptions(chunksize = 1e+09)

project_dir = "/home/chrisgraywolf/analysis_desktop/PA_matching/"
setwd(project_dir)

######################################## land-related

elev = raster("data_input/elev.tif")
mask = raster("data_input/mask.tif")
land = raster("data_input/land.tif")

######################################## species richness

threatened = readAll(raster("data_processed/threatened.tif"))
non_threatened = readAll(raster("data_processed/non-threatened.tif"))

######################################## forest rasters

loss = mosaic(readAll(raster("data_input/loss-0000000000-0000000000.tif")),
              readAll(raster("data_input/loss-0000000000-0000023296.tif")), fun=mean)

cover = mosaic(readAll(raster("data_input/cover-0000000000-0000000000.tif")),
              readAll(raster("data_input/cover-0000000000-0000023296.tif")), fun=mean)      

cover_loss = mosaic(readAll(raster("data_input/cover_loss_v2-0000000000-0000000000.tif")),
              readAll(raster("data_input/cover_loss_v2-0000000000-0000023296.tif")), fun=mean)      

######################################## protected areas (dl from WDPA)

PAs = read_sf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp")
PAs_rast = fasterize(PAs, raster=elev, field="WDPAID")

######################################## ecoregions https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world

eco_url = "https://c402277.ssl.cf1.rackcdn.com/publications/15/files/original/official_teow.zip?1349272619"
system(paste0("wget --content-disposition \"", eco_url, "\" -nc --tries=0 -P ", project_dir, "data_input"))
system("unzip data_input/official_teow.zip?1349272619 -d data_input")
file.remove("data_input/official_teow.zip?1349272619")

eco_regions = read_sf("data_input/official/wwf_terr_ecos.shp")
eco_rast = fasterize(eco_regions, elev, "ECO_ID")


