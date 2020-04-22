require(raster)
require(sf)
require(fasterize)
require(rgdal)
require(gdalUtils)
require(ggplot2)
require(RColorBrewer)
require(viridis)
require(rasterVis)
require(landscapetools)
require(velox)
require(rnaturalearth)
require(ncdf4)

rasterOptions(maxmemory = 1e+09) # want to handle all cells of 1-km global rasters in memory
rasterOptions(chunksize = 1e+09)

project_dir = "/home/chrisgraywolf/shared/analysis/PA_matching/"
setwd(project_dir)

######################################## topography (slope)

elev = raster("data_processed/rasters/elev.tif")
slope = terrain(elev)

writeRaster(slope, "data_processed/rasters/slope.tif")

######################################## species richness

gdalwarp(srcfile="data_processed/threatened.tif",
         dstfile="data_processed/threatened.vrt",
         t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
         tr=c(1000,1000),
         of="VRT")
gdal_translate("data_processed/threatened.vrt","data_processed/rasters/threatened.tif",of="GTiff",co="COMPRESS=LZW")

gdalwarp(srcfile="data_processed/non-threatened.tif",
         dstfile="data_processed/non-threatened.vrt",
         t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
         tr=c(1000,1000),         
         of="VRT")
gdal_translate("data_processed/non-threatened.vrt","data_processed/rasters/non_threatened.tif",of="GTiff",co="COMPRESS=LZW")

######################################## Curtis et al. data (drivers)

drivers_url = "https://science.sciencemag.org/highwire/filestream/715492/field_highwire_adjunct_files/2/aau3445-Data-S3.tif"
system(paste0("wget --content-disposition \"", drivers_url, "\" -nc --tries=0 -P ", project_dir, "data_input"))

drivers = raster("data_input/aau3445-Data-S3.tif")
crs(drivers) = CRS("+proj=igh +towgs84=0,0,0") # pretty sure this should be interupted form (o.w. reprojected version is wrong)
writeRaster(drivers, "data_input/drivers.tif")

gdalwarp(srcfile="data_input/drivers.tif",
         dstfile="data_processed/drivers_proj.vrt",
         t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
         tr=c(1000,1000),
         te=extent(elev)[c(1,3,2,4)],
         of="VRT")

gdal_translate("data_processed/drivers_proj.vrt",
               "data_processed/rasters/drivers.tif",of="GTiff",co="COMPRESS=LZW")

######################################## lossyear

gdalwarp(srcfile="data_input/lossyear.tif",
                 "data_input/lossyear_proj.vrt",
                 tr=c(1000,1000),
                 of="VRT",
                 t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
                 te=extent(elev)[c(1,3,2,4)])

gdal_translate("data_input/lossyear_proj.vrt",
               "data_processed/rasters/lossyear.tif",
               of="GTiff",
               co=c("COMPRESS=LZW","BIGTIFF=YES"))

######################################## forest, travel time, HPD, etc.

for(rast_name in c("cover_loss","loss","cover","gain","pop_dens","travel_time","lossyear")){
        gdalbuildvrt(gdalfile=list.files("data_input",pattern=paste0("^",rast_name,"-.*tif$"),full.names=TRUE),
               output.vrt=paste0("data_input/", rast_name, ".vrt"))

        gdalwarp(srcfile=paste0("data_input/", rast_name, ".vrt"),
                 dstfile=paste0("data_input/", rast_name, "_proj.vrt"),
                 tr=c(1000,1000),
                 of="VRT",
                 t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
                 te=extent(elev)[c(1,3,2,4)])

        gdal_translate(paste0("data_input/", rast_name, "_proj.vrt"),paste0("data_processed/rasters/", rast_name, ".tif"),
                       of="GTiff",
                       co=c("COMPRESS=LZW","BIGTIFF=YES"))
}

######################################## forest carbon

carbon_URL = "https://daac.ornl.gov/bundle/Global_Biomass_1950-2010_1296.zip"
system(paste0("wget --content-disposition \"", carbon_URL, "\" -nc --tries=0 -P ", project_dir, "data_input"))
system("unzip data_input/Global_Biomass_1950-2010_1296.zip -d data_input")
file.remove("data_input/Global_Biomass_1950-2010_1296.zip")

#  Aboveground forest tree biomass (tonnes carbon / 1-degree grid cell)
carbon = stack("data_input/Global_Biomass_1950-2010_1296/data/historical_global_1-degree_forest_biomass.nc4", varname="AboveGroundBiomass")
carbon = carbon[["X50.5"]] # 50.5 years after 1/1/1950 -> 2000
carbon = carbon / area(carbon) # divide by km_2 per 1-degree pixel area
writeRaster(carbon, "data_input/carbon.tif")

gdalwarp(srcfile="data_input/carbon.tif",
                 "data_input/carbon_proj.vrt",
                 tr=c(1000,1000),
                 of="VRT",
                 t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
                 te=extent(elev)[c(1,3,2,4)])

gdal_translate("data_input/carbon_proj.vrt",
               "data_processed/rasters/carbon.tif",
               of="GTiff",
               co=c("COMPRESS=LZW","BIGTIFF=YES"))

######################################## protected areas (dl from WDPA)

PAs = read_sf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp")
PAs$ID = as.numeric(factor(PAs$WDPAID)) # cannot use WDPAID directly!  see below..
PAs_proj = st_transform(PAs, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
PAs_rast = fasterize(PAs_proj, raster=elev, field="ID") # could prioritize by IUCN_CAT or...
writeRaster(PAs_rast, "data_processed/rasters/PAs.tif") # WDPAID seems to default to int, but values exceed max allowed (also lose too much precision w/ float..)
PAs_tab = PAs
PAs_tab$geometry = NULL
write.csv(data.frame(PAs_tab), "data_processed/PAs_tab.csv", row.names=FALSE)

PAs_binary = fasterize(PAs_proj, raster=elev) # could prioritize by IUCN_CAT or...
PAs_binary[is.na(PAs_binary[])] = 0
writeRaster(PAs_binary, "data_processed/PAs_binary.tif")
system("gdal_proximity.py data_processed/PAs_binary.tif data_processed/PAs_dist.tif -distunits GEO")
PAs_dist = raster("data_processed/PAs_dist.tif")
PAs_buffer = PA_dists <= 10e3
writeRaster(PAs_buffer, "data_processed/rasters/PAs_buffer.tif")
file.remove("data_processed/PAs_dist.tif")

######################################## ecoregions and forest biomes https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world

eco_url = "https://c402277.ssl.cf1.rackcdn.com/publications/15/files/original/official_teow.zip?1349272619"
system(paste0("wget --content-disposition \"", eco_url, "\" -nc --tries=0 -P ", project_dir, "data_input"))
system("unzip data_input/official_teow.zip?1349272619 -d data_input")
file.remove("data_input/official_teow.zip?1349272619")

eco_regions = read_sf("data_input/official/wwf_terr_ecos.shp")
eco_regions_proj = st_transform(eco_regions, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")

eco_rast = fasterize(eco_regions_proj, elev, "ECO_ID")
writeRaster(eco_rast, "data_processed/rasters/ecoregions.tif")

forest_biome = fasterize(eco_regions_proj, elev, "BIOME")
forest_biome[] = forest_biome[] %in% c(1:6, 14)
writeRaster(forest_biome, "data_processed/rasters/forest_biome.tif")

######################################## countries

worldmap = ne_download(scale = 10,
                       type = "countries",
                       category = "cultural",
                       destdir = tempdir(),
                       load = TRUE,
                       returnclass = "sf")
worldmap$ID = 1:nrow(worldmap)
worldmap_proj = st_transform(worldmap, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
countries_rast = fasterize(worldmap_proj, elev, "ID")
writeRaster(countries_rast, "data_processed/rasters/countries.tif")

