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

rasterOptions(maxmemory = 1e+09) # want to handle all cells of 1-km global rasters in memory
rasterOptions(chunksize = 1e+09)

project_dir = "/home/chrisgraywolf/analysis_desktop/PA_matching/"
setwd(project_dir)

######################################## topography

gdalwarp(srcfile="data_input/elev.tif",
         dstfile="data_input/elev.vrt",
         t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
         tr=c(1000,1000),
         of="VRT")
gdal_translate("data_input/elev.vrt","data_processed/rasters/elev.tif",of="GTiff",co="COMPRESS=LZW")

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
gdal_translate("data_processed/non-threatened.vrt","data_processed/rasters/non-threatened.tif",of="GTiff",co="COMPRESS=LZW")

######################################## forest and travel time rasters

for(rast_name in c("cover_loss_v2","cover","gain")){
        gdalbuildvrt(gdalfile=paste0("data_input/",rast_name,c("-0000000000-0000000000.tif","-0000000000-0000023296.tif")),
               output.vrt=paste0("data_input/", rast_name, ".vrt"))

        gdalwarp(srcfile=paste0("data_input/", rast_name, ".vrt"),
                 dstfile=paste0("data_input/", rast_name, "_proj.vrt"),
                 tr=c(1000,1000),
                 of="VRT",
                 t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs",
                 te=extent(elev)[c(1,3,2,4)])

        gdal_translate(paste0("data_input/", rast_name, "_proj.vrt"),paste0("data_processed/rasters/", rast_name, ".tif"),of="GTiff",co=c("COMPRESS=LZW","BIGTIFF=YES"))
}

######################################## protected areas (dl from WDPA)

PAs = read_sf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp")
PAs$ID = as.numeric(factor(PAs$WDPAID))
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

######################################## country

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

######################################## quick maps

rast_names = c("elev","slope","threatened","non-threatened",
               "cover","cover_loss","gain","countries",
               "ecoregions","forest_biome","travel_time","PAs")

for(rast_name in rast_names){
        print(rast_name)
        rast = raster(paste0("data_processed/rasters/", rast_name, ".tif"))
        vel = velox(rast)
        vel$aggregate(factor=10, aggtype='median')
        rast = vel$as.RasterLayer()

        if(rast_name == ecoregions){
                rast[] = as.numeric(factor(round(rast[])))
        }

        png(paste0("output/input_rasters/", rast_name, ".png"), width=14, height=7, units="in", res=150)
        print(levelplot(rast, margin=FALSE, main=rast_name, maxpixels = 1e6))
        dev.off()
}








### NOT RUN









rast_names = c("elev","slope","threatened","non-threatened",
               "cover","cover_loss","gain","countries",
               "ecoregions","forest_biome","travel_time","PAs")

pdf("output/input_rasters.pdf", width=14, height=7)
for(rast_name in rast_names){
        print(rast_name)
        rast = raster(paste0("data_processed/rasters/", rast_name, ".tif"))
        vel = velox(rast)
        vel$aggregate(factor=10, aggtype='median')
        rast = vel$as.RayerLayer()

        
        print(levelplot(rast, margin=FALSE, main=rast_name, maxpixels = 1e6))
}
dev.off()





total = raster(paste0("data_processed/rasters/non-threatened.tif")) + raster(paste0("data_processed/rasters/threatened.tif"))
vel = velox(total)
vel$aggregate(factor=10, aggtype='median')
total = vel$as.RasterLayer()

total_pts = data.frame(rasterToPoints(total))

p = ggplot(total_pts, aes(x=x,y=y,fill=layer)) +
        geom_raster() +
        coord_fixed() +
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank()) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_gradientn(colors=magma(100))

pdf("output/forest_species.pdf", width=10, height=5)
print(p)
dev.off()










pdf("output/input_rasters.tif")
plot_list = lapply(rast_names, function(rast_name){
        rast = raster(paste0("data_processed/rasters/", rast_name, ".tif"))
        plot(rast, main=rast_name)
}
dev.off()





coverxloss -> annualzied loss
map grid

try matching








for(rast_name in rast_names){
        print(rast_name)
        rast = raster(paste0("data_processed/rasters/", rast_name, ".tif"))
        plot(rast)


}





