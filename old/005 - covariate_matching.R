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

rasterOptions(maxmemory = 1e+09) # want to handle all cells of 1-km global rasters in memory
rasterOptions(chunksize = 1e+09)

project_dir = "/home/chrisgraywolf/analysis_desktop/PA_matching/"
setwd(project_dir)

########################################

cover = raster("data_processed/rasters/cover.tif")
keep = which(cover[] >= 30)

cover[is.na(cover[])] = 0
keep = cover >= 30 # only retain these cells

rast_names = c("elev","slope","threatened","non-threatened",
               "cover","cover_loss","gain",
               "ecoregions","forest_biome","travel_time","PAs")

rast_data = matrix(NA, sum(keep[]), length(rast_names))
colnames(rast_data) = rast_names

df = NULL
for(rast_name in rast_names){
        print(rast_name)
        rast = raster(paste0("data_processed/rasters/", rast_name, ".tif"))
        rast_data[,rast] = rast[keep]


        vel = velox(rast)
        vel$aggregate(factor=10, aggtype='median')
        rast = vel$as.RayerLayer()










