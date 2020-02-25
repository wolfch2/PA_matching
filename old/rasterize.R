require(sf)
require(fasterize)
require(raster)
require(rgdal)
require(foreach)
require(doMC)

setwd("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/range_maps/")


r = raster(resolution=c(0.01,0.01))

ogrListLayers("BOTW.gdb")

test = read_sf("BOTW.gdb", layer="All_Species")

start = Sys.time()
out = fasterize(test[1,],r) # 1.58 sec
Sys.time() - start



start = Sys.time()
        writeRaster(out,paste0("/home/chrisgraywolf/analysis_desktop/PA_matching/temp/","1",".tif"))
Sys.time() - start
        
        
        
rast_list = foreach(i=1:2) %do% {
        print(i)
        out = fasterize(test[i,],r)

}

for(i in 1:2){
        print(i)
        out = fasterize(test[i,],r)
        writeRaster(out,paste0("/home/chrisgraywolf/analysis_desktop/PA_matching/temp/",i,".tif"))
}







setwd("/home/chrisgraywolf/analysis_desktop/PA_matching/data_input/PAs/")








PAs = read_sf("WDPA_Jan2020-shapefile-polygons.shp")


start = Sys.time()
test = fasterize(PAs, r)
Sys.time() - start # 4.7 sec :O






mammals = read_sf("MAMMALS.shp")

ogrListLayers("BOTW.gdb")

test = read_sf("BOTW.gdb", layer="All_Species")



rast_list = foreach(i=1:10) %do% {
        print(i)
        out = fasterize(test[i,],r)

}



