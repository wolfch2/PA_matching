require(RCurl)
require(rjson)
require(plyr)
require(pbapply)
require(stringr)
require(rgdal)
require(sf)
require(foreach)
require(doMC)
require(raster)
require(sf)
require(fasterize)
require(rgdal)
require(gdalUtils)

setwd("/home/chrisgraywolf/analysis_desktop/PA_matching/")

token = "1f6e9566aa509bcf87bdecb32901e1a7f53d04785749d6afe3db9b3b5b043002"

######################################################## get species

species_data = NULL
page = 0
while(1){
    	print(page)
	page_data = getURL(paste0("http://apiv3.iucnredlist.org/api/v3/species/page/",page,"?token=",token))
	page_data = fromJSON(page_data)
	if(page_data$count == 0)
	    	break
	page_data = do.call("rbind", page_data[["result"]])
	species_data = rbind(species_data, page_data)
	page = page+1
}
species_data = data.frame(species_data)

for(col_name in colnames(species_data)){ # have many nulls in here, which makes this get stored as a list
	species_data[,col_name] = sapply(species_data[,col_name],function(x) ifelse(is.null(x),"",x))
}

species_data = species_data[species_data$class_name %in% c("MAMMALIA","AMPHIBIA","REPTILIA","AVES"),]
species_data = species_data[species_data$infra_rank == "" & species_data$population == "",]

# saveRDS(species_data, "species_data_v1.RDS")

######################################################## habitats

urls = paste0("http://apiv3.iucnredlist.org/api/v3/habitats/species/id/",species_data$taxonid,"?token=", token)
commands = paste0("wget \"", urls, "\" -nc --tries=0 -P /home/chrisgraywolf/analysis_desktop/PA_matching/temp")
write(commands, file = "commands.txt")
system("parallel --jobs 4 < commands.txt")

files = paste0("/home/chrisgraywolf/analysis_desktop/PA_matching/temp/", species_data$taxonid, "?token=", token)

species_data$forest_exclusive = pbsapply(1:nrow(species_data), function(i){
	print(i)
	syn_file = files[i]
	data = fromJSON(file=syn_file)
        codes = unique(sapply(str_split(sapply(data$result,function(x) x$code),"\\."),function(x) x[1]))
        ("1" %in% codes & length(codes) == 1)
})

species_data = species_data[species_data$forest_exclusive,]

file.remove(files)

# saveRDS(species_data, "species_data_v2.RDS")

######################################################## altitude limits

urls = paste0("http://apiv3.iucnredlist.org/api/v3/species/id/",species_data$taxonid,"?token=", token)
commands = paste0("wget \"", urls, "\" -nc --tries=0 -P /home/chrisgraywolf/analysis_desktop/PA_matching/temp")
write(commands, file = "commands.txt")
system("parallel --jobs 4 < commands.txt")

files = paste0("/home/chrisgraywolf/analysis_desktop/PA_matching/temp/", species_data$taxonid, "?token=", token)

species_data$species_info_raw = pblapply(1:nrow(species_data), function(i){
	print(i)
	syn_file = files[i]
	data = fromJSON(file=syn_file)
       	return(data)
})

species_data$lower = pbsapply(species_data$species_info_raw, function(x){
	lower = x$result[[1]]$elevation_lower
	if(is.null(lower)) lower = NA
	return(lower)
})

species_data$upper = pbsapply(species_data$species_info_raw, function(x){
	upper = x$result[[1]]$elevation_upper
	if(is.null(upper)) upper = NA
	return(upper)
})

# saveRDS(species_data, "species_data_v3.RDS")

species_data$species_info_raw = NULL
species_data$forest_exclusive = NULL

file.remove(files)

######################################################## countries

urls = paste0("http://apiv3.iucnredlist.org/api/v3/species/countries/id/",species_data$taxonid,"?token=", token)
commands = paste0("wget \"", urls, "\" -nc --tries=0 -P /home/chrisgraywolf/analysis_desktop/PA_matching/temp")
write(commands, file = "commands.txt")
system("parallel --jobs 4 < commands.txt")

files = paste0("/home/chrisgraywolf/analysis_desktop/PA_matching/temp/", species_data$taxonid, "?token=", token)

species_data$countries = pbsapply(1:nrow(species_data), function(i){
	print(i)
	syn_file = files[i]
	data = fromJSON(file=syn_file)        
	country_tab = data.frame(do.call("rbind",data$result))
      	out = paste(country_tab$code[country_tab$presence == "Extant" & country_tab$origin == "Native"],collapse=",")
	return(out)
})

# saveRDS(species_data, "species_data_v4.RDS")

file.remove(files)

species_data = species_data[! species_data$category %in% c("EW","EX"),]
species_data$group = factor(species_data$category %in% c("VU","EN","CR"), c(FALSE,TRUE), c("non-threatened","threatened"))

write.csv(species_data, "data_processed/species_data.csv", row.names=FALSE)

######################################################## convert bird ranges to shapefiles

ogrListLayers("data_input/range_maps/BOTW.gdb")
bird_ranges = read_sf("data_input/range_maps/BOTW.gdb", layer="All_Species")

names(bird_ranges)[names(bird_ranges) == "SISID"] <- "id_no"
names(bird_ranges)[names(bird_ranges) == "PRESENCE"] = "presence"
names(bird_ranges)[names(bird_ranges) == "ORIGIN"] = "origin"

ids = unique(bird_ranges$id_no)
id_gps = split(ids, cut(ids, quantile(ids,seq(0,1,length=9))))

registerDoMC(8)
foreach(i=1:8, .packages="sf") %dopar%{
        select = bird_ranges[bird_ranges$id_no %in% id_gps[[i]],]
        write_sf(select, paste0("data_input/range_maps/BIRDS_", i, ".shp"))
}

######################################################## build elevation raster (for masking ranges etc.)

gdalbuildvrt(gdalfile=list.files("data_input",pattern="elev.*tif$",full.names=TRUE),
        output.vrt="data_input/elev.vrt")

gdalwarp(srcfile=paste0("data_input/elev.vrt"),
         dstfile=paste0("data_input/elev_proj.vrt"),
         tr=c(1000,1000),
         of="VRT",
         t_srs="+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")

gdal_translate("data_input/elev_proj.vrt",
               "data_processed/rasters/elev.tif",
               of="GTiff",
               co=c("COMPRESS=LZW","BIGTIFF=YES"))

