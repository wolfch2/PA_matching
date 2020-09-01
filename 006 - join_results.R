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
require(foreign)
require(scales)
require(cowplot)
require(rnaturalearth)
require(countrycode)
require(readr)
require(pbapply)
require(stringr)
require(wbstats)
require(reshape2)
require(spNNGP)
require(foreach)
require(MCMC.OTU)
require(spaMM)
require(GpGp)

# https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

quants = function(x) data.frame(ymin=quantile(x,0.25),ymax=quantile(x,0.75))

project_dir = "/home/chrisgraywolf/shared/analysis/PA_matching/"
setwd(project_dir)

######################################## join info. from PA table

PA_df = read.csv("data_processed/data_matched.csv", as.is=TRUE, na.strings=c("","NaN","NA"))
PA_table = read.csv("data_processed/PAs_tab.csv", as.is=TRUE)
PA_table = PA_table[PA_table$ID %in% PA_df$PAs,]

PA_df$GIS_AREA = tapply(PA_table$GIS_AREA, PA_table$ID, sum)[as.character(PA_df$PAs)]
PA_df$IUCN_CAT = tapply(PA_table$IUCN_CAT, PA_table$ID, Mode)[as.character(PA_df$PAs)]
PA_df$GOV_TYPE = tapply(PA_table$GOV_TYPE, PA_table$ID, Mode)[as.character(PA_df$PAs)]
PA_df$ISO3 = tapply(PA_table$PARENT_ISO, PA_table$ID, Mode)[as.character(PA_df$PAs)]
PA_df$NAME = tapply(PA_table$NAME, PA_table$ID, Mode)[as.character(PA_df$PAs)]
PA_df$WDPAID = PA_table$WDPAID[match(PA_df$PAs, PA_table$ID)]

PAs = read_sf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp")
PAs = PAs[match(PA_df$WDPAID, PAs$WDPAID),]
PAs_proj = st_transform(PAs, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
centroids = st_centroid(PAs_proj)
coords_Behr = st_coordinates(centroids)
coords_latlong = st_coordinates(st_transform(centroids,4326))

PA_df$X = coords_Behr[,"X"]
PA_df$Y = coords_Behr[,"Y"]
PA_df$long = coords_latlong[,"X"]
PA_df$lat = coords_latlong[,"Y"]

######################################## add country, GDP, simplified category, etc.

PA_df$ISO3_single = pbsapply(PA_df$ISO3, function(x){
	if(nchar(x) == 3) return(x)
	out = names(sort(table(str_split(x,";")[[1]]),decreasing=TRUE)[1])
	return(out)
})
PA_df$ISO3_single[PA_df$ISO3_single == "TWN"] = "CHN"
PA_df$Continent = countrycode(PA_df$ISO3_single, "iso3c", "continent")
PA_df$Continent[countrycode(PA_df$ISO3_single, "iso3c", "region") %in% c("South America")] = "South America"
PA_df$Continent[PA_df$Continent == "Americas"] = "North America"

PA_df$Country = countrycode(PA_df$ISO3_single, "iso3c", "country.name")

# add developed/developing info from https://unstats.un.org/unsd/methodology/m49/
dev_list = read_csv("/home/chrisgraywolf/shared/analysis/PA_loss/data/developing.csv")
PA_df$developed = PA_df$ISO3_single %in% dev_list$`ISO-alpha3 code`

### ADD GDP
GDP_raw = wb(indicator="NY.GDP.PCAP.KD", startdate=2000, enddate=2018) # GDP per capita, PPP (current international $)
GDP = aggregate(value ~ iso3c, GDP_raw, function(x) mean(x,na.rm=TRUE))
PA_df$GDP = GDP[match(PA_df$ISO3_single, GDP$iso3c),"value"]

### ADD SIMPLIFIED CAT; update: also add cat_pure and indigenous
PA_df$cat_simple = 'Unknown'
PA_df$cat_simple[PA_df$IUCN_CAT %in% c('Ia','Ib','II','III','IV')] = 'Strict'
PA_df$cat_simple[PA_df$IUCN_CAT %in% c('V','VI')] = 'Nonstrict'

######################################## export

PA_df$group = "omit"
PA_df$group[PA_df$STATUS_YR <= 2000] = "main"
PA_df$group[PA_df$STATUS_YR %in% 2002:2017] = "change"

PA_df$matched = ! is.na(PA_df$n_control)

#test = PA_df[PA_df$group == "main",]
#sum(is.na(test$Control_loss)) # hmm.. something seems wrong -> look into PA 12621
#test = PA_df[PA_df$group == "main" & ! is.na(PA_df$n_control),]
#sum(is.na(test$Control_loss)) # nvm.. we're okay

saveRDS(PA_df, "data_processed/PA_df.RDS")
