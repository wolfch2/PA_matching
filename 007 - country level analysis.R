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
require(spgwr)
require(parallel)
require(tidyverse)
require(doMC)
require(plyr)
require(stars)
require(lme4)
require(spmoran)
require(grid)
require(ggforce)
require(RCurl)
require(jsonlite)
require(data.table)
require(dtplyr)
require(unix)
require(smoothr)
require(ggrepel)

rlimit_as(58e9) # soft (adjustable limit)

project_dir = "/home/chrisgraywolf/shared/analysis/PA_matching/"
setwd(project_dir)

######################################## set up country-level data

countries_sf = ne_load(scale = 10, # use ne_download first time
                       type = "countries",
                       category = "cultural",
                       destdir = tempdir(),
                       returnclass = "sf")
countries_sf$ID = 1:nrow(countries_sf)

total_species = read.csv("data_processed/species_data.csv") %>% # LC-CR forest exclusive
        filter(countries != "") %>%
        mutate(countries = as.character(countries)) %>%
        pull(countries) %>%   
        str_split(",") %>%
        unlist %>%
        table

countries_sf$total = as.numeric(total_species[countries_sf$ISO_A2])
countries_sf$total[is.na(countries_sf$total)] = 0

### add info from rasters

countries_rast = readAll(raster("data_processed/rasters/countries.tif"))
cover = readAll(raster("data_processed/rasters/cover.tif"))
carbon = readAll(raster("data_processed/rasters/carbon.tif"))
cover_loss = readAll(raster("data_processed/rasters/cover_loss.tif"))

PAs = readAll(raster("data_processed/rasters/PAs.tif"))
PAs[is.na(PAs[])] = 0
PA_df = read.dbf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.dbf", as.is=TRUE)
PA_df$ID = as.numeric(factor(PA_df$WDPAID))
PAs[PAs[] %in% PA_df$ID[PA_df$STATUS == "Proposed"]] = 0
PAs[PAs[] > 0] = 1

dt = data.table(cover=cover[],
                PAs=PAs[],
                countries=countries_rast[],
                carbon=carbon[],
                cover_loss=cover_loss[])
cover = PAs = countries = carbon = cover_loss = NULL; gc();
dt = dt[! is.na(countries),]; gc();

prot_df = data.frame(dt[, .(protected_raw=mean(PAs),
                            prop_forest=mean(cover >= 30),
                            area_forest=sum(cover >= 30),
                            cover=sum(cover),
                            cover_loss=sum(cover_loss, na.rm=TRUE),
                            carbon=sum(carbon, na.rm=TRUE)), by=c('countries'),])
prot_df$loss = 1 - ((prot_df$cover - prot_df$cover_loss)/prot_df$cover)^(1/18)

countries_sf = merge(countries_sf, prot_df, by.x="ID", by.y="countries", all.x=TRUE)

###

eff_df = readRDS("data_processed/PA_df.RDS") %>%
        filter(matched) %>%
        filter(group == "main") %>%
        group_by(ISO3_single) %>%
        dplyr::summarise(effectiveness = mean(Control_loss)/mean(PA_loss),
                         n = length(PA_loss)) %>%
        data.frame

setdiff(eff_df$ISO3_single, countries_sf$ISO_A3)
countries_sf$ISO_A3[countries_sf$NAME == "Norway"] = "NOR"
countries_sf$ISO_A3[countries_sf$NAME == "France"] = "FRA"

countries_sf = merge(countries_sf, eff_df, by.x="ISO_A3", by.y="ISO3_single", all.x=TRUE)

countries_sf$protected = countries_sf$protected_raw * countries_sf$effectiveness
countries_sf$adj_threat = countries_sf$total / countries_sf$protected
countries_sf$log_carbon = log(countries_sf$carbon, base=10)
countries_sf$adj_carbon = countries_sf$log_carbon / countries_sf$protected
countries_sf$adj_loss = countries_sf$loss / countries_sf$protected

saveRDS(countries_sf, "temp/countries_sf.RDS")
countries_sf = readRDS("temp/countries_sf.RDS")

country_data = countries_sf %>%
        filter(n >= 15 & total >= 5 & area_forest >= 100^2)

cor(country_data$protected, country_data$protected_raw) # 0.54 (paper text)
cor(country_data$total, country_data$log_carbon) # 0.37 (paper text)

##################### mapping etc.

robin = CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
bbox = st_sf(geometry=st_as_sfc(st_bbox())) %>% densify(100)

bbox_robin = st_transform(bbox, robin)
countries_sf_robin = st_transform(countries_sf, robin)
country_data_robin = st_transform(country_data, robin)

ex = st_bbox(country_data_robin)

##################### total (adjusted)

col_mult_adj = max(country_data$total)/max(country_data$protected) # for color scale
brks = c(100,1000,10000)

bg_adj = expand.grid(protected=seq(0,
	                           1.05*max(country_data$protected,na.rm=TRUE),
                                   length=1e2),
	             total=seq(1e-10,
                               1.05*max(country_data$total,na.rm=TRUE),length=1e2))
bg_adj$adj_threat = bg_adj$total / bg_adj$protected
col_lim = atan(range(bg_adj$adj_threat)/col_mult_adj)

p_bot = ggplot() +
	geom_sf(data=bbox_robin, fill="#ccccff") +
	geom_sf(data=countries_sf, fill="#888888", colour=NA) +	
	geom_sf(data=country_data, aes(fill=atan(adj_threat/col_mult_adj)), color=NA) +
	geom_sf(data=countries_sf, fill=NA, colour="black", size=0.05) +	
	theme_bw() +
	coord_sf(xlim=ex[c(1,3)], ylim=ex[c(2,4)]) +
	scale_fill_gradientn(colours=plasma(50)[15:50],
                             breaks=atan(brks/col_mult_adj),
                             labels=brks,
                             limits=col_lim) +
	theme(axis.line=element_blank(),
              plot.margin=margin(3,3,3,20),              
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              axis.title=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              legend.position="bottom",
              legend.key.width=unit(2.5,"lines"),
              legend.margin=margin(3,11,3,3),
              legend.background = element_rect(fill=NULL,color="black")) +
	guides(fill=guide_colorbar(title="Adjusted threat  \nindex", nbin=300))

p_top = ggplot(country_data, aes(y=total,x=protected)) +
	geom_raster(data=bg_adj, aes(fill=atan(adj_threat/col_mult_adj))) +
	geom_point() +
	scale_fill_gradientn(colours=plasma(50)[15:50],
                             breaks=atan(brks/col_mult_adj),
                             labels=brks,
                             limits=col_lim) +
	theme_bw() +
	xlab("(Proportion of forested area protected) * (PA effectiveness)") +
	theme(axis.text=element_text(color="black"),
              plot.margin=margin(3,3,3,20),
	      legend.position="bottom") +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	ylab("Number of forest vertebrates") +
	guides(fill=FALSE) +
        geom_text_repel(data = country_data[country_data$total >= 300 | country_data$protected > 1.5,],
                        aes(label=NAME),
                        size=3)

all = plot_grid(p_top, p_bot, ncol=1, labels=c("A.","B."), rel_heights=c(2,1.5), hjust=0)

png("output/adj_threat.png", width=6, height=9, units="in", res=400)
all
dev.off()

pdf("output/adj_threat.pdf", width=6, height=9)
all
dev.off()

##################### carbon and deforestation plots
# https://stackoverflow.com/questions/46058055/r-ggplot2-for-loop-plots-same-data
# https://stackoverflow.com/questions/49183067/trying-to-make-a-list-of-ggplot-objects-in-a-for-loop-all-items-in-list-are-wri

col_mult_adj_left = max(country_data$loss)/max(country_data$protected) # for color scale
brks = c(0.01,0.03,0.07,1)

bg_adj = expand.grid(protected=seq(0,
	                           1.05*max(country_data$protected,na.rm=TRUE),
                                   length=1e2),
	             loss=seq(0.95*min(country_data$loss,na.rm=TRUE),
                                    1.05*max(country_data$loss,na.rm=TRUE),
                                    length=1e2))

bg_adj$adj_loss = bg_adj$loss / bg_adj$protected
col_lim = atan(range(bg_adj$adj_loss)/col_mult_adj_left)

p_left = ggplot(country_data, aes(y=loss,x=protected)) +
	geom_raster(data=bg_adj, aes(fill=atan(adj_loss/col_mult_adj_left))) +
	geom_point() +
	scale_fill_gradientn(colours=plasma(50)[15:50],
                             breaks=atan(brks/col_mult_adj_left),
                             labels=brks,
                             limits=col_lim) +
	theme_bw() +
	xlab("(Proportion of forested area protected) * (PA effectiveness)") +
	theme(axis.text=element_text(color="black"),
              plot.margin=margin(3,3,3,20),
              legend.key.width=unit(2.5,"lines"),              
	      legend.position="bottom") +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	ylab("Annual deforestation rate (%)") +
        geom_text_repel(data=country_data[country_data$loss > 0.017 | country_data$protected > 1,], aes(label=NAME)) +
	guides(fill=guide_colorbar(title="Adjusted threat  \nindex (forest loss)", nbin=300))

###

col_mult_adj_right = max(country_data$log_carbon)/max(country_data$protected) # for color scale
brks = c(15,25,50,100)

bg_adj = expand.grid(protected=seq(0,
	                           1.05*max(country_data$protected,na.rm=TRUE),
                                   length=1e2),
	             log_carbon=seq(0.95*min(country_data$log_carbon,na.rm=TRUE),
                                    1.05*max(country_data$log_carbon,na.rm=TRUE),
                                    length=1e2))

bg_adj$adj_carbon = bg_adj$log_carbon / bg_adj$protected
col_lim = atan(range(bg_adj$adj_carbon)/col_mult_adj_right)

p_right = ggplot(country_data, aes(y=log_carbon,x=protected)) +
	geom_raster(data=bg_adj, aes(fill=atan(adj_carbon/col_mult_adj_right))) +
	geom_point() +
	scale_fill_gradientn(colours=plasma(50)[15:50],
                             breaks=atan(brks/col_mult_adj_right),
                             labels=brks,
                             limits=col_lim) +
	theme_bw() +
	xlab("(Proportion of forested area protected) * (PA effectiveness)") +
	theme(axis.text=element_text(color="black"),
              plot.margin=margin(3,3,3,20),
              legend.key.width=unit(2.5,"lines"),              
	      legend.position="bottom") +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	ylab("Total aboveground forest carbon (gt)") +
        geom_text_repel(data=country_data[country_data$log_carbon < 8.2 | 
                        country_data$log_carbon > 10 | country_data$protected > 1,], aes(label=NAME)) +
	guides(fill=guide_colorbar(title="Adjusted threat  \nindex (carbon)", nbin=300))

all = plot_grid(p_left, p_right, ncol=2, labels=c("A.","B."))

png("output/adj_loss_carbon.png", width=12, height=8, units="in", res=400)
all
dev.off()

pdf("output/adj_loss_carbon.pdf", width=12, height=8)
all
dev.off()

##################### write table

country_tab = data.frame(country_data)
country_tab$geometry = NULL
write.csv(country_tab, "output/country_tab.csv")

# paper text
dim(country_tab)
sum(country_tab$protected_raw >= 0.17)
mean(country_tab$protected_raw >= 0.17)
sum(country_tab$protected_raw >= 0.5)

country_tab %>%
        select(NAME, total, adj_threat) %>%
        filter(total >= 500) %>%
        arrange(desc(adj_threat))

country_tab %>%
        select(NAME, loss, protected_raw) %>%
        arrange(desc(loss))

country_tab %>%
        select(NAME, carbon, adj_carbon) %>%
        arrange(desc(carbon))

