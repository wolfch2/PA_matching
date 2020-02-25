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

project_dir = "/home/chrisgraywolf/analysis_desktop/PA_matching/"
setwd(project_dir)

######################################## set up data

df = readRDS("data_processed/PA_df.RDS")

df_mod = na.omit(df[,c("loss","control_loss",
	"GDP","lat","long","cat_simple","GIS_AREA"),])
df_mod$abs_lat = abs(df_mod$lat); df_mod$log_GDP = log(df_mod$GDP, base=10); df_mod$log_GIS_AREA = log(df_mod$GIS_AREA, base=10)
df_mod = df_mod[df_mod$cat_simple %in% c("Strict","Nonstrict"),]

worldmap = ne_download(scale = 110,
                       type = "countries",
                       category = "cultural",
                       destdir = tempdir(),
                       load = TRUE,
                       returnclass = "sf")

######################################## gwr model

worldmap$value = 1
land_stars  = st_rasterize(worldmap["value"], deltax=1, deltay=1, options = "ALL_TOUCHED=TRUE")
land_pts = st_as_sf(st_xy2sfc(land_stars, as_points = TRUE))

PAs = st_as_sf(df_mod, coords = c('long', 'lat'), crs=4326)
PAs_union = st_union(PAs)

registerDoMC(24)
land_pts$dist = unlist(foreach(i=1:nrow(land_pts)) %dopar% {as.numeric(st_distance(land_pts[i,],PAs_union))})
land_pts = land_pts[land_pts$dist / 1e3 < 250,]

mod = gwr(loss ~ control_loss + log_GDP + cat_simple + log_GIS_AREA,
    data=df_mod,
    adapt=0.05,
    fit.points=st_coordinates(land_pts),
    coords=cbind(df_mod$long,df_mod$lat),
    longlat=TRUE)

mod_fitted = mod$SDF[,-c(1,2)] %>%
        data.frame %>%
        select(- optional) %>%
        melt(id.vars=c("X","Y")) %>%
        mutate(variable = revalue(variable, c("control_loss"="Deforestation rate\n(control)",
                                              "log_GDP"="GDP per capita\n(log transformed)",
                                              "log_GIS_AREA"="Reserve area in km2\n(log transformed)",
                                              "cat_simpleStrict"="Strict protection")))

plot_list = lapply(split(mod_fitted, mod_fitted$variable), function(fitted){
        #clamp = quantile(fitted$value, c(0.25,0.75))
        #fitted$value[fitted$value < clamp[1]] = clamp[1]
        #fitted$value[fitted$value > clamp[2]] = clamp[2]

        ggplot(fitted) +
                geom_raster(aes(x=X,y=Y,fill=value)) +
                scale_fill_gradientn(colors=magma(20),
                                     guide=guide_colorbar(title=fitted$variable[1]),
                                     values=rescale(quantile(fitted$value,seq(0,1,length=1e2)))) +
                geom_sf(data=worldmap, fill=NA, color="black") +
                theme_bw() +
                theme(axis.title=element_blank(),
                      panel.grid.major=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      plot.title=element_text(hjust=0.5),
                      legend.position=c(0.00,0.00),
                      legend.justification=c(0,0),
                      legend.background=element_rect(color="black",fill="white")) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=range(fitted$Y), expand=c(0,0))
})

all = plot_grid(plotlist=plot_list)

pdf("output/gwr.pdf", width=17, height=6*1.1)
all
dev.off()

