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
require(ggeasy) # remotes::install_github("jonocarroll/ggeasy")

project_dir = "/home/chrisgraywolf/shared/analysis/PA_matching/"
setwd(project_dir)

######################################## set up data

df_mod = readRDS("data_processed/PA_df.RDS") %>%
        filter(matched) %>%
        filter(group == "main") %>%
        filter(cat_simple %in% c("Strict","Nonstrict")) %>%
        filter(STATUS_YR > 0) %>%
        mutate(age = 2020 - STATUS_YR,
               travel_time = log(1 + travel_time),
               GDP = log(GDP),
               GIS_AREA = log(GIS_AREA),
               PA_loss = log(1e-7 + PA_loss),
               Control_loss = log(1e-7 + Control_loss)) %>%
        select(lat,long,PA_loss,Control_loss,pop_dens,travel_time,age,GDP,cat_simple,GIS_AREA, threatened, non_threatened)

worldmap = ne_download(scale = 110,
                       type = "countries",
                       category = "cultural",
                       destdir = tempdir(),
                       load = TRUE,
                       returnclass = "sf")

######################################## fit main model, plot results

X = model.matrix(PA_loss ~ Control_loss + pop_dens + travel_time + age + GDP + cat_simple + GIS_AREA, data=df_mod)

start = Sys.time()
set.seed(0)
mod = besf_vc(df_mod$PA_loss, x=X,
              coords=cbind(df_mod$lat, df_mod$long),
              covmodel="gau", # fails to run w/ "exp"!
              maxiter=100) # ~20 sec.
Sys.time() - start

# https://rdrr.io/cran/spmoran/man/besf_vc.html
mod$s # note second row: all below 0.15 except for intercept (extremely high)
mod$vc # everything but GDP effect was varying

rbind(mod$b_vc[1,],mod$bse_vc[1,],mod$p_vc[1,]) # note results for GDP

coef = mod$b_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%
        select(- c(GDP, X.Intercept.)) %>%
        melt(id.vars=c("lat","long")) %>%
        mutate(variable = revalue(variable, c("Control_loss"="Background rate",
                                              "pop_dens"="Population density",
                                              "travel_time"="Travel time",
                                              "age"="Reserve age",
                                              "cat_simpleStrict"="Strict protection",
                                              "GIS_AREA"="Reserve area")))

SE = mod$bse_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%
        select(- c(GDP, X.Intercept.)) %>%
        melt(id.vars=c("lat","long")) %>%
        mutate(variable = revalue(variable, c("Control_loss"="Background rate",
                                              "pop_dens"="Population density",
                                              "travel_time"="Travel time",
                                              "age"="Reserve age",
                                              "cat_simpleStrict"="Strict protection",
                                              "GIS_AREA"="Reserve area")))

p_val = mod$p_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%        
        select(- c(GDP, X.Intercept.)) %>%
        reshape2::melt(id.vars=c("lat","long")) %>%
        mutate(variable = revalue(variable, c("Control_loss"="Background rate",
                                              "pop_dens"="Population density",
                                              "travel_time"="Travel time",
                                              "age"="Reserve age",
                                              "cat_simpleStrict"="Strict protection",
                                              "GIS_AREA"="Reserve area")),
              value = p.adjust(value,"fdr"))

# https://stackoverflow.com/questions/55922441/expand-argument-in-scale-color-gradient-is-ignored
p_ramp = c(rev(colorRampPalette(brewer.pal(9,"Blues"))(1e3)),
           rev(colorRampPalette(brewer.pal(9,"Reds"))(1e3)))

plot_list = lapply(unique(coef$variable), function(var){
        p_0 = textGrob(var,gp=gpar(fontsize=20), rot=90)
        u = max(abs(range(coef[coef$variable == var,]$value)))
        p_1 = ggplot(coef[coef$variable == var,]) +
                geom_point(aes(x=long,y=lat,color=value)) +
                scale_color_gradientn(colors=rev(brewer.pal(11,"PiYG")),
                                      limits=c(-u,u),
                                      guide=guide_colorbar(title=NULL,
                                                           nbin=1000)) +
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
                scale_y_continuous(limits=range(df_mod$lat), expand=c(0.05,0.05))
        p_2 = ggplot(SE[SE$variable == var,]) +
                geom_point(aes(x=long,y=lat,color=value)) +
                scale_color_gradientn(colors=brewer.pal(9,"BuPu"),
                                      guide=guide_colorbar(title=NULL,
                                                           nbin=1000)) +
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
                scale_y_continuous(limits=range(df_mod$lat), expand=c(0.05,0.05))
        p_3 = ggplot(p_val[p_val$variable == var,]) +
                geom_point(aes(x=long,y=lat,color=value)) +
                scale_color_gradientn(colors=p_ramp,
                                      values=c(seq(0,0.05^(1/2),length=1e3),seq(0.05^(1/2),1,length=1e3)),
                                      trans=power_trans(1/2),
                                      limits=c(0,1),
                                      breaks=c(0,0.05,0.25,0.5,1),
                                      guide=guide_colorbar(title=NULL,
                                                           nbin=1000)) +
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
                scale_y_continuous(limits=range(df_mod$lat), expand=c(0.05,0.05))
        return(list(p_0, p_1, p_2, p_3))
})
plot_list = unlist(plot_list, recursive=FALSE)

plot_list_titles = c(lapply(c("","Coefficient","Standard error","p-value"), function(x) textGrob(x,gp=gpar(fontsize=20))),
              plot_list)

all = plot_grid(plotlist=plot_list_titles,ncol=4,rel_widths=c(0.6,rep(5,3)),rel_heights=c(0.35,rep(3,length(unique(coef$variable)))))

png("output/svc.png", width=18, height=16, units="in", res=200)
all
dev.off()

pdf("output/svc.pdf", width=18, height=16)
all
dev.off()

plot_list_single = list(plot_list_titles[[22]] + ggtitle("Coefficient") + easy_center_title(),
                        plot_list_titles[[23]] + ggtitle("Standard error") + easy_center_title(),
                        plot_list_titles[[24]] + ggtitle("p-value") + easy_center_title())  # pull out category plots for main paper

single = plot_grid(plotlist=plot_list_single,ncol=1)

png("output/svc_single.png", width=1.1*5, height=7, units="in", res=500)
single
dev.off()

pdf("output/svc_single.pdf", width=1.1*5, height=7)
single
dev.off()

######################################## species richness models

set.seed(0)
mod_threatened = besf_vc(df_mod$PA_loss, x=model.matrix(PA_loss ~ Control_loss + threatened, data=df_mod),
              coords=cbind(df_mod$lat, df_mod$long),
              covmodel="gau",
              maxiter=100)
set.seed(0)
mod_non_threatened = besf_vc(df_mod$PA_loss, x=model.matrix(PA_loss ~ Control_loss + non_threatened, data=df_mod),
              coords=cbind(df_mod$lat, df_mod$long),
              covmodel="gau",
              maxiter=100)

mod_threatened$vc # not SV
mod_threatened$b_vc$threatened[1] # coef
mod_threatened$bse_vc$threatened[1] # SE
mod_threatened$p_vc$threatened[1] # p-val

mod_non_threatened$vc

coef = mod_non_threatened$b_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%
        select(- c(Control_loss, X.Intercept.))
SE = mod_non_threatened$bse_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%
        select(- c(Control_loss, X.Intercept.))
p_val = mod_non_threatened$p_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%
        select(- c(Control_loss, X.Intercept.)) %>%
        mutate(non_threatened=p.adjust(non_threatened,"fdr"))

u = coef$non_threatened %>% range %>% abs %>% max

p_coef = ggplot(coef) +
                geom_point(aes(x=long,y=lat,color=non_threatened)) +
                scale_color_gradientn(colors=rev(brewer.pal(11,"PiYG")),
                                      limits=c(-u,u),
                                      guide=guide_colorbar(title="Coefficient",
                                                           nbin=1000)) +
                geom_sf(data=worldmap, fill=NA, color="black") +
                theme_bw() +
                theme(axis.title=element_blank(),
                      panel.grid.major=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      plot.title=element_text(hjust=0.5),
                      legend.position="right",
                      legend.background=element_rect(color="black",fill="white")) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=range(df_mod$lat), expand=c(0.05,0.05))

p_SE = ggplot(SE) +
                geom_point(aes(x=long,y=lat,color=non_threatened)) +
                scale_color_gradientn(colors=brewer.pal(9,"BuPu"),
                                      guide=guide_colorbar(title="Standard\nerror",
                                                           nbin=1000)) +
                geom_sf(data=worldmap, fill=NA, color="black") +
                theme_bw() +
                theme(axis.title=element_blank(),
                      panel.grid.major=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      plot.title=element_text(hjust=0.5),
                      legend.position="right",
                      legend.background=element_rect(color="black",fill="white")) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=range(df_mod$lat), expand=c(0.05,0.05))

p_p_val = ggplot(p_val) +
                geom_point(aes(x=long,y=lat,color=non_threatened)) +
                scale_color_gradientn(colors=p_ramp,
                                      values=c(seq(0,0.05^(1/2),length=1e3),seq(0.05^(1/2),1,length=1e3)),
                                      trans=power_trans(1/2),
                                      limits=c(0,1),
                                      breaks=c(0,0.05,0.25,0.5,1),
                                      guide=guide_colorbar(title="p-value",
                                                           nbin=1000)) +
                geom_sf(data=worldmap, fill=NA, color="black") +
                theme_bw() +
                theme(axis.title=element_blank(),
                      panel.grid.major=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      plot.title=element_text(hjust=0.5),
                      legend.position="right",
                      legend.background=element_rect(color="black",fill="white")) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(limits=range(df_mod$lat), expand=c(0.05,0.05))

p_non_threatened = plot_grid(plotlist=list(p_coef,p_SE,p_p_val), ncol=1, align="v")

png("output/non_threatened.png", width=7, height=7, units="in", res=300)
p_non_threatened
dev.off()

pdf("output/non_threatened.pdf", width=7, height=7)
p_non_threatened
dev.off()

### map species richness for context

titles = c("threatened"="Threatened forest\nspecies richness",
           "non_threatened"="Non-threatened forest\nspecies richness")

plot_list = lapply(names(titles), function(rast_name){
        rast = velox(raster(paste0("data_processed/rasters/",rast_name,".tif")))
        rast$aggregate(10,aggtype="median")
        rast = rast$as.RasterLayer()
        rast[rast[] %in% 0] = NA
        rast = projectRaster(rast, res=0.1, crs="+proj=longlat +datum=WGS84 +no_defs", method="ngb")
        pts = data.frame(rasterToPoints(rast))
        
        p = ggplot(pts, aes(x=x,y=y,fill=layer)) +
                geom_raster() +
                scale_fill_gradientn(colors=brewer.pal(11,"Spectral"),
                                     limits=c(0,NA),
                                     guide=guide_colorbar(title=titles[rast_name],
                                                           nbin=1000)) +
                geom_sf(data=worldmap, aes(x=NULL,y=NULL), fill=NA, color="black", size=0.2) +
                theme_bw() +
                theme(axis.title=element_blank(),
                      panel.grid.major=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      plot.title=element_text(hjust=0.5),
                      legend.position=c(0,0),
                      legend.justification=c(0,0),
                      legend.background=element_rect(color="black",fill="white")) +
                coord_sf(xlim=c(-180,180), ylim=c(-55.95,72.75)) # extents of non-threatened
        return(p)
})

p_richness = plot_grid(plotlist=plot_list[2:1], ncol=1, align="hv")

png("output/richness.png", width=9.5, height=7, units="in", res=400)
p_richness
dev.off()

pdf("output/richness.pdf", width=9.5, height=7)
p_richness
dev.off()

