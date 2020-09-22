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
              #covmodel="gau", # fails to run w/ "exp"! // UPDATE -> works in new version
              maxiter=100) # ~20 sec. // UPDATE -> around 3-4 minutes
Sys.time() - start

# https://rdrr.io/cran/spmoran/man/besf_vc.html
mod$s # note second row: all below 0.15 except for intercept (extremely high)
mod$vc # everything but GDP effect was varying

rbind(mod$b_vc[1,],mod$bse_vc[1,],mod$p_vc[1,]) # note results for GDP

coef = mod$b_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%
        select(- c(travel_time, X.Intercept.)) %>%
        melt(id.vars=c("lat","long")) %>%
        mutate(variable = revalue(variable, c("Control_loss"="Background rate",
                                              "pop_dens"="Population density",
                                              "age"="Reserve age",
                                              "GDP"="GDP per capita",                                              
                                              "cat_simpleStrict"="Strict protection",
                                              "GIS_AREA"="Reserve area")))

SE = mod$bse_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%
        select(- c(travel_time, X.Intercept.)) %>%
        melt(id.vars=c("lat","long")) %>%
        mutate(variable = revalue(variable, c("Control_loss"="Background rate",
                                              "pop_dens"="Population density",
                                              "age"="Reserve age",
                                              "GDP"="GDP per capita",                                              
                                              "cat_simpleStrict"="Strict protection",
                                              "GIS_AREA"="Reserve area")))

p_val = mod$p_vc %>%
        data.frame(lat=df_mod$lat, long=df_mod$long) %>%        
        select(- c(travel_time, X.Intercept.)) %>%
        reshape2::melt(id.vars=c("lat","long")) %>%
        mutate(variable = revalue(variable, c("Control_loss"="Background rate",
                                              "pop_dens"="Population density",
                                              "age"="Reserve age",
                                              "GDP"="GDP per capita",                                              
                                              "cat_simpleStrict"="Strict protection",
                                              "GIS_AREA"="Reserve area")),
              value = p.adjust(value,"fdr"))

# https://stackoverflow.com/questions/55922441/expand-argument-in-scale-color-gradient-is-ignored
p_ramp = c(rev(colorRampPalette(brewer.pal(9,"Blues"))(1e3)),
           rev(colorRampPalette(brewer.pal(9,"Reds"))(1e3)))

plot_list = lapply(unique(coef$variable), function(var){
        p_0 = textGrob(var,gp=gpar(fontsize=20), rot=90)
        u = max(abs(range(coef[coef$variable == var,]$value)))
        p_1 = ggplot(coef[coef$variable == var & p_val[p_val$variable == var,"value"] < 0.05,]) +
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

### simplified version for main paper

var_text = lapply(c("","Background rate","Reserve area"), function(x) textGrob(x,gp=gpar(fontsize=20)))
type_text = lapply(c("Coefficient","Standard error","p-value"), function(x) textGrob(x,gp=gpar(fontsize=20),rot=90))

plot_list_simple = c(var_text, type_text[1], plot_list[c(2,22)], type_text[2], plot_list[c(3,23)], type_text[3], plot_list[c(4,24)])

all_small = plot_grid(plotlist=plot_list_simple,ncol=3,rel_heights=c(0.9,rep(5,3)),rel_widths=c(0.2,rep(3,2)))

png("output/svc_small.png", height=7.5, width=12, units="in", res=500)
all_small
dev.off()

pdf("output/svc_small.pdf", height=7.5, width=12)
all_small
dev.off()

######################################## species richness models

set.seed(0)
mod_threatened = besf_vc(df_mod$PA_loss, x=model.matrix(PA_loss ~ Control_loss + threatened, data=df_mod),
              coords=cbind(df_mod$lat, df_mod$long),
              # covmodel="gau",
              maxiter=100)
set.seed(0)
mod_non_threatened = besf_vc(df_mod$PA_loss, x=model.matrix(PA_loss ~ Control_loss + non_threatened, data=df_mod),
              coords=cbind(df_mod$lat, df_mod$long),
              # covmodel="gau",
              maxiter=100)

mod_threatened$vc # not SV
mod_threatened$b_vc$threatened[1] # coef
mod_threatened$bse_vc$threatened[1] # SE
mod_threatened$p_vc$threatened[1] # p-val

mod_non_threatened$vc # not SV
mod_non_threatened$b_vc$non_threatened[1] # coef
mod_non_threatened$bse_vc$non_threatened[1] # SE
mod_non_threatened$p_vc$non_threatened[1] # p-val
