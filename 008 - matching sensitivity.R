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
require(spmoran)
require(grid)
require(tidyverse)
require(plyr)
require(ggforce)
require(rbounds)

# https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

quants = function(x) data.frame(ymin=quantile(x,0.25),ymax=quantile(x,0.75))

project_dir = "/home/chrisgraywolf/shared/analysis/PA_matching/"
setwd(project_dir)

######################################## join info. from PA table

PA_df = read.csv("data_processed/data_matched_10_class.csv", as.is=TRUE, na.strings=c("","NaN","NA"))
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

saveRDS(PA_df, "data_processed/PA_df_10_class.RDS")

######################################## sensitivity analysis

df = readRDS("data_processed/PA_df.RDS") %>%
        dplyr::filter(group == "main" & matched) # 18,171 obs.

df_10 = readRDS("data_processed/PA_df_10_class.RDS") %>%
        dplyr::filter(group == "main" & matched) # 13,291 obs.

inside = df_10$PA_loss
outside = df_10$Control_loss
round(100*(1-mean(inside)/mean(outside)),3) # 42.74% compared to 41.12% w/ 5 classes

df = rbind(data.frame(df,n_class=5), data.frame(df_10,n_class=10))

developed = df[df$developed,]
developed$Continent = "Higher GDP"
developing = df[! df$developed,]
developing$Continent = "Lower GDP"
df_dev = rbind(df, developed, developing)

df_small = reshape2::melt(df_dev, measure.vars=c("PA_loss","Control_loss"))
df_small$variable = factor(df_small$variable,
	levels=rev(c("PA_loss","Control_loss")),
	labels=rev(c("PA deforestation rate","Control deforestation rate")))
df_small$Continent = factor(df_small$Continent, levels=unique(rev(c("Higher GDP","Lower GDP",sort(unique(df$Continent))))))
df_small$gp = paste(df_small$Continent,df_small$variable,df_small$cat_simple,df_small$n_class)

p = ggplot(df_small[table(df_small$gp)[df_small$gp] >= 10 & df_small$cat_simple != "Unknown",],
           aes(x=Continent, y=value*100, color=variable, shape=variable)) +
  stat_summary(fun.data = quants, geom = "errorbar",
               position=position_dodge(width=2/3),width=0) +
  stat_summary(fun.y=median, geom = "point",position=position_dodge(width=2/3)) +
  facet_grid(n_class ~ cat_simple) +
  scale_color_manual(values=brewer.pal(3,"Set1")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(color="black"),
        axis.text.y=element_text(face=c(rep("plain",6),rep("bold",2))),
        axis.title.y=element_blank(),
		strip.text.y=element_text(angle=0),
        panel.grid.major.x=element_line(color="lightgray"),
        legend.position="bottom",
        legend.box="vertical",
        legend.background=element_rect(color="black")) +
  ylab("Annual forest loss rate (%)") +
  guides(color=guide_legend(title=NULL,reverse=TRUE),
         shape=guide_legend(title=NULL,reverse=TRUE))

png("output/sensitivity_loss_box_no_unk.png", width=7, height=7, units="in", res=300)
p
dev.off()

### SVC model

######################################## set up data

df_mod = readRDS("data_processed/PA_df_10_class.RDS") %>%
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
        dplyr::select(lat,long,PA_loss,Control_loss,pop_dens,travel_time,age,GDP,cat_simple,GIS_AREA, threatened, non_threatened)

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
              maxiter=100) # around 2 min
Sys.time() - start

# https://rdrr.io/cran/spmoran/man/besf_vc.html
mod$s
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
                                      labels = scales::number_format(accuracy = 0.01),
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

png("output/sensitivity_svc.png", width=18, height=16, units="in", res=200)
all
dev.off()

######################################## Rosenbaum bounds

df_mod = readRDS("data_processed/PA_df.RDS") %>%
        filter(matched) %>%
        filter(group == "main") %>%
        mutate(PA_loss = log(1e-7 + PA_loss),
               Control_loss = log(1e-7 + Control_loss))

hlsens(df_mod$PA_loss, df_mod$Control_loss)

# Gamma Lower bound Upper bound
#     1     -1.9378   -1.937800
#     2     -4.4378   -0.637840
#     3     -5.3378   -0.237840
#     4     -5.7378    0.062158
#     5     -6.0378    0.262160
#     6     -6.2378    0.462160