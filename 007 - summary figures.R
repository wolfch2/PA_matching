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
require(sf)
require(smoothr)

quants = function(x) data.frame(ymin=quantile(x,0.25),ymax=quantile(x,0.75))
se_fun = function(x) data.frame(ymin=mean(x) - sd(x)/sqrt(length(x)),
                            ymax=mean(x) + sd(x)/sqrt(length(x)))
fmt_dcimals <- function(decimals=0){ function(x) formatC(x,format="fg",digits=1,big.mark=",")} # only need 1 signif. digit..

project_dir = "/home/chrisgraywolf/shared/analysis/PA_matching/"
setwd(project_dir)

######################################## set up data

df = readRDS("data_processed/PA_df.RDS")

######################################## numbers for paper text

### overall totals and proportions

addmargins(table(matched=df$matched, cat=df$cat_simple, gp=df$group))
round(prop.table(addmargins(table(matched=df$matched, cat=df$cat_simple, gp=df$group),2),2:3),3)

### adjust total % of forest protected (for paper text...)

countries = raster("data_processed/rasters/countries.tif")
cover = raster("data_processed/rasters/cover.tif") %>%
        mask(countries)
PA_rast = read_sf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp") %>%
        dplyr::filter(STATUS != "Proposed") %>%
        st_transform(st_crs(cover)) %>%
        fasterize(cover)

tab = table(cover[] >= 30, PA_rast[], useNA="always")
prop.table(tab,1) # 15.7% of forest is protected

df = readRDS("data_processed/PA_df.RDS") %>%
        dplyr::filter(group == "main" & matched)

# absolute terms - approximate total forest loss per year on average:
(100*sum(df$GIS_AREA * df$cover_loss/100)/18)/1e6 # (/100 to convert to proportion of pixel)

sum(df$GIS_AREA)

mean(df$PA_loss == 0)

inside = df$PA_loss
outside = df$Control_loss
round(100 * mean(inside), 3) # 0.620% mean deforestation rate inside PAs
round(100 * mean(outside), 3) # 1.053% mean deforestation rate inside control areas

### adjust total % of forest protected (for paper text...)
prop.table(tab,1)["TRUE","1"] * (1 - mean(inside)/mean(outside)) # 6.45% protected after adjusting
1 - mean(inside)/mean(outside) # reduction
0.17 / (1 - mean(inside)/mean(outside)) # would need 41.34% to match 17% w/ no loss

### simple comparison of rates
round(100*(1-mean(inside)/mean(outside)),3) # 41.12% reduction in deforestation rate
round(100*median(outside-inside),3) # median additive difference of 0.194% in annual deforestation rate
IQR(outside-inside)
median(outside)
sd(outside)
median(inside)
sd(inside)
median(outside)/median(inside)

######################################## summary figureumber of protected or control areass

developed = df[df$developed,]
developed$Continent = "Higher GDP"
developing = df[! df$developed,]
developing$Continent = "Lower GDP"
df_dev = rbind(df, developed, developing)

### main comparison (boxplot) (UPDATED)

df_small = reshape2::melt(df_dev, measure.vars=c("PA_loss","Control_loss"))
df_small$variable = factor(df_small$variable,
	levels=rev(c("PA_loss","Control_loss")),
	labels=rev(c("PA deforestation rate","Control deforestation rate")))
df_small$Continent = factor(df_small$Continent, levels=unique(rev(c("Higher GDP","Lower GDP",sort(unique(df$Continent))))))
df_small$gp = paste(df_small$Continent,df_small$variable,df_small$cat_simple)

# paper text
temp = df_small[table(df_small$gp)[df_small$gp] >= 10,]
temp = temp[temp$variable == "PA deforestation rate",]
temp$y = 100*temp$value
tab = melt(tapply(temp$y, list(temp$Continent, temp$cat_simple), function(x) quantile(x,0.5)))
tab[order(tab$value, decreasing=TRUE),]

p = ggplot(df_small[table(df_small$gp)[df_small$gp] >= 10 & df_small$cat_simple != "Unknown",],
           aes(x=Continent, y=value*100, color=variable, shape=variable)) +
  stat_summary(fun.data = quants, geom = "errorbar",
               position=position_dodge(width=2/3),width=0) +
  stat_summary(fun.y=median, geom = "point",position=position_dodge(width=2/3)) +
  facet_wrap(~ cat_simple, nrow=2) +
  scale_color_manual(values=brewer.pal(3,"Set1")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(color="black"),
        axis.text.y=element_text(face=c(rep("plain",6),rep("bold",2))),
        axis.title.y=element_blank(),
        panel.grid.major.x=element_line(color="lightgray"),
        legend.position="bottom",
        legend.box="vertical",
        legend.background=element_rect(color="black")) +
  ylab("Annual forest loss rate (%)") +
  guides(color=guide_legend(title=NULL,reverse=TRUE),
         shape=guide_legend(title=NULL,reverse=TRUE))

pdf("output/loss_box_no_unk.pdf", width=7, height=7)
p
dev.off()

png("output/loss_box_no_unk.png", width=7, height=7, units="in", res=300)
p
dev.off()

### PAPER TEXT (extra)
main = df_small[table(df_small$gp)[df_small$gp] >= 10,]
main$value = 100*main$value
(value = tapply(main$value, list(main$Continent,main$cat_simple,main$variable), median))
mean(df_small$value > 0.005)
mean(df_small$value > 0.01)

### net loss (updated)

df_small = melt(df_dev, measure.vars=c("PA_loss_gain","Control_loss_gain"))
df_small$variable = factor(df_small$variable,
	levels=rev(c("PA_loss_gain","Control_loss_gain")),
	labels=rev(c("Net forest loss inside PA","Net forest loss inside control pixels")))
df_small$Continent = factor(df_small$Continent, levels=unique(rev(c("Higher GDP","Lower GDP",sort(unique(df$Continent))))))
df_small$gp = paste(df_small$Continent,df_small$variable,df_small$cat_simple)

p = ggplot(df_small[table(df_small$gp)[df_small$gp] >= 10 & df_small$cat_simple != "Unknown",],
           aes(x=Continent, y=100*value, color=variable, shape=variable)) +
	stat_summary(fun.data = quants, geom = "errorbar",
		position=position_dodge(width=2/3),width=0) +
	stat_summary(fun.y=median, geom = "point",position=position_dodge(width=2/3)) +
        facet_wrap(~ cat_simple, ncol=1) +
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
	ylab("Net annual forest loss rate (%)") +
	guides(color=guide_legend(title=NULL,reverse=TRUE),
		shape=guide_legend(title=NULL,reverse=TRUE))

pdf("output/loss_box_net.pdf", width=7, height=7)
p
dev.off()

png("output/loss_box_net.png", width=7, height=7, units="in", res=300)
p
dev.off()

### change comparison (boxplot)

df = readRDS("data_processed/PA_df.RDS") %>%
        dplyr::filter(group == "change" & matched) %>%        
        mutate(loss_delta_inside = PA_loss_after - PA_loss_before,
               loss_delta_outside = Control_loss_after - Control_loss_before) %>%
        melt(measure.vars=c("loss_delta_inside","loss_delta_outside")) %>%
        mutate(variable = factor(variable, 
                                 levels=rev(c("loss_delta_inside","loss_delta_outside")),
	                         labels=rev(c("Protected areas","Control areas"))),
               Continent = factor(Continent,
                                  levels=unique(rev(c("Higher GDP","Lower GDP",sort(unique(Continent)))))),
               gp = paste(Continent,variable,cat_simple))

developed = df[df$developed,]
developed$Continent = "Higher GDP"
developing = df[! df$developed,]
developing$Continent = "Lower GDP"
df_small = rbind(df, developed, developing)

p = ggplot(df_small[table(df_small$gp)[df_small$gp] >= 10 & df_small$cat_simple != "Unknown",],
		aes(x=Continent, y=100*value, color=variable, shape=variable)) +
        geom_hline(yintercept=0, linetype="dashed", alpha=0.5) +
	stat_summary(fun.data=se_fun, geom = "errorbar",
		position=position_dodge(width=2/3),width=0) +
	stat_summary(fun.y=mean, geom = "point",position=position_dodge(width=2/3)) + 
        facet_wrap(~ cat_simple, ncol=1) +
	scale_color_manual(values=brewer.pal(3,"Set1")) +
	coord_flip() +
	theme_bw() +
	theme(axis.text=element_text(color="black"),
	      	axis.text.y=element_text(face=c(rep("plain",6),rep("bold",2))),
		axis.title.y=element_blank(),
		panel.grid.major.x=element_line(color="lightgray"),
		legend.position="bottom",
		legend.box="vertical",
		legend.background=element_rect(color="black")) +
	ylab("Additive change in annual forest loss rate (%)") +
	guides(color=guide_legend(title=NULL,reverse=TRUE),
	       shape=guide_legend(title=NULL,reverse=TRUE))

pdf("output/change_comp.pdf", width=7, height=7)
p
dev.off()

png("output/change_comp.png", width=7, height=7, units="in", res=400)
p
dev.off()

# paper text
round(100*mean(df$PA_loss_before), 3) 
round(100*mean(df$PA_loss_after), 3) 
P = df$PA_loss_after - df$PA_loss_before
C = df$Control_loss_after - df$Control_loss_before
round(100*mean(P), 3)
round(100*mean(C), 3)
100*sd(P)/sqrt(length(P))
100*sd(C)/sqrt(length(C))
round(100*mean(df$Control_loss_before), 3) 
round(100*mean(df$Control_loss_after), 3) 
round(100*tapply(df_small$value, list(df_small$variable, df_small$Continent), mean),3)

### summarize basic discrete data (UPDATED)

df_long = readRDS("data_processed/PA_df.RDS") %>%
        dplyr::filter(group == "main") %>%
        dplyr::filter(matched) %>%
        dplyr::count(IUCN_CAT, cat_simple, Continent) %>%
        mutate(IUCN_CAT = factor(IUCN_CAT,
                                 levels=rev(c("Ia","Ib","II","III","IV","V","VI",
                                          "Not Reported","Not Assigned","Not Applicable"))),
               cat_simple = factor(cat_simple, levels=c("Strict","Nonstrict","Unknown")),
               Continent = factor(Continent, levels=rev(sort(unique(Continent)))))

p = ggplot(df_long, aes(x=IUCN_CAT, y=n/sum(df_long$n), fill=Continent)) +
	geom_bar(stat="identity", position="stack") +
        facet_grid(cat_simple ~ ., scales="free", space="free") +
	scale_fill_manual(values=rev(brewer.pal(6,"Set1"))) +
	coord_flip() +
	theme_bw() +
	theme(axis.text=element_text(color="black"),
		axis.title.y=element_blank(),
		strip.text.y=element_text(angle=0),
		legend.background=element_rect(color="black"),
		legend.position="bottom") +
	scale_y_continuous(labels=percent) +
	ylab("Percentage of protected areas") +
	guides(fill=guide_legend(reverse=TRUE,title="Continent",nrow=2,byrow=TRUE))

pdf("output/discrete.pdf", width=6, height=6)
p
dev.off()

png("output/discrete.png", width=6, height=6, units="in", res=400)
p
dev.off()

### loss histograms

df = readRDS("data_processed/PA_df.RDS") %>%
        dplyr::filter(group == "main" & matched) %>%
        select(cat_simple, Continent, PA_loss, Control_loss) %>%
        melt(id.vars=c("cat_simple","Continent")) %>%
        mutate(variable=factor(variable,
                               levels=c("PA_loss","Control_loss"),
                               labels=c("Protected areas","Control areas")))

mean(df$value < 0.05) # cover more than 98% by stopping at 5% loss rate

p = ggplot(df, aes(x=100*value,y=..count..,fill=Continent)) +
	geom_histogram(boundary=0,color="black") +
	scale_fill_manual(values=brewer.pal(6,"Set1")) +
	facet_grid(cat_simple ~ variable, scales="free_y") +
	ylab("Count") +
	xlab("Annual deforestation rate (%)") +
	theme_bw() +
	theme(axis.text=element_text(color="black"),
		legend.background=element_rect(color="black"),
		legend.position="bottom") +
	guides(fill=guide_legend(title=element_blank(), nrow=1)) +
        scale_x_continuous(trans="log1p", labels = fmt_dcimals(0), limits=c(0,5), breaks=c(0,1,2,5,10,25,50,100))

pdf("output/loss_hist.pdf", width=6.5, height=6.5)
p
dev.off()

png("output/loss_hist.png", width=6.5, height=6.5, units="in", res=300)
p
dev.off()

p = ggplot(df, aes(x=100*value,y=..count..)) +
	geom_histogram(boundary=0,color="black",fill="lightgray") +
	facet_grid(cat_simple ~ variable, scales="free_y") +
	ylab("Number of protected or control areas") +
	xlab("Annual deforestation rate (%)") +
	theme_bw() +
	theme(axis.text=element_text(color="black"),
		legend.background=element_rect(color="black"),
		legend.position="bottom") +
	guides(fill=guide_legend(title=element_blank(), nrow=1)) +
        scale_x_continuous(trans="log1p", labels = fmt_dcimals(0), limits=c(0,5), breaks=c(0,1,2,5,10,25,50,100))

pdf("output/loss_hist_bw.pdf", width=6.5, height=6.5)
p
dev.off()

png("output/loss_hist_bw.png", width=6.5, height=6.5, units="in", res=300)
p
dev.off()

# paper text
100*tapply(df$value, df$variable, median)
100*tapply(df$value, df$variable, mean)

### basic PA map

worldmap = ne_download(scale = 50,
                       type = "land",
                       category = "physical",
                       destdir = tempdir(),
                       load = TRUE,
                       returnclass = "sf") %>%
        st_transform(54030) %>%
        st_union

df = readRDS("data_processed/PA_df.RDS") %>%
        dplyr::filter(group == "main" & matched)

PAs = read_sf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp") %>%
        filter(WDPAID %in% df$WDPAID) %>%
        st_transform(54030) %>%
        st_simplify(dTolerance=5000) %>%
        mutate(cat_simple = df$cat_simple[match(WDPAID, df$WDPAID)]) %>%
        group_by(cat_simple)

PAs = dplyr::summarise(PAs) # had issues with piping

border = raster() %>%
        st_bbox %>%
        st_as_sfc %>%
        densify(100) %>%
        st_transform(54030)

gdalwarp(srcfile="data_processed/rasters/cover.tif",
         dstfile="data_processed/cover_low.tif",
         t_srs="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
         tr=c(10000,10000))

gdalwarp(srcfile="data_processed/rasters/cover_loss.tif",
         dstfile="data_processed/cover_loss_low.tif",
         t_srs="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
         tr=c(10000,10000))

C = raster("data_processed/cover_low.tif")
CL = raster("data_processed/cover_loss_low.tif")
D = -(((C - CL)/C)^(1/18) - 1)
D[which(D[] > 0.1)] = 0.1
D[which(C[] < 30)] = NA
D_pts = data.frame(rasterToPoints(D))

p_top = ggplot() +
        geom_sf(data=PAs, aes(fill=cat_simple), color=NA) +
        geom_sf(data=worldmap, fill=NA, color="black", size=0.2) +
        geom_sf(data=border, fill=NA, color="black", size=0.2) +
        scale_x_continuous(limits=st_bbox(PAs)[c(1,3)], expand=c(0.05,0.05)) +
        scale_y_continuous(limits=st_bbox(PAs)[c(2,4)], expand=c(0.05,0.05)) +
        scale_fill_manual(values=brewer.pal(9,"Set1")[c(2,4,9)],
                           guide=guide_legend(title="Reserve category",
                                              override.aes=list(color="black"))) +
        theme_bw() +
        theme(axis.title=element_blank(),
                panel.grid.major=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                plot.margin=margin(5,2,5,23),
                plot.title=element_text(hjust=0.5),
                legend.position=c(0,0),
                legend.justification=c(0,0),
                legend.background=element_rect(color="black",fill="white"))

p_bot = ggplot() +
        geom_raster(data=D_pts, aes(x=x,y=y,fill=100*layer)) +
        geom_sf(data=worldmap, fill=NA, color="black", size=0.2) +
        geom_sf(data=border, fill=NA, color="black", size=0.2) +
        scale_x_continuous(limits=st_bbox(PAs)[c(1,3)], expand=c(0.05,0.05)) +
        scale_y_continuous(limits=st_bbox(PAs)[c(2,4)], expand=c(0.05,0.05)) +
        scale_fill_gradientn(colors=brewer.pal(9,"Reds"),
                             breaks=seq(0,10,length=5),
                             labels=paste0(c(rep("",4),"≥"), seq(0,10,length=5)),
                             guide=guide_colorbar(title="Annual\ndeforestation\nrate (%)",
                                                  nbin=300)) +
        theme_bw() +
        theme(axis.title=element_blank(),
                panel.grid.major=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                plot.margin=margin(5,2,5,23),                
                plot.title=element_text(hjust=0.5),
                legend.position=c(0,0),
                legend.justification=c(0,0),
                legend.background=element_rect(color="black",fill="white"))

all = plot_grid(p_top, p_bot, ncol=1, labels=c("A.","B."))

png("output/PA_map.png", width=.9*10, height=.9*9, units="in", res=300)
all
dev.off()

pdf("output/PA_map.pdf", width=.9*10, height=.9*9)
all
dev.off()

