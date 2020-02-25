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

project_dir = "/home/chrisgraywolf/analysis_desktop/PA_matching/"
setwd(project_dir)

########################################

new = read.csv("data_processed/data_matched.csv", as.is=TRUE)
PA_table = read.csv("data_processed/PAs_tab.csv", as.is=TRUE)
PA_table = PA_table[PA_table$ID %in% new$PAs,]

new$GIS_AREA = tapply(PA_table$GIS_AREA, PA_table$ID, sum)[as.character(new$PAs)]
new$IUCN_CAT = tapply(PA_table$IUCN_CAT, PA_table$ID, Mode)[as.character(new$PAs)]
new$GOV_TYPE = tapply(PA_table$GOV_TYPE, PA_table$ID, Mode)[as.character(new$PAs)]
new$ISO3 = tapply(PA_table$PARENT_ISO, PA_table$ID, Mode)[as.character(new$PAs)]
new$NAME = tapply(PA_table$NAME, PA_table$ID, Mode)[as.character(new$PAs)]
new$WDPAID = PA_table$WDPAID[match(new$PAs, PA_table$ID)]

PAs = read_sf("data_input/PAs/WDPA_Jan2020-shapefile-polygons.shp")
PAs = PAs[match(new$WDPAID, PAs$WDPAID),]
PAs_proj = st_transform(PAs, "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
centroids = st_centroid(PAs_proj)
coords_Behr = st_coordinates(centroids)
coords_latlong = st_coordinates(st_transform(centroids,4326))

new$X = coords_Behr[,"X"]
new$Y = coords_Behr[,"Y"]
new$long = coords_latlong[,"X"]
new$lat = coords_latlong[,"Y"]

# new = merge(new, PA_table, by.x="PAs", by.y="ID", all.x=TRUE, all.y=FALSE)
# new$PAs = PA_table$WDPAID[match(new$PAs,PA_table$ID)]

########################################

worldmap = ne_download(scale = 10,
                       type = "countries",
                       category = "cultural",
                       destdir = tempdir(),
                       load = TRUE,
                       returnclass = "sf")

p = ggplot(new, aes(x=long, y=lat, color=n_control)) +
        geom_point(size=0.002, shape=16) +
        coord_fixed() +
        scale_color_gradientn(colors=magma(20), trans="log", breaks=10^(-1:6), labels=function(x) formatC(x, format="fg", big.mark=","),
                              guide=guide_legend(title="Control points",nrow=1,override.aes=list(size=2))) +
        geom_sf(data=worldmap, aes(x=NULL,y=NULL), fill=NA, color="black", size=0.2) +
        scale_x_continuous(limits=range(new$long)) +
        scale_y_continuous(limits=range(new$lat)) +        
        theme_bw() +
        theme(axis.ticks=element_blank(),
              axis.title=element_blank(),
              legend.position="bottom",
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.text=element_blank()) 

png("output/control_map.png", width=7, height=2.75, units="in", res=450)
p
dev.off()

######################################## basic summary

round(100*mean(new$loss),2)
round(100*mean(new$control_loss),2)
mean(new$loss) / mean(new$control_loss)

######################################## continent etc.

### ADD CONTINENT AND COUNTRY
df = new
df$ISO3_single = pbsapply(as.character(df$ISO3), function(x){
	if(nchar(x) == 3) return(x)
	out = names(sort(table(str_split(x,";")[[1]]),decreasing=TRUE)[1])
	return(out)
})
df$Continent = countrycode(df$ISO3_single, "iso3c", "continent")
df$Continent[df$ISO3_single == "TWN"] = "Asia" # fix a bunch of NA
df$Continent[countrycode(df$ISO3_single, "iso3c", "region")
	%in% c("South America")] = "South America"
df$Continent[df$Continent == "Americas"] = "North America"
table(df$ISO3_single[is.na(df$Continent)])
# ^-- just little islands and AA left: https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3
df$Country = countrycode(df$ISO3_single, "iso3c", "country.name")

### ADD GDP (perhaps mean of 2001-2016 GDP/PC is better given dates of loss, or maybe just 2001 GDP/PC?); also, can have multiple countries sep. by ";" (for ISO3, not ISO3_single)
# wbsearch(pattern = "GDP per capita")
# wbsearch("GDP")
GDP_raw = wb(indicator="NY.GDP.PCAP.PP.CD", startdate=2000, enddate=2018) # GDP per capita, PPP (current international $)
GDP = aggregate(value ~ iso3c, GDP_raw, function(x) mean(x,na.rm=TRUE))
# data.frame(iso3c=c("TWN","ALA"), value=(GDP$value[GDP$iso3c == "CHN"],GDP$value[GDP$iso3c == "FIN"]))
df$GDP = GDP[match(df$ISO3_single, GDP$iso3c),"value"]
# table(is.na(df$GDP),df$space_sel+1) # need to fix 86 PAs..; or not -- it's just 86/25267=0.3% and we'd like our GDP numbers to match those in covar map!
# sort(table(df$ISO3_single[is.na(df$GDP) & df$space_sel]))

### ADD SIMPLIFIED CAT; update: also add cat_pure and indigenous
df$cat_simple = 'Unknown'
df$cat_simple[df$IUCN_CAT %in% c('Ia','Ib','II','III','IV')] = 'Strict'
df$cat_simple[df$IUCN_CAT %in% c('V','VI')] = 'Nonstrict'
df$cat_simple[df$GOV_TYPE %in% 'Indigenous peoples'] = 'Indigenous'
df$cat_pure = 'Unknown'
df$cat_pure[df$IUCN_CAT %in% c('Ia','Ib','II','III','IV')] = 'Strict'
df$cat_pure[df$IUCN_CAT %in% c('V','VI')] = 'Nonstrict'
df$Indigenous = 'Not Indigenous'
df$Indigenous[df$GOV_TYPE %in% 'Indigenous peoples'] = "Indigenous"
df$Indigenous[df$GOV_TYPE %in% 'Not Reported'] = "Unknown"

saveRDS(df, "output/df_modeling.RDS")

######################################## basic plot TODO - investigate why so few Indigenous reserves... what step caused them to be lost??  rasterization??

dev_list = read_csv("/home/chrisgraywolf/shared/analysis/PA_loss/data/developing.csv")
df$developed = df$ISO3 %in% dev_list$`ISO-alpha3 code`
developed = df[df$developed,]
developed$Continent = "Developed"
developing = df[! df$developed,]
developing$Continent = "Developing"
df_dev = rbind(df, developed, developing)

df_small = melt(df_dev, measure.vars=c("loss","control_loss"))
df_small$variable = factor(df_small$variable,
	levels=rev(c("loss","control_loss")),
	labels=rev(c("Forest loss inside PA","Forest loss outside PA")))
df_small$Continent = factor(df_small$Continent, levels=unique(rev(c("Developed","Developing",sort(unique(df$Continent))))))
df_small$gp = paste(df_small$Continent,df_small$variable,df_small$cat_simple)

df_small = df_small[df_small$cat_simple != "Unknown",]
df_small$cat_simple = factor(df_small$cat_simple, levels=c("Nonstrict","Indigenous","Strict"))

p = ggplot(df_small[table(df_small$gp)[df_small$gp] >= 10,], aes(x=Continent, y=value*100, color=variable, shape=variable)) +
  stat_summary(fun.data = quants, geom = "errorbar",
               position=position_dodge(width=2/3),width=0) +
  stat_summary(fun.y=median, geom = "point",position=position_dodge(width=2/3)) +
  facet_wrap(~ cat_simple, nrow=2) +
  scale_color_manual(values=brewer.pal(3,"Set1")) +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(color="black"),
        axis.text.y=element_text(face=c(rep("plain",6),rep("bold",2))),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title.y=element_blank(),
        strip.text.y=element_text(angle=0),
        panel.grid.major.x=element_line(color="lightgray"),
        legend.position="bottom",
        legend.box="vertical",
        legend.background=element_rect(color="black")) +
  ylab("Annual forest loss rate (%)") +
  guides(color=guide_legend(title=NULL,reverse=TRUE),
         shape=guide_legend(title=NULL,reverse=TRUE))

png("output/loss_box_new.png", width=7, height=7, units="in", res=300)
p
dev.off()

######################################## simple model!

df_mod = na.omit(df[,c("loss","control_loss",
	"GDP","Continent","lat","long","cat_simple","GIS_AREA"),])
df_mod$abs_lat = abs(df_mod$lat); df_mod$log_GDP = log(df_mod$GDP); df_mod$log_GIS_AREA = log(df_mod$GIS_AREA)
sd_vec = apply(df_mod,2,sd)
GDP_range = range(df_mod$GDP)
df_mod[,sapply(df_mod, class) != "character"] = scale(df_mod[,sapply(df_mod, class) != "character"])[,]




start = Sys.time()
mod <- fitme(loss ~ abs_lat + control_loss + log_GDP + cat_simple + log_GIS_AREA + Continent + Matern(1 | long + lat),
		data=df_mod[1:2000,], family = "gaussian")
Sys.time() - start

# confint(mod,"Elevation")


n2 <- 20
n <- n1*n2
locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
covparms <- c(2,0.1,1/2,0)
y <- 7 + fast_Gp_sim(covparms, "matern_isotropic", locs)


X = model.matrix(~ abs_lat + control_loss + log_GDP + cat_simple + log_GIS_AREA + Continent, df_mod)
y = df_mod$loss
locs = cbind(df_mod$lat,df_mod$long)

fit = fit_model(y, locs, X, "matern_isotropic")
summary(fit)



start = Sys.time()
set.seed(2)
samps = foreach(i=1:4, .combine="rbind", .packages="spNNGP") %do% { # not dopar since we use n.omp.threads > 1
	mod_Bayes = spNNGP(loss ~ abs_lat + control_loss + log_GDP + cat_simple + log_GIS_AREA + Continent,
		data=df_mod,
		coords=cbind(df_mod$long,df_mod$lat),
		cov.model="exponential", # segfaults w/ "matern" and nu (unif.)
		starting=list("phi"=rexp(1), "sigma.sq"=rexp(1), "tau.sq"=rexp(1)),
		tuning=list("phi"=1, "sigma.sq"=1, "tau.sq"=1),
		priors=list("phi.Unif"=c(1e-3, 1e3), "sigma.sq.IG"=c(0.001, 0.001), "tau.sq.IG"=c(0.001, 0.001)),
		n.samples=3500,
		n.omp.threads = 2)
	out = data.frame(mod_Bayes$p.beta.samples)
	out$chain = i
	out$sample = 1:nrow(out)
	return(out)
}
saveRDS(samps, "output/samps_std.RDS")
Sys.time() - start # 43 minutes

###################################### check tuning params. (proposal vars. ok)

trace = melt(samps, id.vars=c("chain","sample"))
pdf("output/trace_std.pdf")
for(i in 1:4)
print(ggplot(trace[trace$chain == i,], aes(x=sample, y=value)) +
	geom_line() +
	facet_wrap(~ variable,scales="free")) +
	theme_bw()
for(i in 1:4)
print(ggplot(trace[trace$chain == i & trace$sample >= 1001,], aes(x=sample, y=value)) +
	geom_line() +
	facet_wrap(~ variable,scales="free")) +
	theme_bw()
dev.off()
samps_burn = samps[samps$sample > 1000,]
samps_list = as.mcmc.list(lapply(split(samps_burn, samps_burn$chain),as.mcmc))
gelman.diag(samps_list, multivariate=FALSE, transform=TRUE, autoburnin=FALSE)$psrf # looks good - all upper limits < 1.02
max(gelman.diag(samps_list, multivariate=FALSE, transform=TRUE, autoburnin=FALSE)$psrf[,2],na.rm=T)

summary(lm(loss ~ abs_lat + control_loss + log_GDP + log_GIS_AREA + cat_simple + Continent,data=df_mod[df_mod$cat_simple != "Indigenous",]))

###################################### paper text backtrans. summary table

summary_tab = t(sapply(colnames(samps_burn), function(colname){
	x = samps_burn[,colname]
	out = quantile(x,c(.025,.5,.975))
	sd_out = out * sd_vec["rel_loss"] / sd_vec[colname] # https://stats.stackexchange.com/questions/269452/formula-for-standardized-regression-coefficientsderivation-and-intuition
	names(sd_out) = paste0(names(sd_out),"_backtrans")
	out = c(out, sd_out, pval=mcmc.pval(x, testlim = 0, sided = 2, ptype="mcmc"))
	return(out)
}))
100*summary_tab*log(2)

# y = ... + beta * log(x)
# ... + beta * log(k * x)
# ... + beta * log(k) + beta * log(x)
# double x -> y increases additively by beta*log(2)
# mult. x by k -> y increases additively by beta*log(k)
# https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faqhow-do-i-interpret-a-regression-model-when-some-variables-are-log-transformed/

# compare lowest/highest GDP effects in model dataset:
df_mod_country = na.omit(df_main[df_main$space_sel,c("Country","rel_loss","rel_loss_5000","GDP","Continent","lat","long","cat_simple","GIS_AREA"),])
low = df_mod_country$Country[df_mod_country$GDP %in% min(df_mod_country$GDP,na.rm=T)][1]
high = df_mod_country$Country[df_mod_country$GDP %in% max(df_mod_country$GDP,na.rm=T)][1]
k = max(df_mod_country$GDP,na.rm=T)/min(df_mod_country$GDP,na.rm=T) # mult. increase from lowest -> highest GDP
(100*summary_tab*log(1/k))["log_GDP","50%_backtrans"]
(summary_tab*log(1/k))["log_GDP","50%_backtrans"] / sd(df_mod_country$rel_loss) # put in terms of response var. s.d. for context

# compare lowest/highest area effects in model dataset:
df_mod_area = na.omit(df_main[df_main$space_sel,c("Country","NAME","rel_loss","rel_loss_5000","GDP","Continent","lat","long","cat_simple","GIS_AREA"),])
low = df_mod_area$GIS_AREA[df_mod_area$GIS_AREA %in% min(df_mod_area$GIS_AREA,na.rm=T)][1] # 0.01 km^2 (by design)
df_mod_area[df_mod_area$GIS_AREA %in% min(df_mod_area$GIS_AREA,na.rm=T),c("Country","NAME","GIS_AREA")] # "VEP nr.129098" in Estonia
high = df_mod_area$GIS_AREA[df_mod_area$GIS_AREA %in% max(df_mod_area$GIS_AREA,na.rm=T)][1] # 80,570
df_mod_area[df_mod_area$GIS_AREA %in% max(df_mod_area$GIS_AREA,na.rm=T),c("Country","NAME","GIS_AREA")] # "Alto Rio Negro" in Brazil
k = max(df_mod_area$GIS_AREA,na.rm=T)/min(df_mod_area$GIS_AREA,na.rm=T) # mult. increase from lowest -> highest GDP
(100*summary_tab*log(1/k))["log_GIS_AREA","50%_backtrans"]
(summary_tab*log(1/k))["log_GIS_AREA","50%_backtrans"] / sd(df_mod_country$rel_loss) # put in terms of response var. s.d. for context

# v-- rel_loss_5000 was not log transformed
summary_tab["rel_loss_5000",] * 100 * 0.1 # multiply by 0.1 since an increase of 1 in loss_5000 is an additive increase of 100%..;
summary_tab["abs_lat",] * 100 * 10 # effect of 10 deg lat. increase

###################################### summary table

summary_tab = t(apply(samps_burn,2,function(x){
	out = c(quantile(x,c(.025,.5,.975)),pval=mcmc.pval(x, testlim = 0, sided = 2, ptype="mcmc"))
	return(out)
}))
summary_tab
write.csv(summary_tab, "output/summary_all.csv", row.names=FALSE)

samps = samps_burn
samps$cat_simpleIndigenous = samps$ContinentAfrica = 0

Continent = paste0("Continent",gsub(" ",".",unique(df_mod$Continent)))
cat_simple = paste0("cat_simple",gsub(" ",".",unique(df_mod$cat_simple)))
comparisons = rbind(expand.grid(cat_simple,cat_simple,stringsAsFactors=FALSE),
	expand.grid(Continent,Continent,stringsAsFactors=FALSE))
comparisons = comparisons[comparisons[,1] != comparisons[,2],]

comparisons = data.frame(comparisons,do.call("rbind",lapply(1:nrow(comparisons), function(i){
	diff = samps[,comparisons[i,1]] - samps[,comparisons[i,2]]
	out = data.frame(
		lower = quantile(diff, 0.025),
		median = median(diff),
		upper = quantile(diff, 0.975),
		p_raw=mcmc.pval(diff, testlim = 0, sided = 2, ptype="mcmc"),
		negative = mean(diff) < 0) # drop (after finding adj. pvals.) for simplifity
	return(out)
})))
comparisons$p_adj = p.adjust(comparisons$p_raw, method="fdr")
write.csv(comparisons[! comparisons$negative,], "output/comparisons_all.csv", row.names=FALSE)
comparisons[! comparisons$negative,]
comparisons[comparisons$p_adj < 0.05 & ! comparisons$negative,]






