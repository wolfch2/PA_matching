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

project_dir = "/home/chrisgraywolf/analysis_desktop/PA_matching/"
setwd(project_dir)

########################################

new = read.csv("data_processed/data_matched.csv", as.is=TRUE)
PA_table = read.csv("data_processed/PAs_tab.csv")
new$PAs = PA_table$WDPAID[match(new$PAs,PA_table$ID)]

######################################## join, compare etc. w/ old

df = readRDS("/home/chrisgraywolf/shared/analysis/PA_loss//data/WDPA_post_v2.RDS")
df = df[df$space_sel,]

df$inside = 1-(1-df[df$space_sel,]$rel_loss)^(1/16)
df$outside = 1-(1-df[df$space_sel,]$rel_loss_5000)^(1/16)

df$control_loss = new$control_loss[match(df$WDPAID, new$PAs)]
df$loss = new$loss[match(df$WDPAID, new$PAs)]

r1 = round(cor(na.omit(cbind(df$outside, df$control_loss)), method="spearman")[1,2],2)
p1 = ggplot(df, aes(x=100*outside, y=100*control_loss)) +
        geom_point(size=0.5, alpha=0.5) +
        scale_x_log10(labels=function(x) formatC(x, format="fg")) +
        scale_y_log10(labels=function(x) formatC(x, format="fg")) +
        annotation_logticks() +
        theme_bw() +
        theme(axis.ticks=element_line(color="black"),
              axis.text=element_text(color="black"),
              plot.title=element_text(hjust=0.5),
              panel.border=element_rect(color="black")) +
        #geom_density_2d() +
        geom_abline(a=0,b=1,color="red",size=1.5) +
        stat_smooth(method="lm") +
        xlab("% Forest loss outside (old)") +
        ylab("% Forest loss outside (new)") +
        ggtitle(paste0("Control vs. Buffer forest loss rates (rho = ",r1,")"))

r2 = round(cor(na.omit(cbind(df$inside, df$loss)), method="spearman")[1,2],2)
p2 = ggplot(df, aes(x=100*inside, y=100*loss)) +
        geom_point(size=0.5, alpha=0.5) +
        scale_x_log10(labels=function(x) formatC(x, format="fg")) +
        scale_y_log10(labels=function(x) formatC(x, format="fg")) +
        annotation_logticks() +
        theme_bw() +
        theme(axis.ticks=element_line(color="black"),
              axis.text=element_text(color="black"),
              plot.title=element_text(hjust=0.5),
              panel.border=element_rect(color="black")) +
        #geom_density_2d() +
        geom_abline(a=0,b=1,color="red",size=1.5) +
        stat_smooth(method="lm") +
        xlab("% Forest loss inside (old)") +
        ylab("% Forest loss inside (new)") +
        ggtitle(paste0("Updated vs. old PA forest loss rates (rho = ",r2,")"))

r3 = round(cor(na.omit(cbind(df$inside, df$outside)), method="spearman")[1,2],2)
p3 = ggplot(df, aes(x=100*inside, y=100*outside)) +
        geom_point(size=0.5, alpha=0.5) +
        scale_x_log10(labels=function(x) formatC(x, format="fg")) +
        scale_y_log10(labels=function(x) formatC(x, format="fg")) +
        annotation_logticks() +
        theme_bw() +
        theme(axis.ticks=element_line(color="black"),
              axis.text=element_text(color="black"),
              plot.title=element_text(hjust=0.5),
              panel.border=element_rect(color="black")) +
        #geom_density_2d() +
        geom_abline(a=0,b=1,color="red",size=1.5) +
        stat_smooth(method="lm") +
        xlab("% Forest loss inside (old)") +
        ylab("% Forest loss outside (old)") +
        ggtitle(paste0("Old forest loss rates (rho = ",r3,")"))

r4 = round(cor(na.omit(cbind(df$loss, df$control_loss)), method="spearman")[1,2],2)
p4 = ggplot(df, aes(x=100*loss, y=100*control_loss)) +
        geom_point(size=0.5, alpha=0.5) +
        scale_x_log10(labels=function(x) formatC(x, format="fg")) +
        scale_y_log10(labels=function(x) formatC(x, format="fg")) +
        annotation_logticks() +
        theme_bw() +
        theme(axis.ticks=element_line(color="black"),
              axis.text=element_text(color="black"),
              plot.title=element_text(hjust=0.5),
              panel.border=element_rect(color="black")) +
        #geom_density_2d() +
        geom_abline(a=0,b=1,color="red",size=1.5) +
        stat_smooth(method="lm") +
        xlab("% Forest loss inside (new)") +
        ylab("% Forest loss outside (new)") +
        ggtitle(paste0("New forest loss rates (rho = ",r4,")"))

png("output/old_new_loss_comp.png", width=12, height=12, units="in", res=300)
plot_grid(p1,p2,p3,p4,nrow=2,labels=LETTERS[1:4])
dev.off()

