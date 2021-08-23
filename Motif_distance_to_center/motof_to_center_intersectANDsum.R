args = commandArgs(trailingOnly=TRUE)
library("ggplot2")
####Add table here which contains chromosome peak.start peak.end gene_id peak.name peak.to.gene.association (original data contained annotation, too)
###Data has headers
low<-read.delim( file = "minus.center.txt", header = TRUE, sep = "\t") ####change to read.table in case your data does not contain annotation
low$feature_minus <- paste0(low$V5, "_minus") 
#head(low)
control<-read.delim( file = "plus.center.txt", header = TRUE, sep = "\t") ####change to read.table in case your data does not contain annotation
control$feature_plus <- paste0(control$V5, "_plus") 
#head(control)

###Load motif data, merge with center information; merge by peak name of lowCu or controlCu peaks and calculate motif distance to center
##Be aware that you need to have determined the genomic location of the motif start and end first!
motif <- read.delim( file = args[1], row.names = NULL, header = FALSE, sep = "\t" )

motif_low <- merge(motif,low, by.x ="V9", by.y = "feature_minus", all.x = T)
#head(motif_low)
motif_low_control <- merge(motif_low,control, by.x ="V9", by.y = "feature_plus", all.x = T)
#head(motif_low_control)
motif_low_control = as.data.frame(t(apply(motif_low_control,1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )))
motif_low_control_cut <- motif_low_control[,c(1,2,8,25)]
head(motif_low_control_cut)

####distance calculation
center_plot = motif_low_control_cut
colnames(center_plot) <- c("peak_name","AGI", "motif_start", "center")
center_plot$motif_start = as.numeric(center_plot$motif_start)
center_plot$center = as.numeric(center_plot$center)
str(center_plot)
center_plot$distance = center_plot$motif_start - center_plot$center
head(center_plot)

write.table(center_plot,file = args[2], sep = "\t")

####Setup the plotting area######
#+ xlim(-500,500) 
###############################################################################
plot = ggplot(center_plot, aes(x=distance)) + geom_histogram(aes(y=..density..), binwidth = 70, colour="black", fill="grey") + geom_density(colour = "red") + theme_classic(base_size = 30) + scale_x_continuous(limits = c(-500,500), breaks = c(-500,-400,-300,-200,-100, 0,100,200,300,400, 500), labels=c("-500","-400","-300","-200","-100", "center","100","200","300","400","500"))
plot = plot + labs( title = "Motif position relative to peak center", y = "Density (10-3)", x = "Distance to peak center (nt)") + scale_y_continuous(expand = c(0, 0))
plot( plot )
png( args[3] , height = 800, width = 1024 )
plot( plot )
dev.off()
pdf( args[4], height = 9, width = 13, useDingbats=FALSE)
plot( plot )
dev.off()


#center_motif =merge(center_d, motif, by.x = "V6" , by.y = "V1", )
#head(center_motif)
#center_plot <- center_motif[c(1, 9, 15)]
#colnames(center_plot) <- c("AGI", "center", "motif_start")
#head(center_plot)
#center_plot$distance = center_plot$motif_start - center_plot$center
#head(center_plot)
