args = commandArgs(trailingOnly=TRUE)
library("ggplot2")
####Add table here which contains chromosome peak.start peak.end gene_id peak.name peak.to.gene.association (original data contained annotation, too)
###Data has headers
d<-read.delim( file = args[1], row.names = NULL, header = FALSE, sep = "\t") ####change to read.table in case your data does not contain annotation

####calculate peak center
#head(d)
d$center <- d$V2 + (d$V4/2)
d$center <- as.integer(format(round(d$center, 0), nsmall = 0))
center_d <- d[c(1:8, 14)]
#head(center_d)
write.table(center_d,file = args[2], sep = "\t") ###saves center data under second argument name

###Load motif data, merge with center information and calculate motif distance to center
##Be aware that you need to have determined the genomic location of the motif start and end first!
motif <- read.delim( file = args[3], row.names = NULL, header = FALSE, sep = "\t" )
#head(motif)
center_motif =merge(center_d, motif, by.x = "V6" , by.y = "V1", )
#head(center_motif)
center_plot <- center_motif[c(1, 9, 15)]
colnames(center_plot) <- c("AGI", "center", "motif_start")
#head(center_plot)
center_plot$distance = center_plot$motif_start - center_plot$center
#head(center_plot)
write.table(center_plot,file = "center_d.txt", sep = "\t")
####Setup the plotting area######


#+ xlim(-500,500) 
###############################################################################
#plot = ggplot(center_plot, aes(x=distance)) + geom_density() + theme_classic() + scale_x_continuous(limits = c(-500,500), breaks = c(-500, 0, 500), labels=c("-500", "center", "500"))
plot = ggplot(center_plot, aes(x=distance)) + geom_histogram(aes(y=..density..), binwidth = 70, colour="black", fill="grey") + geom_density(colour = "red") + theme_classic(base_size = 30) + scale_x_continuous(limits = c(-500,500), breaks = c(-500,-400,-300,-200,-100, 0,100,200,300,400, 500), labels=c("-500","-400","-300","-200","-100", "center","100","200","300","400","500"))
plot = plot + labs( title = "Motif position relative to peak center", y = "Density (10-3)", x = "Distance to peak center (nt)") + scale_y_continuous(expand = c(0, 0))
plot( plot )
png( args[4] , height = 800, width = 1024 )
plot( plot )
dev.off()
pdf( args[5], height = 9, width = 13, useDingbats=FALSE)
plot( plot )
dev.off()