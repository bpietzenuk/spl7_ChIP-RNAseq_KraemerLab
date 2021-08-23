options(scipen=999)
args = commandArgs(trailingOnly=TRUE)
library( ggplot2 )

#original data entry
d<-read.csv( file = args[1] )


thePlot = ggplot( data.frame( x = d$x, y = d$y ), aes( x = x, y = y ) ) 

frect = data.frame( x1 = c(0), x2 = c(2000), y1 = c(0), y2 = c(1) )


##### change the color of the block by setting fill !!!! 
thePlot = thePlot + geom_rect( mapping = aes( xmin = 0, xmax = 2000, ymin = 0, ymax = max( d$y )  ), fill = "lightgrey", alpha = 0.05  )
thePlot = thePlot + geom_line( size = 0.8 )
#print( data.frame( x = c( d$x, reverse( d$x ) ), y = c( d$ll, reverse( d$hh ) ) ) )

#set to FALSE if you dont want it 
plotPolyGon = FALSE 
polyGonFillColor = rgb( 0.9, 0.9, 1 ) 
polyGonLineColor = "skyblue4" 
#print( data.frame( x = c( d$x, rev( d$x ) ), y = c( d$ll, rev( d$hh ) ) ) )

if ( plotPolyGon == TRUE ) {
	thePlot = thePlot + geom_polygon( data = data.frame( x = c( d$x, rev( d$x ) ), y = c( d$ll, rev( d$hh ) ) ), aes( x = x, y = y ), fill = polyGonFillColor, alpha = 0.6 )
	thePlot = thePlot + geom_line( data = data.frame( x = d$x, y = d$mn ), aes( x = x, y = y ), col = polyGonLineColor, linetype = 2 )
}
#thePlot = thePlot + geom_line( size = 1.1 )

#axis labels 
thePlot = thePlot + labs( title = "Relative distance SPL7 ChIPseq peaks to TSS", y = "Frequency", x = "Distance nt")


######### THIS SECTION IS SENSITIVE AND SPeCIFIC TO THE DATA 
breaksLeft = seq( -3000,0, by = 500 )
breaksRight = seq( 2000, 3000, by = 500 )

ticksLeft = as.character( breaksLeft )  
ticksLeft[ length( ticksLeft ) ] = "TSS"
ticksRight = as.character( breaksRight - 2000 )
ticksRight[ 1 ] = "TTS"

breaks = c( breaksLeft, breaksRight )
ticks = c( ticksLeft, ticksRight )

thePlot = thePlot + scale_x_continuous( breaks = breaks, labels = ticks )
############################################################################

thePlot = thePlot + theme_classic(base_size = 30) 
thePlot = thePlot + theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(expand = c(0, 0))
plot(thePlot)
png(args[2] , height = 800, width = 1024)
plot( thePlot )
dev.off()
pdf( args[3], height = 9, width = 17, useDingbats=FALSE)
plot( thePlot )
dev.off()