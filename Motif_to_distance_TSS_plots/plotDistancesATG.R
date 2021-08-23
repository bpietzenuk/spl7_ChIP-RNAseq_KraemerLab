args = commandArgs(trailingOnly=TRUE)
d<-read.table( args[1] )

getLow <- function( X ) { 
  s = X[ order( X ) ]
  s[ round( 0.025 * length( s ) ) ]
} 

getHi <- function( X ) { 
  s = X[ order( X ) ]
  s[ round( 0.975 * length( s ) ) ]
} 


rn = seq( -3000, 3000, by = 100 )
hst = hist( d$V11, breaks = rn, plot = F )

x = hst$mids 
yo = hst$density * diff( hst$breaks )

#png( "ATGrelDist.png" )
#plot( x, yo, type = "l", xlab = "Distance to ATG", ylab = "Frequency" )
#abline( v = 0, col = "red" )
#abline( v = 2000, col = "red" )
#dev.off() 

#write.csv( data.frame( x = x, y = y ), file = "ATGrelDist.csv" )
ss = 1000
pbaseMatrix = matrix( nrow = length( x ), ncol = ss ) 

q = 0
for ( i in 1: ss ) {
    waka = as.integer( system2( "wc", args = c( "-l", paste( args[2],"/dd_", i, sep = "" ), "| awk '{ print $1}'" ) , stdout = TRUE ) )
    print( paste( ":", waka, sep = " " ) )
    if ( waka > 100 ) {  
    d<-read.table( paste( args[2],"/dd_", i, sep = "" ) )
    hst = hist( d$V11, breaks = rn, plot = F ) 
    y = hst$density * diff( hst$breaks ) 
    q=q + 1
    pbaseMatrix[ , q ] = y
    } else { print( "skipping" ) }
}
baseMatrix = matrix( nrow = length( x ), ncol = q ) 

for ( i in 1:q ) {
  baseMatrix[ , i ] = pbaseMatrix[ ,i ] 
} 

#baseMatrix = baseMatrix[ !is.na( baseMatrix ), ] 
print( "here" )

mn = apply( baseMatrix,1,mean ) 
print( "here2" ) 
ll = apply( baseMatrix,1,getLow ) 
print( "here3" )
hh = apply( baseMatrix,1,getHi ) 
print( "here4" )
#lines( x, mn, col = "red" ) 
#lines( x, ll, col = "blue" ) 
#lines( x, hh, col = "blue" ) 
#dev.off()
print( x ) 
print( yo )
print( mn ) 
print( ll ) 
print( hh )
write.csv( data.frame( x = x, y = yo, mn = mn, ll = ll, hh = hh ), file = args[3] ) 

