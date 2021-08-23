import sys 

while 1:
	line = sys.stdin.readline() 
	if not line:
		break 
	qline = line.strip()
	if len( qline ) > 0:
		H = qline.split() 

		W = H[0].split( "^" ) 
		
		P = [ W[1] + "-" + W[2] ]
		P.append( H[-1] ) 

		print "\t".join( P ) 
