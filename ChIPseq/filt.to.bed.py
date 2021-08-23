import sys 

while 1:
	line = sys.stdin.readline() 
	if not line:
		break 
	qline = line.strip() 
	if len( qline ) > 0:
		H = qline.split( "^" ) 

		O = [] 
		O.append( H[0] ) 
		O.append( H[3] ) 
		O.append( H[4] ) 
		O.append( H[1] + "-" + H[2] ) 
		O.append( "1" ) 

		print "\t".join( O ) 
