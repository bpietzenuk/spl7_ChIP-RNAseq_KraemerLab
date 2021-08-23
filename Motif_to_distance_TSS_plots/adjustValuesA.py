import sys 

while 1:
	line = sys.stdin.readline() 
	if not line:
		break 
	qline = line.strip() 
	if len( qline ) > 0:
		H = qline.split() 
		#print H
		if H[3] == "ATG":
			H.append( str( -1 * int( H[2] ) ) )
		elif H[3] == "TERM":
			H.append( str( 2000 + int( H[2] ) ) )
		else:
			H.append( H[2] ) 

		print "\t".join( H ) 
