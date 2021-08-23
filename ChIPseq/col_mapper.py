#!/usr/bin/python
import sys 

def load_table( source ):
	L = {}
	fin = open( source ) 
	while 1:
		line = fin.readline() 
		if not line:
			break
		qline = line.strip() 
		if len( qline ) > 0:
			H = qline.split( "\t" ) 
			#print H
			L[ H[0].split( "." )[0] ] = "\t".join( H[1:] ) 
	fin.close() 

	return L 

if len( sys.argv ) != 3:
	sys.stderr.write( "table-1, target-col\n" ) 
	sys.exit()
else:
	L = load_table( sys.argv[1] ) 

	while 1:
		line  = sys.stdin.readline() 
		if not line:
			break
		qline = line.strip() 
		if len( qline ) > 0:
			H = qline.split( "\t" ) 
			#print H			
			if L.has_key( H[ int( sys.argv[2] ) ].split( ".")[0] ):
				H.append( L[ H[ int( sys.argv[2] ) ].split(".")[0] ] )
			else:
				H.append( "" ) 

			sys.stdout.write( "\t".join( H ) + "\n" )


