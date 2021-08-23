import sys 
import primClusters 
import primLight_V3

L = {} 
sys.stderr.write( "staring\n" ) 
for f in sys.argv[1:]:
	sys.stderr.write( f + "\n" ) 
	fin = open( f ) 
	while 1:
		line = fin.readline() 
		if not line:
			break 
		qline = line.strip() 
		if len( qline ) > 0:
			H = qline.split() 

			if not L.has_key( H[0] ):
				L[ H[0] ] = [] 
			L[ H[0] ].append( primLight_V3.Prim(H[1],H[2] ) + "^" + H[3] )

	fin.close() 


for key in L.keys():
	ww = primClusters.clusterPrims_backwards( L[key] )

	for www in ww:
		print key + "\t" +  www  
 
