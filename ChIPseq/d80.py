import sys 
import primLight_V3

def ov( p1, p2 ):
	l = float( min( primLight_V3.PrimL( p1 ), primLight_V3.PrimL( p2 ) ) )

	osize = float( primLight_V3.primOverlapSize( p1, p2 ) )

	if ( ( osize / l ) >= 0.8 ):
		return True 

	return False 

def fov( r1, r2 ):
	for rr1 in r1:
		for rr2 in r2:
			if ov( rr1, rr2 ):
				return True 

	return False 	


L = {} 
while 1:
	line = sys.stdin.readline() 
	if not line:
		break 
	qline = line.strip()
	if len( qline ) > 0:
		H = qline.split( "\t" )
		if not L.has_key( H[0] ):
			L[ H[0] ] = []  
		
		L[ H[0] ]+= H[-1].split( "," ) 


for key in L.keys():
	sys.stderr.write( "new chromosome: %s\n" % ( key ) )
	R = []
	U = L[ key ] 
	U = primLight_V3.primSort( U )  
	for ll in U:
		#F = -1 
		#for i in range( len( R ) ):
		#	pass
		#if F == -1:
		#	R+= [ ll ] 
		R.append( [ ll ] ) 
	
	changes = True 
	while changes:
		sys.stderr.write( "new round\n" )
		changes = False
		R2 = [] 
		for r1 in R:
			F = -1
			for i in range( len( R2 ) ):
				if fov( r1, R2[ i ] ):
					F = i 
					break 

			if F == -1:
				R2.append( r1 ) 
			else:
				changes = True 
				R2[ F ]+= r1 
				sys.stderr.write( "merged:\n" )

		R = R2 

	for rr in R:
		print key + "\t" + ",".join( rr )

