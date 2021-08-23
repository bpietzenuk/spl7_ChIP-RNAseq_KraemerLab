import sys 

L = {} 
fin = open( sys.argv[1] ) 
#fin.readline()
while 1:
	line = fin.readline() 
	if not line:
		break 
	qline = line.strip() 
	if len( qline ) > 0:
		L[ qline ] = True 
fin.close() 

sys.stdin.readline()
while 1:
	line = sys.stdin.readline() 
	if not line:
		break 
	qline = line.strip() 
	if len( qline ) > 0:
		H = qline.split() 
		#print H
		if L.has_key( H[5] ):
			print H[4] 
