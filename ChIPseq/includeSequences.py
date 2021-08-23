import sys 

L = {} 
fin = open( sys.argv[1] ) 
while 1:
	line = fin.readline() 
	if not line:
		break 
	qline = line.strip() 
	if len( qline ) > 0:
		H = qline.split() 

		L[ H[0] ] = True 
fin.close() 

cpp = False 
while 1:
	line = sys.stdin.readline() 
	if not line:
		break 
	qline = line.strip() 
	if len( qline ) > 0:
		if qline[0] == ">":
			if L.has_key( qline[1:] ):
				cpp = True 
			else:
				cpp = False 

		if cpp:
			print qline 
