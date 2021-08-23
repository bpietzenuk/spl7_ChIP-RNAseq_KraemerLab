import sys 

while 1:
    line = sys.stdin.readline() 
    if not line:
        break 
    qline = line.strip() 
    if len( qline ) > 0:
        H = qline.split( "\t" ) 

        S = int( H[1] ) 
        E = int( H[2] ) 

        M = ( S + E ) / 2 
        H[1] = str( M ) 
        H[2] = str( M ) 

        print "\t".join( H ) 
