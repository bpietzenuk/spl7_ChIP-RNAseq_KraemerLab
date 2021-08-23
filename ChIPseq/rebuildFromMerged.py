import sys 

while 1:
    line = sys.stdin.readline() 
    if not line:
        break 
    qline = line.strip() 
    if len( qline ) > 0:
        H = qline.split( "^" ) 

        P = [] 
        P.append( H[0] ) 
        P.append( H[3] ) 
        P.append( H[4] ) 
        P.append( H[1] + "-" + H[2] ) 
        P.append( "2" ) 
        
        print "\t".join( P ) 
