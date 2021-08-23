import sys 
import primLight_V3 

while 1:
    line = sys.stdin.readline() 
    if not line:
        break 
    qline = line.strip() 
    if len( qline ) > 0:
        H = qline.split()[-1].split( "," ) 

        for i in range( len( H ) ):
            for j in range( i + 1, len( H ) ):
                u = H[i].split( "^" )[-1] 
                u2 = H[j].split( "^" )[-1] 

                if u != u2:
                    l1 = primLight_V3.PrimL( H[i] ) 
                    l2 = primLight_V3.PrimL( H[j] ) 

                    s = primLight_V3.primOverlapSize( H[i], H[j] ) 

                    print float( s ) / float( min( l1, l2 ) )
