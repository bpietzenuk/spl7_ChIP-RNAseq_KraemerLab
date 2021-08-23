import primLight_V3

def mergePrims( p1, p2 ):
    S1, E1 = primLight_V3.PrimC( p1 )
    S2, E2 = primLight_V3.PrimC( p2 )

    if primLight_V3.PrimO( p1 ) == "+":
        return primLight_V3.Prim( min( S1, S2 ), max( E1, E2 ) )
    else:
        return primLight_V3.Prim( max( S1, S2 ), min( E1, E2 ) )
    

def cleanClusters( clusters ):
    for i in range( len( clusters ) ):
        H = clusters[i].split( "," )
        R = [] 
        for j in range( 1, len( H ) ):
            R.append( H[j].lstrip( "0" ) )
        clusters[i] = ",".join( R )

    return clusters

def clusterPrims_backwards( primList ):
    if len( primList ) < 1:
        return []

    

    #sort primList
    primList = primLight_V3.primSort( primList )

    clusters = []
    for prim in primList:
        #print prim
        clusters.append( "^".join( prim.split( "^" )[0:2] ) + "," + prim )

    if len( clusters ) < 2:
        return cleanClusters( clusters )  


    changes = 1
    while changes > 0:
        H = len( clusters ) - 1
        changes = 0
        while H > 0:
            #print clusters[H]
            cLine = clusters[H].split( "," )[0]
            pLine = clusters[H-1].split( "," )[0]

            S,E = primLight_V3.PrimC( cLine )
            S1, E1 = primLight_V3.PrimC( pLine )
            
            Q = False
            if abs( E - S1 ) == 1:
                Q = True
            elif abs( E1 - S ) == 1:
                Q = True 
            if primLight_V3.primOverlap( cLine, pLine ) == True or Q:
                merged = mergePrims( cLine, pLine )
                w = clusters[H-1].split( "," )
                w[0] = merged
                clusters[H-1] = ",".join( w )

                w2 = clusters[ H ].split( "," )
                clusters[ H - 1 ] = clusters[ H-1 ] + "," + ",".join( w2[1:] )
                del clusters[H]
                changes = changes + 1
            H = H - 1

        
    return  cleanClusters( clusters )


def clusterPrims_backwards_restrain( primList, minL ):
    if len( primList ) < 1:
        return []

    

    #sort primList
    primList = primLight_V3.primSort( primList )

    clusters = []
    for prim in primList:
        #print prim
        clusters.append( "^".join( prim.split( "^" )[0:2] ) + "," + prim )

    if len( clusters ) < 2:
        return cleanClusters( clusters )  


    changes = 1
    while changes > 0:
        H = len( clusters ) - 1
        changes = 0
        while H > 0:
            #print clusters[H]
            cLine = clusters[H].split( "," )[0]
            pLine = clusters[H-1].split( "," )[0]

            S,E = primLight_V3.PrimC( cLine )
            S1, E1 = primLight_V3.PrimC( pLine )
            
            Q = False
            if abs( E - S1 ) == 1:
                Q = True
            elif abs( E1 - S ) == 1:
                Q = True 
            if primLight_V3.primOverlapSize( cLine, pLine ) >= minL or Q:
                merged = mergePrims( cLine, pLine )
                w = clusters[H-1].split( "," )
                w[0] = merged
                clusters[H-1] = ",".join( w )

                w2 = clusters[ H ].split( "," )
                clusters[ H - 1 ] = clusters[ H-1 ] + "," + ",".join( w2[1:] )
                del clusters[H]
                changes = changes + 1
            H = H - 1

        
    return  clusters 



