#primlight specifically for mutation engine

def Prim( S, E ):
    """ build prim """
    return str( S ) + "^" + str( E )


def PrimS( prim ):
    """ return S """
    return int( prim.split("^")[0] )

def PrimE( prim ):
    """ return E """
    return int( prim.split("^")[1] )

def PrimC( prim ):
    """ return both """
    return PrimS( prim ), PrimE( prim )

def PrimL( prim ):
    S, E = PrimC( prim )
    L = 1 + abs( E - S )
    return L

def dirSwap( S, E ):
    if S > E:
        tmp = S
        S = E
        E = tmp
    return S, E


def PrimO( prim ):
    if PrimS( prim ) > PrimE( prim ):
        return "-"
    else:
        return "+"

def fusePrims( prim1, prim2):
    """ expected same O for fusion """

    return Prim( PrimS( prim1 ), PrimE( prim2 ) )


def pointInPrim( P, prim ):
    S, E = PrimC( prim )
    S, E = dirSwap( S, E )

    if P >= S and P <= E:
        return True
    else:
        return False

def primOverlap( prim1, prim2 ):
    S1, E1 = PrimC( prim1 )
    S2, E2 = PrimC( prim2 )

    if pointInPrim( S1 , prim2 ) == True:
        return True
    elif pointInPrim( E1, prim2 ) == True:
        return True
    elif pointInPrim( S2, prim1 ) == True:
        return True
    elif pointInPrim( E2, prim1 ) == True:
        return True
    else:
        return False

def primInPrimFull( prim1, prim2 ):
    S, E = PrimC( prim1 )

    if pointInPrim( S, prim2 ) == True and pointInPrim( E, prim2 ) == True:
        return True
    else:
        return False

def primOverlapSize( prim1, prim2 ):
    if primOverlap( prim1, prim2 ) == False:
        return 0

    S1, E1 = PrimC( prim1 )
    S1, E1 = dirSwap( S1, E1 )
    S2, E2 = PrimC( prim2 )
    S2, E2 = dirSwap( S2, E2 )

    return 1 + ( min( E1, E2 ) - max( S1, S2 ) )

def splitPrim( prim, junc ):
    jS, jE = PrimC( junc )

    lS = PrimS( prim )
    lE = jS

    rS = jE
    rE = PrimE( prim )

    return Prim( lS, lE ), Prim( rS, rE )


def masterPrim( prims ):
    return Prim( PrimS( prims[0] ), PrimE( prims[-1] ) )

def masterPrimF( prims ):
    U = [] 
    for prim in prims:
        U.append( PrimS( prim ) )
        U.append( PrimE( prim ) ) 

    return Prim( min( U ), max( U ) )

def intronToJunc( intron ):
    S, E = PrimC( intron )
    if E > S:
        S = S - 1
        E = E + 1
    else:
        S = S + 1
        E = E - 1
    return Prim( S, E )

def intronsToJuncs( introns ):
    juncs = []
    for intron in introns:
        juncs.append( intronToJunc( intron ) )
    return juncs

def juncToIntron( junc ): #new 3-12-2010
    S, E = PrimC( junc )
    if S <= E:
        S = S + 1
        E = E - 1
    else:
        S = S - 1
        E = E + 1
    return Prim( S, E )

def juncsToIntrons( juncs ): #new 3-12-2010
    introns = []
    for jun in juncs:
        intron.append( intronToJunc( junc ) )
    return introns


def fullInFilter( prims, MasterP ):
    rP = []
    #print prims
    for prim in prims:
        #print prim, MasterP
        if primInPrimFull( prim, MasterP ) == True:
            rP.append( prim )

    return rP

def _primToSplit( prims, junc ):
    u = -1
    #print prims
    for i in range( len( prims ) ):
        if primInPrimFull( junc, prims[i] ) == True:
            u = i
            break
    return u

def splitMasterPrim( Mprim, juncs ):
    #pre-filter
    err = 0
    prims = [ Mprim ]


    Fjuncs = fullInFilter( juncs, Mprim )
    #print len( Fjuncs ) 

    #I removed the (commented out these lines)
    
    #if len( Fjuncs ) != len( juncs ):
    #    return prims, -1 #indicated error
    
    for junc in Fjuncs:
        p = _primToSplit( prims, junc )

        if p > -1:
            lprim, rprim = splitPrim( prims[p], junc )
            #insertion

            prims[p] = lprim

            if p == len( prims ) - 1:
                prims.append( rprim )
            else:
                prims.insert( p + 1, rprim )
        else:
            err+=1 #conflict

    return prims, err #splitted 


def swapPrim( prim ):
    return Prim( PrimE( prim ), PrimS( prim ) )

def swapPrims( prims ):
    prims.reverse() #rever list

    rL = []

    for prim in prims:
        rL.append( swapPrim( prim ) )

    return rL 

def mxS( prims ):
    mx = 0
    for prim in prims:
        mx = max( mx, len( prim.split( "^" )[0] ) )
    return mx

def primPad( prims ):
    #print prims 
    mx = mxS( prims )

    for i in range( len( prims ) ):
        prims[i] = "0" * ( mx - len( prims[i].split( "^" )[0] ) ) + prims[i] 
    return prims

def primSort( prims ):
    padded = primPad( prims )
    padded.sort()
    S, E = PrimC( padded[0] )
    if S > E:
        padded.reverse()
    return padded
