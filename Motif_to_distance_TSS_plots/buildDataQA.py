#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 10:14:25 2018

@author: edouardsevering
"""

import sys 
import re 
import primLight_V3


def primExtent( prim, S, E ):
    SS, EE = primLight_V3.PrimC( prim )
    
    if SS > EE:
        SS = SS + S
        EE = EE - E 
    else:
        SS = SS - S 
        EE = EE + E
        
    return primLight_V3.Prim( SS, EE )


def loadCDSdata( source, S, E ):
    data = {} 
    fin = open( source )
    while 1:
        line = fin.readline() 
        if not line:
            break 
        qline = line.strip() 
        if len( qline ) > 0:
            if qline[0] == "#":
                continue
            
            H = qline.split( "\t" )
            if H[2] == "gene":
                np = re.search( "ID=(.+?);", H[-1] + ";" )
                if np:
                    wo =  np.groups(1)[0].split( "," )[0]
                    #print wo
                    if not data.has_key( H[0] ):
                        data[ H[0] ] = {} 
                    if not data[ H[0] ].has_key( wo ):
                        data[ H[0] ][ wo ] = []
                            
                    if H[6] == "+":
                        data[ H[0] ][ wo ].append( primLight_V3.Prim( H[3], H[4] ) )
                    else:
                        data[ H[0] ][ wo ].append( primLight_V3.Prim( H[4], H[3] ) )
                            
    fin.close() 
    
    fStructure = {}
    for chrom in data.keys():
        fStructure[ chrom ] = []
        for nm in data[ chrom ].keys():
            masterPrim = primLight_V3.masterPrim( primLight_V3.primSort( data[ chrom ][nm] ) )
            masterPrimE = primExtent( masterPrim, S, E ) + "^" + nm 
            
            fStructure[ chrom ].append( masterPrimE )
    
    return fStructure 

def firstFullInSide( structure, ch, S, E ):
    M = ( int( S ) + int( E ) ) / 2 
    
    r = []
    for item in structure[ ch ]:
        if primLight_V3.pointInPrim( M, item ) == True:
            r.append( item ) 
    
    return r


def determineDistance( item, S, E ):
    M = ( int( S ) + int( E ) ) / 2
    
    vl = None 
    if primLight_V3.PrimO( item ) == "+":
        if primLight_V3.PrimS( item ) > M:
            return primLight_V3.PrimS( item ) - M, "ATG"
        elif primLight_V3.PrimE( item ) < M:
            return M - primLight_V3.PrimE( item ), "TERM"
        else:
            return None, None 
    
    if primLight_V3.PrimS( item ) < M:
        return M - primLight_V3.PrimS( item ), "ATG"
    elif primLight_V3.PrimE( item ) > M:
        return primLight_V3.PrimE( item ) - M, "TERM"
    else:
        return None, None 
    
def closestItem( structure, ch, S, E, mxS, mxE ):
    #M = ( int( S ) + int( E ) ) / 2 
    
    #dst = None
    #loc = None 
    #itm = None 
    r = [ ]
    for item in structure[ ch ]:
        tDst, tLoc = determineDistance( item, S, E )
        if tDst != None:
            if tLoc == "TERM" and tDst <= mxE:
                #if dst == None or dst > tDst:
                    #dst = tDst
                    #loc = tLoc
                    #itm = item
                r.append( [ tDst, tLoc, item ] )
            else:
                if tLoc == "ATG" and tDst <= mxS:
                    #if dst == None or dst > tDst:
                    r.append( [ tDst, tLoc, item ] )

    
    if len( r ) > 0:
        for i in range( len( r ) ):
            if i == 0:
                qq = r[ 0 ] 
                mn = r[ 0 ][ 0 ]
            else:
                if mn > r[ i ][ 0 ]: 
                    mn = r[i][ 0 ] 
                    qq = r[ i ] 

        #print "here", qq
        return [ qq ] 

    return r
    

        
myStructure = loadCDSdata( sys.argv[1], 0, 0 )

sl = int( sys.argv[2] ) 
el = int( sys.argv[3] )

while 1:
    line = sys.stdin.readline() 
    if not line:
        break 
    qline = line.strip() 
    if len( qline ) > 0:
        H = qline.split() 
        M = ( int( H[1] ) + int( H[2 ]) ) / 2
        seen = {}
        insides = firstFullInSide( myStructure, H[0], H[1], H[2] )
        for inside in insides:
            seen[ inside.split( "^" )[-1] ] = True
            #print inside.split( "^" )
            if primLight_V3.PrimO( inside ) == "+":
                D = M - primLight_V3.PrimS( inside )
            else:
                D = primLight_V3.PrimS( inside ) - M
                
            D2 = 2000 * float( D ) / primLight_V3.PrimL( inside )
            print "inside\t" + inside + "\t" + str( int( D2 ) ) + "\tINS\t" + str( M ) + "\t" + qline 
#"""i
            break #only one inside 

        if len( insides ) > 0:
            continue 
        #else:
        r = closestItem( myStructure, H[0],H[1],H[2], sl, el )
        for rr in r:#if dst != None:
                #print rr
                if not seen.has_key( rr[-1].split( "^" )[-1] ):
            
                    print "range\t" + rr[-1] + "\t" + str( rr[0] ) + "\t" + rr[1] + "\t" + str( M ) + "\t" + qline
#"""
