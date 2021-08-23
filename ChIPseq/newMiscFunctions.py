import sys 

def readFastaFull( source ):
	L = {} 
	header = ""
	sequence = "" 

	fin = open( source ) 
	while 1:
		line = fin.readline() 
		if not line:
			break 
		qline = line.strip() 
		if len( qline ) > 0:
			if qline[0] == ">":
				if header != "":
					L[ header ] = sequence 
				header = qline.replace( ">", "" ).split()[0] 
				sequence = "" 
			elif header != "":
				sequence+=qline.upper() 
	fin.close() 

	if header != "":
		L[ header ] = sequence 

	return L 

def complement( ch ):
	if ch == "A":
		ch = "T" 
	elif ch == "C":
		ch = "G" 
	elif ch == "G":
		ch = "C" 
	elif ch == "T":
		ch = "A" 

	return ch 

def revComplement( ss ):
	l = [ complement( X ) for X in list( ss ) ]
	l.reverse() 
	return "".join( l ) 


def _forceRange( P, S, E ):
	S, E = min( S, E ), max( S, E ) 
	return min( max( P, S ), E )  

def getSegment( ch, S, E, data, strand = "+" ):
	

	S, E = min( S, E ), max( S, E ) 
	
	if not data.has_key( ch ):
		return "" 

	S = _forceRange( S, 0, len( data[ ch ] ) - 1 ) 
	E = _forceRange( E, 0, len( data[ ch ] ) - 1 ) 

	ss = data[ ch ][ S: E + 1 ] 

	if strand == "+":
		return ss 
	else:
		return revComplement( ss ) 

def nextSequence( fin ):
	header = ""
	sequence = "" 

	while 1:
		otell = fin.tell() 
		line = fin.readline() 
		if not line:
			break 
		qline = line.strip() 
		if len( qline ) > 0:
			if qline[0] == ">":
				if header != "":
					fin.seek( otell ) 
					break 
				header = qline.replace( ">", "" ).split()[0] 

				sequence = "" 
			elif header != "":
				sequence+=qline 

	return header, sequence 


