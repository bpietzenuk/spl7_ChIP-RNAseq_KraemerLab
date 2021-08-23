import sys 
sys.path.append( "/netscratch/dep_coupland/grp_coupland/bioinformatics/edouard/tools" ) 
import newMiscFunctions

if len( sys.argv ) != 2:
	sys.stderr.write( "Usage FastaIn < regions > output\n" ) 
	sys.exit() 

sys.stderr.write( "reading fasta\n" ) 
L = newMiscFunctions.readFastaFull( sys.argv[1] ) 
sys.stderr.write( "done reading\n" ) 
sys.stderr.write( "start scanning \n" ) 
while 1:
	line = sys.stdin.readline() 
	if not line:
		break 
	qline = line.strip() 
	if len( qline ) > 0:
		H = qline.split( "\t" ) 

		if len( H ) == 1:
			H = qline.split( "," )

		S = int( H[3] ) - 1 #+ int( sys.argv[4] ) 
		E = int( H[4] ) - 1 #+ int( sys.argv[4] ) 

		S, E = min( S, E ), max( S, E ) 

		S = S #- int( sys.argv[2] ) 
		E = E #+ int( sys.argv[3] )
 
		print ">" + H[0] 
		print newMiscFunctions.getSegment( H[1], S, E, L, strand = H[4] ) 

sys.stderr.write( "Done scanning\n" ) 
