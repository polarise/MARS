#!/usr/bin/env python
from __future__ import division
import sys
import random

# Script to collate sample names between both platforms

"""
Here's what this script does:
0. Assume the name of a sample is given by the prefix of the experiment ID
1. For each row in the details file
1.1 Split the sample name and take all except the last part as a key for dict()
1.2 If it is TrueSeq (RNA-Seq) experiment record as such; similar for Affymetrix Expression
1.3 Also record the tissue (as part of the key)
2. Open an appendable file for each tissue
2.1 If both platforms are present write the names in the tissue file; pick at random if multiple
Done!
"""

f = open( sys.argv[1] ) # The file beginning with GTEx_Analysis*
sample_names = dict()
for row in f:
	L = row.strip().split( '\t' )
	if L[0] == 'SAMPID': continue
	
	sample_name = L[0].split( '-' )
	key = ( "-".join( sample_name[:-1] ), L[6] )
	
	if key not in sample_names:
		if L[13] == 'TrueSeq.v1':
			sample_names[ key ] = { 'rnaseq': [ L[0] ], 'array': [] }
		elif L[13] == 'Affymetrix Expression':
			sample_names[ key ] = { 'array': [ L[0] ], 'rnaseq': [] }
	else:
		if L[13] == 'TrueSeq.v1':
			sample_names[ key ][ 'rnaseq' ] += [ L[0] ]
		elif L[13] == 'Affymetrix Expression':
			sample_names[ key ][ 'array' ] += [ L[0] ]
			
f.close()

total = 0
for s in sample_names:
	ofn = "_".join( s[1].split( ' ' )) + ".txt"
	f = open( "sample_names_by_tissue/" + ofn, 'a' )		
	if len( sample_names[ s ][ 'rnaseq' ] ) > 0 and len( sample_names[ s ][ 'array' ] ) > 0:
		print >> f, random.choice( sample_names[ s ][ 'rnaseq' ] ) + "\t" + random.choice( sample_names[ s ][ 'array' ] )
	f.close()
	
#	print "Root:       %s" % s[0]
#	print "Tissue:     %s" % s[1]
#	
#	total += len( sample_names[ s ][ 'rnaseq' ] ) + len( sample_names[ s ][ 'array' ] )
#	if len( sample_names[ s ][ 'rnaseq' ] ) == 0:
#		print "RNA-Seq:    None" 
#	else:
#		print "RNA-Seq:    %s" % ",".join( sample_names[ s ][ 'rnaseq' ] )
#	
#	if len( sample_names[ s ][ 'array' ] ) == 0:
#		print "Microarray: None"		
#	else:
#		print "Microarray: %s" % ",".join( sample_names[ s ][ 'array' ] )
#	
#	print
print >> sys.stderr, total
