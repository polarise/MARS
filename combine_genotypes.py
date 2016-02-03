#!/usr/bin/env python
"""
combines encoded genotypes across populations
"""
from __future__ import division
import sys

data1_D = dict()
f = open( sys.argv[1] )
for row_s in f:
	L = row_s.strip().split( '\t' )
	data1_D[L[0]] = L[1:]
f.close()

f = open( sys.argv[2] )
for row_s in f:
	L = row_s.strip().split( '\t' )
	if L[0] in data1_D:
		l_L = data1_D[L[0]]
		print "\t".join( [ L[0], l_L[0], l_L[1], str( int( l_L[2] ) + int( L[3] )), l_L[3]+","+L[4], l_L[4]+","+L[5] ] )
	else:
		print >> sys.stderr, "Missing SNP %s between populations." % L[0]
f.close()
