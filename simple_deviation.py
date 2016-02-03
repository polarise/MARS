#!/usr/bin/env python
from __future__ import division
import sys
import scipy

try:
	f1 = sys.argv[1]
	f2 = sys.argv[2]
#	fdr = float( sys.argv[3] )
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <f1> <f2>"
	sys.exit( 0 )
	
#try:
#	assert 0 < fdr <= 1
#except:
#	raise ValueError( "Invalid value for 'fdr': it should be in (0,1]" )
#	sys.exit( 0 )

for fdr in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]:
	f1_data = dict()
	for row in open( f1 ):
		L = row.strip().split( '\t' )
		if float( L[-2] ) <= fdr:
			if L[-1] == 'Inf' or L[-1] == '-Inf': continue
			f1_data[L[0]] = float( L[-1] )

	dev = 0
	f2_data = dict()
	for row in open( f2 ):
		L = row.strip().split( '\t' )
		if L[-1] == 'Inf' or L[-1] == '-Inf': continue
		try:
			comp = f1_data[L[0]]
		except KeyError:
			continue
		dev += ( float( L[-1] ) - comp )**2
	
	print "%s\t%s" % ( fdr, scipy.sqrt( dev ))

