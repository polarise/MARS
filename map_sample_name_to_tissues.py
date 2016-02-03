#!/usr/bin/env python
from __future__ import division
import sys
import re

unwatned = [ 'S', 'h' ]

try:
	ann_fn = sys.argv[1]
	samples_fn = sys.argv[2]
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <annotation_fn> <samples_fn>"
	sys.exit( 0 )

samples = dict()
f = open( samples_fn )
samples = [ row.strip().split( '\t' )[1] for row in f if row[0] not in unwatned and row[0] == '*' ]
f.close()

# columns 6 and 7 are the ones of interest
f = open( ann_fn )
sample2tissue = dict()
for row in f:
	if row[0] in unwatned: continue
	L = row.strip().split( '\t' )
	for s in samples:
		if re.search( L[0], s ):
			sample2tissue[s] = L[5], L[6]
f.close()

for s in samples:
	print "\t".join( [ s, sample2tissue[s][0], sample2tissue[s][1] ] )
	


