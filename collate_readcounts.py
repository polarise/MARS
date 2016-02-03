#!/usr/bin/env python
from __future__ import division
import sys

names = [ row.strip().split( '\t' )[0] for row in open( sys.argv[1] ) ]

#fns = "/home/paulk/MARS/1A_raw_readcounts/%s.count"
fns = "/home/paulk/MARS/9A_additional_analyses/1_processed_data/%s.count"

data = dict()
for yri in names:
	f = open( fns % yri )
	c = 0
	for row in f:
		if c > 10: break
		L = row.strip().split( '\t' )
		if L[0] not in data:
			data[L[0]] = [L[1]]
		else:
			data[L[0]] += [L[1]]
		c += 0
	f.close()
	
for d in data:
	print "\t".join( [ d ] + data[d] )

