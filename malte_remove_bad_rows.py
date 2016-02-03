#!/usr/bin/env python
from __future__ import division
import sys

if __name__ == "__main__":
	try:
		fn = sys.argv[1] # TisMix_malte_OOB.txt
	except IndexError:
		print >> sys.stderr, "Usage: ./script.py <TisMix_malte_OOB.txt>"
		sys.exit( 1 )
	
	f = open( fn )
	for row in f:
		if row[0] == 'g':
			print row.strip()
			continue
		L = row.strip().split( "\t" )
		if 'NA' in L[1:]: continue
		if 0 in map( float, L[1:] ): continue
		print "\t".join( L )
	f.close()
