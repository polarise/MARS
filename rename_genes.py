#!/usr/bin/env python
from __future__ import division
import sys

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		sys.exit( 1 )
	
	f = open( fn )
	for row in f:
		if row[0] in [ '#' ]: continue
		if row[0] == 'p': print row.strip(); continue
		L = row.strip().split( '\t' )
		print "\t".join( [ L[0].split( '_' )[0] ] + L[1:] )
	f.close()
