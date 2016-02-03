#!/usr/bin/env python
from __future__ import division
import sys
import cPickle

if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "usage: ./script.py <input-picfile>"
		sys.exit( 0 )
	
	with open( fn ) as f:
		data = cPickle.load( f )
		for d in data:
			print "\t".join( [ d ] + map( str, data[d] ))	