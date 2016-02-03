#!/usr/bin/env python
from __future__ import division
import sys
import scipy

try:
	fn = sys.argv[1]
except IndexError:
	print >> sys.stderr, """\
Script to compute coefficient of variation.

usage: ./script.py <filename>

example: ./script.py GTEx_gene_rpkm.gct"""
	sys.exit( 1 )

unwanted = [ 'N' ]

f = open( fn )
c = 0
for row in f:
	if c > 2000: break
	if row[0] in unwanted: continue
	L = row.strip().split( '\t' )
	vals = map( float, L[2:] )
	if sum( vals ) == 0:
		print >> sys.stderr, "Zero row: %s" % L[0]
		continue
	sd = scipy.std( vals )
	mn = scipy.mean( vals )
	cov = sd/mn
	print "%s\t%s" % ( L[0], cov )	
	c += 1
f.close()
	
