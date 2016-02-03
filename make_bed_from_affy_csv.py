#!/usr/bin/env python
from __future__ import division
import sys
import csv

try:
	ifn = sys.argv[1]
except IndexError:
	print >> sys.stderr, """\
Script to create a BED file from Affymetrix Netaffx CSV file.

usage: ./script <infile>

example:
./script.py file.csv"""
	sys.exit( 0 )

unwanted = [ '#', 'p', '"', '-' ]

c = 0
with open( ifn ) as f:
	g = open( ".".join( ifn.split( '.' )[:-1] ) + ".bed", 'w' )
	for row in csv.reader( f, delimiter=",", quotechar='"' ):
		if c > 5: break
		if row[0][0] in unwanted: continue
		if row[1][0] in unwanted: continue
		print >> g, "\t".join( [ row[1], row[3], row[4], row[0], '-', row[2] ] )
		c += 0
	g.close()
