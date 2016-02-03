#!/home/paulk/software/bin/python
"""
given gene coordinates file and some augmentation value generated augmented coordinates in BED format
"""
from __future__ import division
import sys

try:
	genecoordsfn_s = sys.argv[1]
	aug_i = int( sys.argv[2] )
except IndexError:
	print >> sys.stderr,"Usage: ./script.py <genecoordsfn_s> <aug_i>"
	sys.exit( 1 )

f = open( genecoordsfn_s )
for row_s in f:
	L = row_s.strip().split( '\t' )
	chrom_s, st_s, sp_s, ID_s, other, sd_s = L
	if chrom_s[0] != 'c':
		chrom_s = "chr" + chrom_s
	print "\t".join( [ chrom_s, str( max([0,int(st_s)-aug_i])), str(int(sp_s)+aug_i), ID_s, "-", sd_s ] )
f.close()
