#!/usr/bin/env python
from __future__ import division
import sys

# find duplicates
gene_to_no_probes = dict()
f = open( "/data2/paulk/MARS/YRI_vs_CEU_mine/Ens_to_no_probes.txt" )
for row_s in f:
	L = row_s.strip().split( '\t' )
	if L[0] not in gene_to_no_probes:
		gene_to_no_probes[L[0]] = [int( L[1] )]
	else:
		gene_to_no_probes[L[0]] += [ int( L[1] )]
f.close()

# filter the training/testing data file
f = open( sys.argv[1] )
genes_passed = dict()
for row_s in f:
	L = row_s.strip().split( '\t' )
	if int( L[2] ) == max( gene_to_no_probes[L[0][:15]] ) and L[0][:15] not in genes_passed:
		print "\t".join( [ L[0][:15] ] + L[1:] )
		genes_passed[ L[0][:15] ] = 0
	else:
		print >> sys.stderr, L[0][:15]
f.close()

