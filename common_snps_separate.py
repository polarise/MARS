#!/usr/bin/env python
from __future__ import division
import sys

def ids_F( file_F, p_max_f ):
	ids_S = set()
	ids_D = dict()
	for row_s in file_F:
		if row_s[0] == 's': continue
		L = row_s.strip().split('\t')
		if L[-1] == 'NA': continue
		if float( L[-1] ) < p_max_f:
			ids_S.add( L[0]+"_"+L[1] )
			ids_D[L[0]+"_"+L[1]] = str( float( L[-1] ))
		else:
			break # don't waste time; the file is ordered
	return ids_S,ids_D

verbose = False
detailed = False

# a header
if not detailed:
	print "FDR\trs.total\taf.total\tin.total\tra.count\tri.count\tai.count"

for p_max_f in [ 1e-60, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]:
	# open the rnaseq file and get the rnaseq SNP ids
	f = open( "/home/paulk/MARS/4_association_tests/rnaseq_10k/association_results.txt.qval" )
	rnaseq_ids_S, rnaseq_ids_D = ids_F( f, p_max_f )
	f.close()

	# open the affy file and get the affy SNP ids
	f = open( "/home/paulk/MARS/4_association_tests/affy_comprnaseq/association_results.txt.qval" )
	affy_ids_S, affy_ids_D = ids_F( f, p_max_f )
	f.close()

	# open the file of interest and get the CRF SNP ids
	f = open( "/home/paulk/MARS/4_association_tests/seekarray_comprnaseq/association_results.txt.qval" )
	interest_ids_S, interest_ids_D = ids_F( f, p_max_f )
	f.close()

	# count those in common
	ra_S = rnaseq_ids_S.intersection( affy_ids_S )
	ri_S = rnaseq_ids_S.intersection( interest_ids_S )
	ai_S = affy_ids_S.intersection( interest_ids_S )
	print "%.2f\t%d\t%d\t%d\t%d\t%d\t%d" % ( p_max_f*100, len(rnaseq_ids_S), len(affy_ids_S), len(interest_ids_S), len(ra_S), len(ri_S), len(ai_S) )
