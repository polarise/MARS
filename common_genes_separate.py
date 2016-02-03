#!/usr/bin/env python
from __future__ import division
import sys
import scipy

def ids_F( file_F, p_max_f ):
	ids_S = set()
	ids_D = dict()
	for row_s in file_F:
		if row_s[0] == 'I': continue
		L = row_s.strip().split('\t')
		j = -3
		if L[j] == 'NA': continue
		if float( L[j] ) <= p_max_f:
			ids_S.add( L[0] )
			ids_D[L[0]] = str( float( L[j] ))
		else:
			break # don't waste time; the file is ordered
	return ids_S,ids_D

detailed = False

# a header
if not detailed:
	print "FDR\trs.total\taf.total\tin.total\tra.count\tri.count\tai.count"

#for p_max_f in [ 1e-60, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]:
for p_max_f in [1e-60, 1e-50, 1e-40, 1e-30, 1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3 ] + list(scipy.linspace(2e-3,0.01,20)) + list(scipy.linspace(0.01,0.25,20)) + [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]:
	# open the rnaseq file and get the rnaseq SNP ids
#	f = open( "/home/paulk/MARS/8_pop_diffexpr/0A_raw_common_expr/MIX_common2_rnaseq.txt_simple_diffex" )
#	f = open( "/home/paulk/MARS/8_pop_diffexpr/0A_raw_common_expr/MIX_common_rnaseq.txt_diffex" )
	f = open( "/home/paulk/MARS/6_sex_diffexpr/mine_results/1B_simple_diffext/CEU_common_rnaseq.txt_simple_diffex" )
	rnaseq_ids_S, rnaseq_ids_D = ids_F( f, p_max_f )
	f.close()

	# open the affy file and get the affy SNP ids
#	f = open( "/home/paulk/MARS/8_pop_diffexpr/0A_raw_common_expr/MIX_common2_affy.txt_simple_diffex" )
#	f = open( "/home/paulk/MARS/8_pop_diffexpr/0A_raw_common_expr/MIX_common_affy.txt_diffex" )
	f = open( "/home/paulk/MARS/6_sex_diffexpr/mine_results/1B_simple_diffext/CEU_common_affy.txt_simple_diffex" )
	affy_ids_S, affy_ids_D = ids_F( f, p_max_f )
	f.close()

	# open the file of interest and get the CRF SNP ids
#	f = open( "/home/paulk/MARS/8_pop_diffexpr/0A_raw_common_expr/MIX_common2_cforest.txt_simple_diffex" )
#	f = open( "/home/paulk/MARS/8_pop_diffexpr/0A_raw_common_expr/MIX_common_cforest.txt_diffex" )
	f = open( "/home/paulk/MARS/6_sex_diffexpr/mine_results/1B_simple_diffext/CEU_common_cforest.txt_simple_diffex" )
	interest_ids_S, interest_ids_D = ids_F( f, p_max_f )
	f.close()

	# count those in common
	ra_S = rnaseq_ids_S.intersection( affy_ids_S )
	ri_S = rnaseq_ids_S.intersection( interest_ids_S )
	ai_S = affy_ids_S.intersection( interest_ids_S )
#	print "%f\t%d\t%d\t%d\t%d\t%d\t%d" % ( p_max_f*100, len(rnaseq_ids_S), len(affy_ids_S), len(interest_ids_S), len(ra_S), len(ri_S), len(ai_S) )
	print "%e\t%d\t%d\t%d\t%d\t%d\t%d" % ( p_max_f, len(rnaseq_ids_S), len(affy_ids_S), len(interest_ids_S), len(ra_S), len(ri_S), len(ai_S) )
