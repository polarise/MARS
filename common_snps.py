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
			ids_S.add( L[0] )
			ids_D[L[0]] = str( float( L[-1] ))
		else:
			break # don't waste time; the file is ordered
	return ids_S,ids_D

try:
	method = sys.argv[1]
	chrom = sys.argv[2]
	


#p_max_f = float( sys.argv[3] )
verbose = False
detailed = False

# a header
if not detailed:
	print "FDR\trs.total\taf.total\tin.total\tra.count\tri.count\tai.count"

for p_max_f in [ 1e-60, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]:
	# open the rnaseq file
	f = open( "/home/paulk/MARS/4_association_tests/rnaseq/chr%s/association_results.txt.qval" % chrom )
	rnaseq_ids_S, rnaseq_ids_D = ids_F( f, p_max_f )
	f.close()

	# open the affy file
	f = open( "/home/paulk/MARS/4_association_tests/affy/chr%s/association_results.txt.qval" % chrom )
	affy_ids_S, affy_ids_D = ids_F( f, p_max_f )
	f.close()

	# open the file of interest
	f = open( "/home/paulk/MARS/4_association_tests/seekarray_%s/chr%s/association_results.txt.qval" % ( sys.argv[1], chrom ))
	interest_ids_S, interest_ids_D = ids_F( f, p_max_f )
	f.close()

	# count those in common
	if detailed:
		print "========================================================================="
		print "Chromosome: %s" % sys.argv[2]
		print "Number of SNPs with FDR < %.2f%% found in RNA-Seq: %d" % ( p_max_f*100, len(rnaseq_ids_S) )
		print
		ra_S = rnaseq_ids_S.intersection( affy_ids_S )
		if len( ra_S ) != 0 and verbose:
			print "Common between RNA-Seq and RMA (Affymetrix): %d of %d" % ( len(ra_S), len(affy_ids_S) )
			print "-------------------------------------------------------------------------"
			print "\t".join( [ "snp_id", "RS_FDR", "RMA_FDR" ] )
			for id_s in ra_S:
				print "\t".join( [ id_s, rnaseq_ids_D[id_s], affy_ids_D[id_s] ] )
		elif len( ra_S ) != 0 and not verbose:
			print "Common between RNA-Seq and RMA (Affymetrix): %d of %d" % ( len(ra_S), len(affy_ids_S) )
		else:
			print "No SNPs found in common between RNA-Seq and RMA (Affymetrix)"
		#print

		ri_S = rnaseq_ids_S.intersection( interest_ids_S )
		if len( ri_S ) != 0 and verbose:
			print "Common between RNA-Seq and %s: %d of %d" % ( sys.argv[1].upper(), len(ri_S), len(interest_ids_S) )
			print "-------------------------------------------------------------------------"
			print "\t".join( [ "snp_id", "RS_FDR", "%s_FDR" % sys.argv[1].upper() ] )
			for id_s in ri_S:
				print "\t".join( [ id_s, rnaseq_ids_D[id_s], interest_ids_D[id_s] ] )
		elif len( ri_S ) != 0 and not verbose:
			print "Common between RNA-Seq and %s: %d of %d" % ( sys.argv[1].upper(), len(ri_S), len(interest_ids_S) )
		else:
			print "No SNPs found in common between RNA-Seq and %s" % sys.argv[1].upper()
		#print

		ai_S = affy_ids_S.intersection( interest_ids_S )
		if len( ai_S ) != 0 and verbose:
			print "Common between RMA (Affymetrix) and %s: %d" % ( sys.argv[1].upper(), len(ai_S) )
			print "-------------------------------------------------------------------------"
			print "\t".join( [ "snp_id", "RMA_FDR", "%s_FDR" % sys.argv[1].upper() ] )
			for id_s in ai_S:
				print "\t".join( [ id_s, affy_ids_D[id_s], interest_ids_D[id_s] ] )
		elif len( ai_S ) != 0 and not verbose:
			print "Common between RMA (Affymetrix) and %s: %d" % ( sys.argv[1].upper(), len(ai_S) )	
		else:
			print "No SNPs found in common between RMA (Affymetrix) and %s" % sys.argv[1].upper()
		print
	else:
		ra_S = rnaseq_ids_S.intersection( affy_ids_S )
		ri_S = rnaseq_ids_S.intersection( interest_ids_S )
		ai_S = affy_ids_S.intersection( interest_ids_S )
		print "%.2f\t%d\t%d\t%d\t%d\t%d\t%d" % ( p_max_f*100, len(rnaseq_ids_S), len(affy_ids_S), len(interest_ids_S), len(ra_S), len(ri_S), len(ai_S) )
