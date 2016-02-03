#!/usr/bin/env python
import sys
import re
import random
import multiprocessing

def fitler( expr_fn, fitl_fn, cols, name_map ):
	
	f = open( expr_fn )
	g = open( fitl_fn, 'w' )
	sample_names = list()
	for row in f:
		L = row.strip().split( '\t' )
		if L[1] == 'NA': continue # if 'NA' occurs for this gene, exclude the gene entirely
		if row[0] == 'g':
			sample_names = L[1:]
			print >> g, "\t".join( [ "gene_id" ] + [ name_map[ c ] for c in cols ] )
			continue
		
		expr = dict( zip( sample_names, L[1:] ))
		
		this_row = [ expr[ name_map[c] ] for c in cols ]
		
		print >> g, "\t".join( [ L[0] ] + this_row )
	
	g.close()
	f.close()
	
	return

try:
	no_side = int( sys.argv[1] )
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <no_of_samples per side>"
	sys.exit( 1 )

# get the names of all Whole Blood samples
#f = open( "/data2/paulk/MARS/F_GTEx/Affymetrix_Expression_sample_names_by_tissue/Whole_Blood.txt" )
f = open( "/data2/paulk/MARS/F_GTEx/Affymetrix_Expression_sample_names_by_tissue/Heart_-_Left_Ventricle.txt" )
wb_names = [ row.strip() for row in f ]
f.close()

# get the names of all Muscle Skeletal samples
f = open( "/data2/paulk/MARS/F_GTEx/Affymetrix_Expression_sample_names_by_tissue/Muscle_-_Skeletal.txt" )
ms_names = [ row.strip() for row in f ]
f.close()

# get the names of all samples in the RMA/prediction
f = open( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/GTEx_rma.txt" )
predicted_names = [ nom[:-4].split( '_' )[1] for nom in f.readline().strip().split( '\t' )[1:] ]
f.close()

# get a map of sample names to Affymetrix sample names
f = open( "/data2/paulk/MARS/F_GTEx/Affymetrix_sample_names.txt" )
manes_nap = dict()
for row in f:
	L = row.strip().split( '\t' )
	manes_nap[ L[0] ] = L[1]
f.close()

# find the subset of samples that are Whole Blood
wb_predicted_names = [ n for n in wb_names if n in predicted_names ]

# find the subset of samples that are Muscle Skeletal
ms_predicted_names = [ n for n in ms_names if n in predicted_names ]

# pick 10 WB at random
wb_random = list()
for i in xrange( no_side ):
	nom = random.choice( wb_predicted_names )
	wb_predicted_names.remove( nom )
	wb_random.append( nom )

# pick 10 MS at random
ms_random = list()
for i in xrange( no_side ):
	nom = random.choice( ms_predicted_names )
	ms_predicted_names.remove( nom )
	ms_random.append( nom )

# filter data
print >> sys.stderr, "Starting first filter process ... ",
#p1 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/DE_results_10_samples/GTEx_rma_filtered.txt", wb_random+ms_random, manes_nap, ) )
p1 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/GTEx_rma.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/QN_subset_rnaseq/DE_Heart_Skeletal_"+str( 2*no_side )+"_samples/GTEx_rma_filtered.txt", wb_random+ms_random, manes_nap, ) )
p1.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting second filter process ... ",
#p2 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/GTEx_malte.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/DE_results_10_samples/GTEx_predictions_filtered.txt", wb_random+ms_random, manes_nap, ) )
p2 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/QN_subset_rnaseq/GTEx_malte_QRF_QN_QN.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/QN_subset_rnaseq/DE_Heart_Skeletal_"+str( 2*no_side )+"_samples/GTEx_malte_QRF_filtered.txt", wb_random+ms_random, manes_nap, ) )
p2.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting third filter process ... ",
#p3 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/GTEx_rnaseq.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/DE_results_10_samples/GTEx_rnaseq_filtered.txt", wb_random+ms_random, manes_nap, ) )
p3 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/QN_subset_rnaseq/GTEx_rnaseq_QN.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/QN_subset_rnaseq/DE_Heart_Skeletal_"+str( 2*no_side )+"_samples/GTEx_rnaseq_filtered.txt", wb_random+ms_random, manes_nap, ) )
p3.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting fourth filter process ... ",
#p4 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/GTEx_plier.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/DE_results_10_samples/GTEx_plier_filtered.txt", wb_random+ms_random, manes_nap, ) )
p4 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/GTEx_plier.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2/QN_subset_rnaseq/DE_Heart_Skeletal_"+str( 2*no_side )+"_samples/GTEx_plier_filtered.txt", wb_random+ms_random, manes_nap, ) )
p4.start()
print >> sys.stderr, "OK"

p1.join()
print >> sys.stderr, "Process 1 has returned."
p2.join()
print >> sys.stderr, "Process 2 has returned."
p3.join()
print >> sys.stderr, "Process 3 has returned."
p4.join()
print >> sys.stderr, "Process 4 has returned."


