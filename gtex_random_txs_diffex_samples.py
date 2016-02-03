#!/usr/bin/env python
import sys
import re
import random
import multiprocessing
import gzip

def fitler( expr_fn, fitl_fn, cols, names_map, zipped=False ):
	
	if zipped:
		f = gzip.open( expr_fn )
	else:
		f = open( expr_fn )
	
	g = open( fitl_fn, 'w' )
	sample_names = list()
	for row in f:
		L = row.strip().split( '\t' )
		if row[0] == 't':
			sample_names = L[1:]
#			if zipped:
#				print >> g, "\t".join( [ "tx_id" ] + [ c for c in cols ] )
#			else:
			print >> g, "\t".join( [ "tx_id" ] + [ names_map[ c ] for c in cols ] )
			continue
		
		expr = dict( zip( sample_names, L[1:] ))
		
		if zipped:
			this_row = [ expr[ c ] for c in cols ]
		else:
			this_row = [ expr[ names_map[c] ] for c in cols ]
		
		print >> g, "\t".join( [ L[0] ] + this_row )
	
	g.close()
	f.close()
	
	return

# get the names of all Whole Blood samples
f = open( "/data2/paulk/MARS/F_GTEx/mapped_sample_names_by_tissue/Whole_Blood.txt" )
wb_names = [ row.strip().split( '\t' )[0] for row in f ]
f.close()

# get the names of all Muscle Skeletal samples
f = open( "/data2/paulk/MARS/F_GTEx/mapped_sample_names_by_tissue/Muscle_-_Skeletal.txt" )
ms_names = [ row.strip().split( '\t' )[0] for row in f ]
f.close()

# get the names of all samples in the RMA/prediction
manes_nap = dict()
predicted_names = list()
f = open( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/samples.txt" )
for nom in f:
	if nom[0] == "*":
		L = nom.strip().split( '\t' )
		predicted_names.append( L[0][1: ] )
		manes_nap[ L[0][1:] ] = L[1]
f.close()

# find the subset of samples that are Whole Blood
wb_predicted_names = [ n for n in wb_names if n in predicted_names ]

# find the subset of samples that are Muscle Skeletal
ms_predicted_names = [ n for n in ms_names if n in predicted_names ]

# pick 10 WB at random
wb_random = list()
for i in xrange(5):
	nom = random.choice( wb_predicted_names )
	wb_predicted_names.remove( nom )
	wb_random.append( nom )

# pick 10 MS at random
ms_random = list()
for i in xrange(5):
	nom = random.choice( ms_predicted_names )
	ms_predicted_names.remove( nom )
	ms_random.append( nom )

#print manes_nap
#print
#print wb_random+ms_random
#print
#for a in wb_random+ms_random:
#	print a, manes_nap[ a ]
#sys.exit( 0 )

# filter data

print >> sys.stderr, "Starting first filter process ... ",
p1 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_transcript_rpkm.txt.gz", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_transcript_rpkm_filtered.txt", wb_random+ms_random, manes_nap, ), kwargs={ "zipped": True } )
p1.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting second filter process ... ",
p2 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_thresh_0.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_filtered_thresh_0.txt", wb_random+ms_random, manes_nap, ) )
p2.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting third filter process ... ",
p3 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_thresh_1.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_filtered_thresh_1.txt", wb_random+ms_random, manes_nap, ) )
p3.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting fourth filter process ... ",
p4 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_thresh_2.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_filtered_thresh_2.txt", wb_random+ms_random, manes_nap, ) )
p4.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting fifth filter process ... ",
p5 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_thresh_3.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_filtered_thresh_3.txt", wb_random+ms_random, manes_nap, ) )
p5.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting sixth filter process ... ",
p6 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_thresh_4.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_filtered_thresh_4.txt", wb_random+ms_random, manes_nap, ) )
p6.start()
print >> sys.stderr, "OK"

print >> sys.stderr, "Starting seventh filter process ... ",
p7 = multiprocessing.Process( target=fitler, args=( "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_thresh_5.txt", "/data2/paulk/MARS/F_GTEx/All_Tissues_0.2_OOB/GTEx_txs_predictions_filtered_thresh_5.txt", wb_random+ms_random, manes_nap, ) )
p7.start()
print >> sys.stderr, "OK"

p1.join()
print >> sys.stderr, "Process 1 has returned."
p2.join()
print >> sys.stderr, "Process 2 has returned."
p3.join()
print >> sys.stderr, "Process 3 has returned."
p4.join()
print >> sys.stderr, "Process 4 has returned."
p5.join()
print >> sys.stderr, "Process 5 has returned."
p6.join()
print >> sys.stderr, "Process 6 has returned."
p7.join()
print >> sys.stderr, "Process 7 has returned."

