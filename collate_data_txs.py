#!/usr/bin/env python
from __future__ import division
import sys
import gzip
import re
#import cPickle

unwanted = [ 't', '#' ]

def get_sample_names( f ):
	"""
	process the file 'samples.txt' to extract sample names
	"""
	train_hts_L = list()
	test_hts_L = list()
	train_ma_L = list()
	test_ma_L = list()

	for row in f:
		L = row.strip().split( '\t' )
		if row[0] == 'h': # skip the header
			continue
		if L[0] == 'NA':
			test_ma_L += [L[1]]
		elif L[0][0] == '*':
			test_hts_L += [L[0]]
		else:
			train_hts_L += [L[0]]
			train_ma_L += [L[1]]
	
	return train_hts_L, test_hts_L, train_ma_L, test_ma_L


try:
	id_fn = sys.argv[1]
	gene_to_txs_fn = sys.argv[2]
	train_data_fn = sys.argv[3]
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <sample_names> <gene_to_txs_fn> <train_data_fn>"
	sys.exit( 1 )

if re.search( r"gz$", id_fn ):
	f = gzip.open( id_fn )
elif re.search( r"txt$", id_fn ):
	f = open( id_fn )
else:
	raise ValueError( "Invalid input file: should be .txt or .gz text file" )
	sys.exit( 1 )
	
train_hts_L, test_hts_L, train_ma_L, test_ma_L = get_sample_names( f )
f.close()

samples = train_hts_L + test_hts_L
no_samples = len( samples )

tx_expr = dict()
# get the transcript estimates for all samples
for sample in samples:
	try:
		f = open( "/data1/paulk/Montgomery_Data/output/%s/tophat_out/isoforms.fpkm_tracking" % sample )
	except IOError:
		f = open( "/data1/paulk/%s/isoforms.fpkm_tracking" % sample )
	for row in f:
		if row[0] in unwanted: continue
		L = row.strip().split( '\t' )
		if L[0] not in tx_expr:
			tx_expr[L[0]] = [ L[9] ]
		else:
			tx_expr[L[0]] += [ L[9] ]	
	f.close()

print "tx_id\t%s" % "\t".join( samples )
for tx in tx_expr:
	print "%s\t%s" % ( tx, "\t".join( tx_expr[tx] ) )


sys.exit( 0 )

# map of gene to transcripts
#gene2txs = cPickle.load( open( gene_to_txs_pic_fn ) )

gene2txs = dict()
f = open( gene_to_txs_fn )
for row in f:
	if row[:3] == 'Ens': print >> sys.stderr, "Hit header..."; continue
	L = row.strip().split( '\t' )
	if L[0] not in gene2txs:
		gene2txs[L[0]] = [ L[1] ]
	else:
		gene2txs[L[0]] += [ L[1] ]
f.close()

f = open( train_data_fn )
c = 0
for row in f:
	if c > 5: break
	L = row.strip().split( '\t' )
	
	# get the gene name
	gene = L[0]
	# for each gene get the transcripts associated
	txs = gene2txs[gene]
	# check whether the transcript has an expression estimate
	tx_ids = list()
	present_tx_expr = list()
	for tx in txs:
	# take it if present
		if tx in tx_expr:
	# also, collect the transcript ids
			tx_ids += [tx]
			present_tx_expr += tx_expr[tx]
	# print it all at once
	if len( tx_ids ) > 0:
		print "\t".join( [ gene, ",".join( tx_ids ), str( len( tx_ids ) ), str( no_samples ), L[2], L[3], ",".join( present_tx_expr ), L[5] ] )
	c += 0
