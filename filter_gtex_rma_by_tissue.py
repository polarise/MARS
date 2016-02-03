#!/usr/bin/env python
from __future__ import division
import sys

if __name__ == "__main__":
	try:
		samples_fn = sys.argv[1]
		affy_ids_fn = sys.argv[2]
		rma_fn = sys.argv[3]
	except IndexError:
		print >> sys.stderr, """\
Script to reduce the column of the RMA genes to only those provided in the tissue.

usage: ./script.py <samples_fn> <affy_ids_fn> <rma_fn>

example:
./script.py Whole_Blood/samples.txt Affymetrix_sample_names.txt RMA/GSE45878_series_matrix.txt > Whole_Blood/GTEx_rma.txt"""
		sys.exit( 0 )


	# get the test samples in order
	f = open( samples_fn ) # samples.txt
	affy_samples = list()
	for row in f:
		if row[:3] == 'hts': continue
		L = row.strip().split( '\t' )
		if L[0][0] == '*': affy_samples.append( L[1] )
	f.close()

	f = open( affy_ids_fn ) # Affymetrix_sample_names.txt
	affy_ids = dict()
	for row in f:
		L = row.strip().split( '\t' )
		key = L[1].split( '_' )[0]
		affy_ids[ L[1] ] = key
	f.close()


	c = 0
	f = open( rma_fn ) # Affy RMA expr
	header = [ 'gene_id' ]
	for row in f:
		if c > 5: break
		L = row.strip().split( '\t' )
		if L[0] == 'gene_id':
			sample_names = L[1:]
			for s in affy_samples:
				header += [ s ]
		
			# print header
			print "\t".join( header )
			continue
	
		expr = dict( zip( sample_names, L[1:] ) )
	
		this_row = [ L[0] ]
		for s in affy_samples:
			other_name = affy_ids[ s ]
			this_row += [ expr[ other_name ] ]
	
		print "\t".join( this_row )
		c += 0
	f.close()		
	
	
		
