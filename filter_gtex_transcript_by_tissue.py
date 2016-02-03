#!/usr/bin/env python
from __future__ import division
import sys
import scipy
import gzip

def clean_sample_names( row1 ):
	L = row1.strip().split( '\t' )
	cleaned_sample_names = [ k[:-4] for k in L[4:] ]
	return cleaned_sample_names

if __name__ == "__main__":
	try:
		tissue_fn = sys.argv[1]
		rnaseq_fn = sys.argv[2]
	except IndexError:
		print >> sys.stderr, """\
Script to reduce the column of the RNA-Seq gene to only those provided in the tissue.

usage: ./script.py <tissue_fn> <rnaseq_gene_fn>

example:
mkdir Whole_Blood
./script.py mapped_sample_names_by_tissue/Whole_Blood.txt GTEx_transcript_rpkm.gct > Whole_Blood/GTEx_transcript_rpkm.txt"""
		sys.exit( 1 )
	
	f = open( tissue_fn ) # tissue
	tissue_samples = [ row.strip().split( '\t' )[0] for row in f ]
	f.close()
	
	# some flags
	printed_header = False	# I have printed the output header
	reported = False				# I have already told you which samples are missing

	f = gzip.open( rnaseq_fn ) # rna-seq gene
	c = 0
	for row in f:
		if c > 10: break
		L = row.strip().split( '\t' )
		 
		# how to process header
		if row[:4] == 'Targ':
			sample_names = clean_sample_names( row )
			continue
	
		# how to process row
		this_row = dict( zip( sample_names, L[4:] ))
	
	#	print this_row
	
		# filter
		to_print = list()
		to_print += [ L[0].split( '.' )[0] ]
		header = [ "tx_id" ]
		for s in tissue_samples:
			try:
				expr = this_row[s]
			except KeyError:
				if not reported:
					print >> sys.stderr, "Missing RNA-Seq data for %s" % s
				continue
	
			to_print += [ expr ]
			header += [ s ]

		reported = True
		
		# print the header
		if not printed_header:
			print "\t".join( header )
			printed_header = True
	
		# print the row	
		print "\t".join( to_print )
		c += 0
	
	f.close()	
