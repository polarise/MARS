#!/usr/bin/env python
from __future__ import division
import sys

try:
	tissue_fn = sys.argv[1]
	affy_ids_fn = sys.argv[2]
	affy_expr_fn = sys.argv[3]
except IndexError:
	print >> sys.stderr, """\
Script to reduce the column of the array to only those provided in the tissue.

usage: ./script.py <tissue_fn> <affy_ids_fn> <affy_expr_fn>

example:
# assumes that Whole_Blood directory already exists
./script.py mapped_sample_names_by_tissue/Whole_Blood.txt Affymetrix_sample_names.txt GTEx_probe_intensities.txt > Whole_Blood/GTEx_probe_intensities.txt"""
	sys.exit( 0 )

f = open( sys.argv[1] ) # tissue
tissue_samples = [ row.strip().split( '\t' )[1] for row in f ]
f.close()

# some flags
printed_header = False	# I have printed the output header
reported = False				# I have already told you which samples are missing

f = open( sys.argv[2] ) # Affymetrix id file
affy_ids = dict()
for row in f:
	L = row.strip().split( '\t' )
	affy_ids[ L[0] ] = L[1]
f.close()

unwanted = [ '#' ]

f = open( sys.argv[3] ) # affy probe intensities
for row in f:
	if row[0] in unwanted: continue
	L = row.strip().split( '\t' )
	 
	# hwo to process header
	if row[:8] == 'probe_id':
		sample_names = L[7:]
		continue
	
	# how to process row
	this_row = dict( zip( sample_names, L[7:] ))
	
#	print this_row
	
	# filter
	to_print = list()
	to_print += [ L[0], L[4] ]
	header = [ "probe_id", "probeset_id" ]
	for s in tissue_samples:
		try:
			expr = this_row[ affy_ids[ s ] ] # use the map of affy ids
		except KeyError:
			if not reported:
				print >> sys.stderr, "Missing probe intensity data for %s" % s
			continue
	
		to_print += [ expr ]
		header += [ affy_ids[ s ] ]

	reported = True
		
	# print the header
	if not printed_header:
		print "\t".join( header )
		printed_header = True
	
	# print the row	
	print "\t".join( to_print )
	
f.close()


