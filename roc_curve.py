#!/usr/bin/env python
from __future__ import division
import sys
import scipy

unwanted = ['I', 's']

def gene_filter( gene_list, thresh ):
	result = list()
	for g in gene_list:
		if g[1] <= thresh:
			result += [ g[0] ]
#		else:
#			break # because they are ordered; ignore ties (>:-()
	return result

def g2p( fn, col ):
	f = open( fn )
	names = list()
	data_list = list()
	i = 1
	for row in f:
		if row[0] in unwanted: continue
		L = row.strip().split( '\t' )
		if L[col] == 'NA': continue
		if col >= len( L ):
			raise ValueError( "Value 'col' is larger than the number of columns." )
		names += [L[0]]
#		data_list += [( L[0], float( L[col] ) )]
		data_list += [( L[0], i )]
		i += 1
	f.close()
	return names, data_list

if __name__ == "__main__":
	try:
		rnaseq_fn = sys.argv[1]
		other_fn = sys.argv[2]
		col = int( sys.argv[3] )-1	# 1- to 0-based indexing
	except IndexError:
		print >> sys.stderr, "Usage: ./script.py <rnaseq_fn> <other_fn> <col>\nIt is assumed that the file is sorted by q-values"
		sys.exit( 1 )

	rnaseq, rnaseq_list =  g2p( rnaseq_fn, col )
	other, other_list = g2p( other_fn, col )

#	qvals = [ r[1] for r in rnaseq_list ]
#	qvals.sort()

	if 0 < len( qvals ) <= 500:
		step = 1
	elif 500 < len( qvals ) <= 2000:
		step = 10
	elif 2000 < len( qvals ) <= 10000:
		step = 20
	elif 10000 < len( qvals ) <= 50000:
		step = 25
	elif len( qvals ) > 50000:
		step = 30
	else:
		print >> sys.stderr, "Emtpy list! Exciting..."
		sys.exit( 1 )

	c = 0
	TPR = 0
	FPR = 0
	for q in xrange( 1, len( qvals )-1, step ): 
#		q = qvals[i]
#		q = i
		if c > 100: break
		rnaseq_pos = set( gene_filter( rnaseq_list, q ))	# positives
		rnaseq_compl = set( rnaseq ).difference( rnaseq_pos )	# negatives
	
		other_pos = set( gene_filter( other_list, q ) ) # declared positives
		other_compl = set( other ).difference( other_pos )	# declared negatives
	
		TP = len( rnaseq_pos.intersection( other_pos ))	# number of those declared positive that are really positive
		FP = len( rnaseq_compl.intersection( other_pos ))	# number of those declared positive that are really negative
		TN = len( rnaseq_compl.intersection( other_compl ))	# number of those declared negative that are really negative
		FN = len( rnaseq_pos.intersection( other_compl ))	# number of those declared negative that are really positive
	
		c += 0
		TPR = max([TP/(TP + FN),TPR]) # non-decreasing
		
		FPR = max([FP/(FP + TN),FPR])
	
		print >> sys.stdout, "\t".join( map( str, [FPR, TPR] ))
		
		sys.stdout.flush()
