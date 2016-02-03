#!/home/paulk/software/bin/python
"""
given gene coordinates file and some augmentation value generated augmented coordinates in BED format
"""
from __future__ import division
import sys


try:
	overlappingfn_s, genotypenofn_s, expressionfn_s = sys.argv[1:4]
except ValueError:
	print >> sys.stderr,"Usage: ./script.py <overlappingfn_s> <genotypenofn_s> <expressionfn_s>"
	sys.exit(1)

D1_D = dict()
f = open(overlappingfn_s)
for row_s in f:
	D1_D[ tuple(row_s.strip().split( '\t' )) ] = []
f.close()

D2_D = dict()
f = open(genotypenofn_s)
for row_s in f:		# assume uniqueness; should be unique anyway!
	L = row_s.strip().split( '\t' )
	D2_D[ L[0] ] = L[1:]
f.close()

D3_D = dict()
f = open(expressionfn_s)
for row_s in f:
	if row_s == 'p': continue
	L = row_s.strip().split( '\t' )
	D3_D[ L[0] ] = L[1:]
f.close()

# now combine all data
for d_T in D1_D:
	try:
		genotype_L = D2_D[ d_T[0] ]
		genotype_b = True
	except KeyError:
		genotype_b = False
	
	try:
		expression_L = D3_D[ d_T[1] ] 
		expression_b = True
	except KeyError:
		expression_b = False
	
	if not genotype_b:
		print >> sys.stderr,"Missing genotype for %s @ %s" % d_T
		continue
	elif not expression_b:
		print >> sys.stderr,"Missing expression for %s @ %s" % d_T
		continue
	else:
		print "\t".join( list(d_T) + genotype_L + [ ",".join( expression_L ) ] )
