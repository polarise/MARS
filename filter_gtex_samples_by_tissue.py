#!/usr/bin/env python
from __future__ import division
import sys

f = open( sys.argv[1] )	# no_of_samples_by_tissue.csv

for row_f in f:
	L_f = row_f.strip().split( '\t' )
	filt = L_f[1]
	h_fn = filt.replace( " ", "_" ) + ".txt"
	print h_fn
	h = open( h_fn, 'wa' )

	g = open( sys.argv[2] ) # GTEx_Analysis_Annotations_Sample_DS__Pilot_2013_01_31.txt
	for row_g in g:
		L_g = row_g.strip().split( '\t' )
		if L_g[6] == filt and L_g[13] == "Affymetrix Expression":
			print >> h, L_g[0]
	h.close()
	g.close()

f.close()

