#!/usr/bin/env python
from __future__ import division
import sys

unwanted = [ '#' ]

try:
	id_fn = sys.argv[1]
	ma_fn = sys.argv[2]
	hts_fn = sys.argv[3]
	g2ps_fn = sys.argv[4]
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <id_fn> <ma_fn> <hts_fn> <g2ps_fn>"
	sys.exit( 1 )

# load data
# sample IDs
train_hts_L = list()
test_hts_L = list()
train_ma_L = list()
test_ma_L = list()
f = open( id_fn )
for row in f:
	L = row.strip().split( '\t' )
	if row[0] == 'h':
		continue
	if L[0] == 'NA':
		test_ma_L += [L[1]]
	elif L[0][0] == '*':
		test_hts_L += [L[0]]
	else:
		train_hts_L += [L[0]]
		train_ma_L += [L[1]]
f.close()	

# microarray data
no_samples = 0
ma_names_L = list()
ma_train_D = dict()
ma_test_D = dict()
f = open( ma_fn )
for row in f:
	if row[0] in unwanted: continue
	L = row.strip().split( '\t' )
	if row[0] == 'p':
		ma_names_L = L[2:]
		print "\t".join( [ L[0], L[1] ] + train_ma_L )
		continue
	ps = int( L[1] )
	p = int( L[0] )
	if no_samples == 0:
		no_samples = len( L[2:] )
	
	ints_D = dict( zip( ma_names_L, L[2:] ) )
	
	if ps not in ma_train_D:
		ma_train_D[ps] = dict()
		ma_test_D[ps] = dict()
	
	print "\t".join( [ L[0], L[1] ] + [ ints_D[k] for k in train_ma_L ] )
	#ma_train_D[ps][p] = { k: ints_D[k] for k in train_ma_L }
	#ma_test_D[ps][p] = { k: ints_D[k] for k in test_ma_L }
f.close()

print >> sys.stderr, "Finished reading in microarray data."		
sys.exit( 0 )

# hts data
hts_names_L = list()
hts_train_D = dict()
hts_test_D = dict()
f = open( hts_fn )
for row in f:
	if row[0] in unwanted: continue
	L = row.strip().split( '\t' )
	if row[0] == 'g':
		hts_names_L = L[1:]
		continue
	expr_D = dict( zip( hts_names_L, L[1:] ) )
	
	g = L[0]
	hts_train_D[g] = { k: expr_D[k] for k in train_hts_L }
	
	if len( test_hts_L ) != 0:
		hts_test_D[g] = { k: expr_D[k] for k in test_hts_L }
f.close()
	
# gene to ps
g2ps_D = dict()
f = open( g2ps_fn )
for row in f:
	if row[0] == 'g':
		continue
	L = row.strip().split( '\t' )
	if L[1] == '.':
		continue
	g = L[0]; ps = int( L[1] )
	if g not in g2ps_D:
		g2ps_D[g] = set([ps])
	else:
		g2ps_D[g].add(ps)
f.close()

#for each gene in hts_train_D
for g in hts_train_D:
	# if there are no probesets do not bother
	try:
		pss = list( g2ps_D[g] )
	except KeyError:
		print >> sys.stderr, "Missing:",g
		continue
	TRAIN_HTS = ",".join( [ hts_train_D[g][k] for k in train_hts_L ] )
	TEST_HTS = ",".join( [ hts_test_D[k] for k in test_hts_L ] )
	
	#get the probe data
	PROBES = list()
	TRAIN_INTENSITIES = ''
	TEST_INTENSITIES = ''
	for ps in pss:
		try:
			probes_data = ma_train_D[ps]
		except KeyError:
			print >> sys.stderr, "I don't yet have the probe data for %s" % ps
			continue
		PROBES += probes_data.keys()
		for p in probes_data:
			TRAIN_INTENSITIES += ",".join( [ ma_train_D[ps][p][k] for k in train_ma_L ] ) + ","
			TEST_INTENSITIES += ",".join( [ ma_test_D[ps][p][k] for k in test_ma_L ] ) + ","
	TRAIN_INTENSITIES = TRAIN_INTENSITIES.strip( ',' )
	TEST_INTENSITIES = TEST_INTENSITIES.strip( ',' )
	
	if TRAIN_HTS == "":
		TRAIN_HTS= 'NA'
	if TEST_HTS == "":
		TEST_HTS = 'NA'
	if TRAIN_INTENSITIES == "":
		TRAIN_INTENSITIES = 'NA'
	if TEST_INTENSITIES == "":
		TEST_INTENSITIES = 'NA'
	print "\t".join( map( str, [ g, len( train_ma_L ), len( test_ma_L ), len( PROBES), TRAIN_HTS, TEST_HTS, TRAIN_INTENSITIES, TEST_INTENSITIES ] ))
#print them together



	
