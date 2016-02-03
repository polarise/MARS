#!/usr/bin/env python
import sys
import random
import scipy
import scipy.stats

# generate the gene ids
no_genes = 100
gene_ids = range( 1, no_genes+1 )

# determine the number of txs
no_txs = list()
for i in xrange( no_genes ):
	the_no = scipy.stats.poisson.rvs( 5 )
	if the_no > 0:
		no_txs.append( the_no )
	else:
		no_txs.append( 1 )

# generate the tx ids
tx_names = list()
for no in no_txs:
	tx_names.append( ["".join( [random.choice('abcdefghijklmnopqrstuvwxyz') for i in xrange(4)] ) for j in xrange(no) ] )
	
# no of samples
no_samples_train = 100
no_samples_test = 50
# no of probes
no_probes = [ scipy.stats.poisson.rvs( 100 ) for i in xrange( no_genes )]

# actual probes [optional]
probes = [ "-" for i in xrange( no_genes )]

# rnaseq tx estimates
rs_values_train = list()
rs_values_test = list()
for no in no_txs:
	rs_values_train.append( list( scipy.stats.mielke.rvs( 0.9, 0.9, size=no_samples_train*no )))
	rs_values_test.append( list( scipy.stats.mielke.rvs( 0.9, 0.9, size=no_samples_test*no )))

# probe intensities
probe_values_train = list()
probe_values_test = list()
for no in no_probes:
	probe_values_train.append( list( abs( scipy.stats.norm.rvs( 100, 20, size=no_samples_train*no ))))
	probe_values_test.append( list( abs( scipy.stats.norm.rvs( 100, 20, size=no_samples_test*no ))))

# put it all together
tr = open( "/home/paulk/MARS/7_tx_disambig/0_simulation/train_100_samples.txt", 'w' )
te = open( "/home/paulk/MARS/7_tx_disambig/0_simulation/test_50_samples.txt", 'w' )
for i in xrange( no_genes ):
	print >> tr, "\t".join( map( str, [ gene_ids[i], no_samples_train, no_probes[i], probes[i], ",".join( map( str, rs_values_train[i] )), ",".join( map( str, probe_values_train[i] )) ]))
	print >> te, "\t".join( map( str, [ gene_ids[i], no_samples_test, no_probes[i], probes[i], ",".join( map( str, rs_values_test[i] )), ",".join( map( str, probe_values_test[i] )) ]))
tr.close()
te.close()
