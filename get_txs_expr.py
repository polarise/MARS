#!/home/paulk/software/bin/python
from __future__ import division
from sys import argv,exit,stderr
from cPickle import dump,HIGHEST_PROTOCOL

unwanted = ['#','t']

samples = map(lambda(x):x.strip(),open("YRI_names.txt").readlines())

txs = dict()
for s in xrange(3,11):
#for s in samples:
	f = open("/data2/paulk/RP/output2/%sC/Genic/genes.fpkm_tracking" % s)
#	f = open("/data1/paulk/%s/isoforms.fpkm_tracking" % s)
#	f = open("/data1/paulk/%s/genes.fpkm_tracking" % s)
	for row in f:
		if row[0] in unwanted: continue
		l = row.strip().split('\t')
		tx = l[0]
		if tx not in txs:
			txs[tx,l[6]] = [float(l[9])]
		else:
			txs[tx,l[6]] += [float(l[9])]
	f.close()
	
#f = open("/home/paulk/MARS/2_intermediate_data/RP_txs_expr.pic",'w')
f = open("/home/paulk/MARS/2_intermediate_data/RP_gene_expr.pic",'w')
#f = open("/home/paulk/MARS/2_intermediate_data/YRI_gene_expr.pic",'w')
dump(txs,f,HIGHEST_PROTOCOL)
f.close()
