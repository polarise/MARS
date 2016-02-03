#!/home/paulk/software/bin/python
from __future__ import division
from sys import argv,exit,stderr
from cPickle import dump,HIGHEST_PROTOCOL

unwanted = ['#','p']

ps2p = dict()
probe_expr = dict()
f = open(argv[1])
for row in f:
	if row[0] in unwanted: continue
	l = row.strip().split('\t')
	ps_id = int(l[4])
	p_id = int(l[0])
	probe_expr[p_id] = map(int,l[7:])
	if ps_id not in ps2p:
		ps2p[ps_id] = [p_id]
	else:
		ps2p[ps_id] += [p_id]
f.close()

f = open("probeset_to_probes.pic",'w')
dump(ps2p,f,HIGHEST_PROTOCOL)
f.close()

f = open("probe_expression.pic",'w')
dump(probe_expr,f,HIGHEST_PROTOCOL)
f.close()
