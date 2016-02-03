#!/usr/bin/env python
import random

f = open("/home/paulk/MARS/2_intermediate_data/rma_CEU/rma.summary.Entrez.pauls.txt")
data = dict()
for row in f:
	if row[0] == 'p': print row.strip(); continue
	l = row.strip().split('\t')
	if l[0] not in data:
		data[l[0]] = [l[1:]]
	else:
		data[l[0]] += [l[1:]]
f.close()

for d in data:
	print "\t".join( [d] + random.choice( data[d] ))
