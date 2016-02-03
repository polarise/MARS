#!/home/paulk/software/bin/python
from __future__ import division
import sys
import scipy

data = dict()
f = open(sys.argv[1])
for row in f:
	if row[0] == 'G': continue
	l = row.strip().split('\t')
	data[l[0]] = {'paul':map(float,l[1:])}
f.close()

f = open(sys.argv[2])
for row in f:
	l = row.strip().split('\t')
	if l[0] in data and l[1:].count('NA') == 0:
		data[l[0]]['paulk'] = map(float,l[1:])
	elif l[0] in data and l[1:].count('NA') > 0:
		del data[l[0]]
f.close()

#print data

for d in data:
	D = data[d]
	if 'paul' in D and 'paulk' in D: print "\t".join(map(str,[d,scipy.mean(D['paul']),scipy.mean(D['paulk']),scipy.median(D['paul']),scipy.median(D['paulk'])]))
