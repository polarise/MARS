#!/home/paulk/software/bin/python
from scipy import *
import sys
import rpy2.robjects as R
import rpy2.rinterface

t_test = R.r['t.test']

try:
	pred_fn = sys.argv[1]
	test_fn = sys.argv[2]
except IndexError:
	print >> sys.stderr,"Usage: script.py <predicted> <test>"
	sys.exit(1)

unwanted = ['S','g']

data = dict()
f = open(pred_fn)
for row in f:
	if row[0] in unwanted: continue
	l = row.strip().split('\t')
	data[l[0]] = {'pred':map(float,l[1:])}
f.close()

f = open(test_fn)
for row in f:
	if row[0] in unwanted: continue
	l = row.strip().split('\t')
	data[l[0]]['act'] = map(float,l[1:])
f.close()

# now do some t-tests
for g in data:
	pred = R.FloatVector(data[g]['pred'])
	act = R.FloatVector(data[g]['act'])
	try:
		results = t_test(pred,act)
	except rpy2.rinterface.RRuntimeError:
		results = [["NA"],[],["NA"]]
	print "\t".join(map(str,[g,results[0][0],results[2][0]]))
