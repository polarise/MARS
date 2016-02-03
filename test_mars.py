#!/home/paulk/software/bin/python
from __future__ import division
import sys
import scipy
import numpy
import scipy.linalg
numpy.linalg.lstsq = scipy.linalg.lstsq
import time
import itertools
import Orange

# prepare the predictor names
I = itertools.product('abcdefghijklmnopqrstuvwxyzABCDEF',repeat=2)
J = ["".join(i) for i in I]

f = open(sys.argv[1])
row = f.readline()
f.close()

g,no_samps,no_probes,probe_IDs,fpkm,probes = row.strip().split('\t')

no_samps = int(no_samps)
no_probes = int(no_probes)
fpkm = map(float,fpkm.split(','))
probes = map(float,probes.split(','))

# required Orange constructs
predictors = [Orange.feature.Continuous(name) for name in J[:no_probes]]
response = Orange.feature.Continuous('RS')
domain = Orange.data.Domain(predictors + [response])

probes_matrix = [[probes[no_samps*j + i] for j in xrange(no_probes)] + [fpkm[i]] for i in xrange(no_samps)]
data = Orange.data.Table(domain,probes_matrix)

M = Orange.regression.earth.EarthLearner(data,degree=2)

f = open(sys.argv[2])
row = f.readline()
f.close()

g,no_samps,no_probes,probe_IDs,fpkm,probes = row.strip().split('\t')

no_samps = int(no_samps)
no_probes = int(no_probes)
fpkm = map(float,fpkm.split(','))
probes = map(float,probes.split(','))

predictors = [Orange.feature.Continuous(name) for name in J[:no_probes]]
domain2 = Orange.data.Domain(predictors)

probes_matrix = [[probes[no_samps*j + i] for j in xrange(no_probes)] for i in xrange(no_samps)]
instances = [Orange.data.Instance(domain2,p) for p in probes_matrix]
print map(lambda x:M.predict(x)[0],instances)
#print M
#print probes_matrix[0],len(probes_matrix[0])
#print M.predict(Orange.data.Instance(domain2,probes_matrix[0]))
#print M
#print probes_matrix[1],len(probes_matrix[1])
#print M.predict(Orange.data.Instance(domain2,probes_matrix[1]))

#for preds in probes_matrix:
#	print preds
#	print M
#	print M.predict(Orange.data.Instance(domain2,preds))
print fpkm
print
