#!/home/paulk/software/bin/python
from __future__ import division
import sys
import rpy2.robjects as R

median = R.r['median']
wilcox_test = R.r['wilcox.test']
options = R.r['options']

options(warn=-1)

f = open(sys.argv[1])
c = 0
for row in f:
	if c > 10: break
	l = row.strip().split('\t')
	if l[2].find('NA') < 0:
		true = R.FloatVector(map(float,l[1].split(',')))
		pred = R.FloatVector(map(float,l[2].split(',')))
		result = wilcox_test(true,pred)
		print "\t".join(map(str,[median(true)[0],median(pred)[0],result[2][0]]))
		c += 0
f.close()
