#!/home/paulk/software/bin/python
from __future__ import division
import scipy
import sys

f = open(sys.argv[1])
for row in f:
	l = row.strip().split('\t')
	r = map(float,l[4].split(','))
	print "\t".join(map(str,[l[0],scipy.mean(r),scipy.median(r),scipy.var(r)]))
f.close()
