#!/home/paulk/software/bin/python
from __future__ import division
import sys
import cPickle
import pysam
import os.path
from key_functions import *

try:
	fn = sys.argv[1]
except IndexError:
	print >> sys.stderr,"Usage: %s <dabg_filtered_probe_intensities>" % os.path.basename(sys.argv[0])
	sys.exit(0)

unwanted = ['#','p','t']

probe_intensities = dict()
# for each probe get the intensities and store them in a dictionary
f = open(sys.argv[1])
for row in f:
	if row[0] in unwanted: continue	
	l = row.strip().split('\t')
	probe_intensities[int(l[0])] = map(int,l[7:])
f.close()

f = open("/home/paulk/MARS/2_intermediate_data/%s.pic" % fn[:-4],'w')
cPickle.dump(probe_intensities,f,cPickle.HIGHEST_PROTOCOL)
f.close()