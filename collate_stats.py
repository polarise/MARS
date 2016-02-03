#!/home/paulk/software/bin/python
from __future__ import division
import sys
import scipy

f = open(sys.argv[1])
for row in f:
	l = row.strip().split('\t')
	true_fpkm = map(float,l[1].split(','))
	if l[2].find('NA') < 0:
		pred_fpkm = map(float,l[2].split(','))
	else:
		pred_fpkm = None
	if l[3] == 'nan':
		l[3] == 'NA'
	if pred_fpkm != None: print "\t".join(map(str,[l[0],scipy.mean(true_fpkm),scipy.std(true_fpkm),scipy.mean(pred_fpkm),scipy.std(pred_fpkm),l[3],l[4]]))
	else: print "\t".join(map(str,[l[0],scipy.mean(true_fpkm),scipy.std(true_fpkm),"NA","NA",l[3],l[4]]))
f.close()
