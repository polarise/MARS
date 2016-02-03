from __future__ import division
import operator
import scipy

probes = map(float,l[5].split(','))
probes_matrix = [[probes[52*j + i] for j in xrange(51)] for i in xrange(52)]
R = probes_matrix[:][2]
r_min = min(R)
r_max = max(R)
w = 10.
k = (r_max - r_min)/w
rint = []
i = 0
while i < k:
  rint += [(r_min + i*w,r_min + (i+1)*w)]
  i += 1

bins = {ri:0 for ri in rint}

for r in R:
  done = False
  for ri in bins:
    if ri[0] <= r < ri[1]:
      bins[ri] += 1
      done = True
    if done:
      break

for b in sorted(bins.keys()):
  print b,"\t",bins[b]

rep = max(bins.iteritems(),key=operator.itemgetter(1))
print rep
print scipy.mean(rep[0])
print scipy.median(R)



