from __future__ import division
import sys

pass_dabg = 0
total = 0
f =  open(sys.argv[1])
for row in f:
  if row [0] in ['#','p']: continue
  l = row.strip().split('\t')
  dabgs = map(float,l[5:])
  pass_dabg += sum([1 for d in dabgs if d <= 0.05])
  total += len(dabgs)
f.close()

print pass_dabg
print total
print pass_dabg/total
  

