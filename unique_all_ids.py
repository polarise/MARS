import sys
import random
f = open(sys.argv[1])
c = 0
g = dict()
for row in f:
	if c > 100: break
	gn,mp,en,eg = row.strip().split('\t')
	if mp not in g:
		g[mp] = {0:[gn],1:[en],2:[eg]}
	else:
		g[mp][0] += [gn]
		g[mp][1] += [en]
		g[mp][2] += [eg]
f.close()

# pick one mps at random
for mps in g:
	print "\t".join(map(str,[mps,random.choice(g[mps][0]),random.choice(g[mps][1]),random.choice(g[mps][2])]))
