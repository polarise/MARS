import sys
import random
f = open(sys.argv[1])
c = 0
g = dict()
for row in f:
	if row[0] == '#': continue
	if c > 100: break
	l = row.strip().split(',')
	if row[:2] == '"t': continue
	gene_txt = l[7].strip('"')
	mps = int(l[0].strip('"'))
	if gene_txt != '---':
		gene_name = gene_txt.split(' /// ')[0].split(' // ')[1]
		if gene_name not in g:
			g[gene_name] = [mps]
		else:
			g[gene_name] += [mps]
	c += 0
f.close()

# pick one mps at random
for gene in g:
	print str(random.choice(g[gene]))+"\t"+gene
	
