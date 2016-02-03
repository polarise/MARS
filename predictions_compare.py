import sys
import rpy2.robjects as R

cor_test = R.r['cor.test']
options = R.r['options']
options(warn=-1)

f = open(sys.argv[1])
data = dict()
for row in f:
	if row[0] == 'G': continue
	l = row.strip().split('\t')
	data[l[0]] = {'paul':map(float,l[1:])}
f.close()

f = open(sys.argv[2])
for row in f:
	l = row.strip().split('\t')
	if l[0] in data:
		if l[1].find('NA') < 0 and sum(map(float,l[1].split(','))) > 0:
			data[l[0]]['paulk'] = map(float,l[1].split(','))
		elif l[1].find('NA') >= 0:
			del data[l[0]]
	else:
		pass
f.close()

c = 0
for d in data:
	if c > 5: break
	D = data[d]
	if 'paul' in D and 'paulk' in D and D['paul'][0]:
#		print d,D
		paul = R.FloatVector(D['paul'])
		paulk = R.FloatVector(D['paulk'])
		result = cor_test(paul,paulk,method="pearson")
		if -1 <= result[3][0] <= 1: print "\t".join(map(str,[d,result[0][0],result[2][0],result[3][0]]))
		c += 0
