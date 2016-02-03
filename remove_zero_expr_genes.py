import sys
f = open(sys.argv[1])
for row in f:
	l = row.strip().split('\t')
	if sum(map(float,l[4].split(','))) > 0:
		print row.strip()
f.close()
