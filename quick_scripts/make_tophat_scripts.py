#!/home/paulk/software/bin/python

f = open("data/common_YRI.txt")
data = dict()
for row in f:
	l = row.strip().split('\t')
	if l[0] not in data:
		data[l[0]] = [l[1:6],[l[6]]]
	else:
		data[l[0]][1] += [l[6]]
f.close

for r in data:
	print "\t".join([r]+data[r][0]+[",".join(data[r][1])])
