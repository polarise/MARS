#!/home/paulk/software/bin/python

f = open("combined_reduced_YRI.txt")
for row in f:
	l = row.strip().split('\t')
	print "ln -s /data2/paul/rawExonArrayCels/%s %s.CEL" % (l[3],l[0])
f.close()
