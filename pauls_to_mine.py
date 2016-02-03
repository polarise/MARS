import sys
import itertools

f = open("rnaseq%s.txt" % sys.argv[1])
g = open("exonarray%s.txt" % sys.argv[1])
for frow,grow in itertools.izip(f,g):
	R = frow.strip().split(',')
	E = grow.strip().split(',')
	no = len(R[1:])
	if len(E) == 1:
		
	else:
		print "\t".join([R[0],str(no),str(len(E[1:])/no),"-",",".join(R[1:]),",".join(E[1:])])	
f.close()
g.close()
