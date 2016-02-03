#!/home/paulk/software/bin/python
import sys
from scipy import *
import cPickle

def linear_to_matricular(row):
	l = row.strip().split('\t')
	no_samples = int(l[1])
	no_probes = int(l[2])
	probes = map(int,l[3].split(','))
	probe_expr = map(float,l[5].split(','))
	
	probe_matrix  = list()
	for i in xrange(no_probes):
		probe_matrix.append(probe_expr[i*no_samples:i*no_samples+no_samples])
		
	try:
		probe_array = array(probe_matrix,dtype=int)
	except ValueError:
		print >> sys.stderr,row.strip()
		sys.exit(1)

	return probe_array
	
def filter_matrix(probe_array,probes,no_probes,passed):
	probe_vectors = dict()
	passed_probes = list()
	for i in xrange(no_probes):
		if probes[i] not in passed:
			pass
		else:
			probe_vectors[probes[i]] = probe_array[i,:]	
			passed_probes.append(probes[i])
	
	filtered_matricular = list()
	for p in passed_probes:
		filtered_matricular.append(probe_vectors[p])
		
	filtered_matricular = array(filtered_matricular)
	return passed_probes,filtered_matricular
	
def matricular_to_linear(matricular):
	linear_expr = list()
	for row in matricular:
		for item in row:
			linear_expr.append(item)
	return linear_expr


#f = open("/home/paulk/MARS/data/YRI_intensities_dabg.txt/dabg_passed.txt")
#passed = {int(row.strip()):0 for row in f}
#f.close()

#f = open("/home/paulk/MARS/data/YRI_intensities_dabg.txt/dabg_passed.pic",'w')
#cPickle.dump(passed,f,cPickle.HIGHEST_PROTOCOL)
#f.close()

#f = open("/home/paulk/MARS/data/YRI_intensities_dabg.txt/dabg_not_passed.txt")
#not_passed = {int(row.strip()):0 for row in f}
#f.close()

#f = open("/home/paulk/MARS/data/YRI_intensities_dabg.txt/dabg_not_passed.pic",'w')
#cPickle.dump(not_passed,f,cPickle.HIGHEST_PROTOCOL)
#f.close()

#sys.exit(0)

f = open("/home/paulk/MARS/data/YRI_intensities_dabg.txt/dabg_passed.pic")
passed = cPickle.load(f)
f.close()

f = open("/home/paulk/MARS/data/YRI_intensities_dabg.txt/dabg_not_passed.pic")
not_passed = cPickle.load(f)
f.close()

f = open(sys.argv[1])
c = 0
for row in f:
	if c > 2: break	
	l = row.strip().split('\t')
	no_samples = int(l[1])
	no_probes = int(l[2])
	probes = map(int,l[3].split(','))
	probe_expr = map(float,l[5].split(','))
	probe_array = linear_to_matricular(row)
	
	passed_probes,filtered_matricular = filter_matrix(probe_array,probes,no_probes,passed)
	
#	print no_samples,no_probes
#	print shape(probe_array)
#	print passed_probes
#	print probe_array
	
#	for p in probe_vectors:
#		print p,probe_vectors[p]

#	print filtered_matricular
#	print shape(filtered_matricular)
	linear_expr = matricular_to_linear(filtered_matricular)
#	print linear_expr
	
#	print
	print "\t".join([l[0],l[1],str(len(passed_probes)),",".join(map(str,passed_probes)),l[4],",".join(map(str,linear_expr))])
#	print
	c += 0

f.close()
