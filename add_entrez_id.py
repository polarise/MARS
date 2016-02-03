#!/home/paulk/software/bin/python
import sys

f = open("/data2/paulk/RP/resources/ENSEMBL_to_Entrez_proper.txt")
ENS2Entr = dict()
for row in f:
	l = row.strip().split('\t')
	ENS2Entr[l[0]] = l[1]
f.close()

f = open(sys.argv[1])
for row in f:
	id1 = row.strip().split('\t')[0].split('_')[0]
	try:
		id2 = ENS2Entr[id1]
	except KeyError:
		id2 = "NA"
	print row.strip()+"\t"+id2
f.close()
