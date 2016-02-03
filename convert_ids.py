import sys

# read in ids
f = open(sys.argv[1])
ids = dict()
for row in f:
	l = row.strip().split('\t')
	ids[int(l[0])] = [l[2],l[3]]
f.close()

f_ens = open("rma.summary.Ens.txt",'w')
f_entrez = open("rma.summary.Entrez.txt",'w')

unwanted = ['#']
f = open(sys.argv[2])
header = ''
for row in f:
	if row[0] in unwanted: continue
	if row[0] == 'p':
		header = row.strip()
		print >> f_ens,header
		print >> f_entrez,header
		continue
	l = row.strip().split('\t')
	id = int(l[0])
	exprs = map(float,l[1:])
	try:
		ens,entrez = ids[id]
	except KeyError:
		print >> sys.stderr,id
		continue
	print >> f_ens,"\t".join(map(str,[ens]+exprs))
	print >> f_entrez,"\t".join(map(str,[entrez]+exprs))
f.close

f_ens.close()
f_entrez.close()	
