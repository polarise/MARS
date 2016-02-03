#!/usr/bin/env python
import sys

d = dict()
f = open( sys.argv[1] ) # /home/paulk/MARS/9B_custom_mps/gene_to_HuEx_probesets_Ens72.core.unique.txt
for row in f:
	if row[0] in [ 'g', 'p' ]: continue
	L = row.strip().split( '\t' )
	if L[0] == '': continue
	if L[1] not in d:
		d[L[1]] = [L[0]]
	else:
		d[L[1]] += [L[0]]
f.close()

print "probeset_id\tprobeset_list"
for r in d:
	print "%s\t%s" % ( r, " ".join( set( d[r] )) )



