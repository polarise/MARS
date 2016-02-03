import sys
import csv

f = open( sys.argv[1] )
f_rdr = csv.reader( f, delimiter=",", quotechar='"' )
c = 0
for row in f_rdr:
	if c > 1000: break
	if row[0][0] == '#': continue
	if row[0][0] == 't': continue
	if row[7] == '---': continue
	gene_names = set( [ k.split( ' // ' )[1] for k in row[7].split( ' /// ') ] )
	tx_ids = set( [ k.split( ' // ' )[0] for k in row[7].split( ' /// ') ] )
	if len( gene_names ) > 1: continue
	print "%s\t%s" % ( row[0], list( gene_names )[0] )# , list( tx_ids )
	c += 0
f.close()
