import sys

# Note: pauls data is missing NA11831.CEL, NA11832.CEL, NA12154.CEL, NA12717.CEL

# read in the order of pauls data
f = open(sys.argv[1])	# pauls_header_column.txt
pauls_ids = list()
for row in f:
	pauls_ids.append(row.strip())
f.close()

# read in the new ids
f = open(sys.argv[2]) # combined_reduced_CEU_pauls.txt
new_ids = dict()
for row in f:
	l = row.strip().split('\t')
	new_ids[l[0]] = l[1]	
f.close()

# create a list of the new ids in the order of pauls data
new_ids_in_paul_order = list()
for p in pauls_ids:
	new_ids_in_paul_order.append(new_ids[p])
	
# read in the data and for each row output the data in pauls' order

f = open(sys.argv[3])	# ../2_intermediate_data/rma_CEU/rma.summary.Ens.txt OR ../2_intermediate_data/rma_CEU/rma.summary.Entrez.tx
for row in f:
	if row[0] == 'p':
		col_names = row.strip().split('\t')[1:]
		print "\t".join(["probeset_id"] + [n+".CEL" for n in new_ids_in_paul_order])
		continue
	l = row.strip().split('\t')
	
	# get by name
	for c in col_names:
		row_data = dict(zip(col_names,l[1:]))
		
	# now use the new ordering
	l_new = [row_data[k+".CEL"] for k in new_ids_in_paul_order]
#	print len(l_new)
	
	# print out the newly-ordered data
	print "\t".join(map(str,[l[0]] + l_new))
f.close()

