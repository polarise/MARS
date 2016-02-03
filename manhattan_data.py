#!/usr/bin/env python
import gzip

# incremental lengths
lengths_D = {1:249250621,2:243199373,3:198022430,4:191154276,5:180915260,6:171115067,7:159138663,8:146364022,9:141213431,10:135534747,11:135006516,12:133851895,13:115169878,14:107349540,15:102531392,16:90354753,17:81195210,18:78077248,19:59128983,20:63025520,21:48129895,22:51304566}

X_len_i = 155270560

def get_lengths_before_F( chrom_i ):
	add_length_i = 0
	for i in range(1,23):
		if i < chrom_i:
			add_length_i += lengths_D[i]
		else:
			break
	return add_length_i
			
all_lengths_i = get_lengths_before_F( 23 )
add_len_X = all_lengths_i
add_len_Y	= all_lengths_i + X_len_i

# get snp coordinates
for chrom_i in range(1,23) + ['X', 'Y']:
#	chrom_i = 12
	if chrom_i == 'X':
		add_len_i = add_len_X
	elif chrom_i == 'Y':
		add_len_i = add_len_Y
	else:
		add_len_i = get_lengths_before_F( chrom_i )

	f = gzip.open( "/home/paulk/RP/resources/variation/snp_coords.chr%s.txt.gz" % chrom_i )
	coords_D = dict()
	for row_s in f:
		if row_s[:3] == 'rs#': continue
		L = row_s.strip().split('\t')
		coords_D[L[0]] = int( L[1] )
	f.close()

	# get the p-values
	root_s = "/home/paulk/MARS/4_association_tests/seekarray_mars/chr%s/association_results.txt.qval"
	f = open( root_s % chrom_i )
	for row_s in f:
		if row_s[0] == 's': continue
		L = row_s.strip().split('\t')
		pos_i = coords_D[L[0]]
		print "%s\t%s" % (pos_i + add_len_i,L[-1])
	f.close()

