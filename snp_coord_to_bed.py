#!/home/paulk/software/bin/python
"""
given a genotypes file outputs the coordinates as BED
BED format: chr|st|sp|name|NULL|sd
"""
import sys
import gzip

f = gzip.open( sys.argv[1] )
for row_s in f:
	if row_s[:3] == 'rs#': continue
	name_s, alleles_s, chrom_s, st_s, sd_s = row_s.strip().split( ' ' )[:5]
	print "\t".join( [ chrom_s, st_s, st_s, name_s, "-", sd_s ] )
f.close()
