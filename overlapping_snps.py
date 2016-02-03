#!/home/paulk/software/bin/python
"""
given BED files for SNPs and genes gets SNPs that overlap each gene
"""
from __future__ import division
import sys
import gzip 	# for gzip.open since the bed files are now zipped
import pysam	# for pysam.tabix

def process_genes_F( genes_I ):
	"""
	given fetch results extracts the genes
	"""
	genes_found_L = [ gene_s.split( '\t' )[3] for gene_s in genes_I ]
	return genes_found_L,len(genes_found_L)

try:
	snpfn_s, genefn_s = sys.argv[1:3]
except ValueError:
	print >> sys.stderr,"Usage. ./script.py <snpfn_s> <genefn_s>"
	sys.exit(1)
	
snp_p = gzip.open(snpfn_s)
genetabix_p = pysam.Tabixfile(genefn_s)

c_i = 0
for row_s in snp_p:
	if c_i > 10: break
	L = row_s.strip().split( '\t' )
	
	# make a region string
	snp_region_s = L[0] + ":" + L[1] + "-" + L[2]
	genes_I = genetabix_p.fetch( region = snp_region_s )
	
	genes_found_L, no_genes_i = process_genes_F( genes_I )
	
	if not no_genes_i:
		print >> sys.stderr,"Missing genes: %s" % L[3]
		continue
	
	for gf_s in genes_found_L:
		print "\t".join( [ L[3], gf_s ] )
	c_i += 0
snp_p.close()
