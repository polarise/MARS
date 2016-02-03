#!/usr/bin/env python
"""
script to perform association tests for significant SNPs found with RNA-seq
"""
import sys

# generic path
path_to_rnaseq_s = "/home/paulk/MARS/4_association_tests/rnaseq_10k/chr%s/%s"

# files
for i in range(1,23) + ['X','Y']:
	chrom_file_s = path_to_rnaseq_s % ( i, "association_results.txt.qval" )
	genotypes_file_s = path_to_rnaseq_s % ( i, "converted_genotypes.txt" )
	expression_file_s = sys.argv[1]	## handle this better!

	# get the SNPs that are signficicant at a certain level together with their genes
	snps_genes_D = dict()
	f = open( chrom_file_s )
	for row_s in f:
		if row_s[0] == 's': continue
		L = row_s.strip().split( '\t' )
		# apply a filter: we only want those that show some measure of significance
		if float( L[-1] ) < 1:
			snps_genes_D[L[0]] = L[1:4]
		else:
			break
	f.close()	

	# get the genotype for that SNP
	genotypes_D = dict()
	f = open( genotypes_file_s )
	for row_s in f:
		L = row_s.strip().split( '\t' )
		genotypes_D[L[0]] = L[3],L[4]
	f.close()

	# get the expression (Affy RMA or Seekized)
	expr_D = dict()
	f = open( expression_file_s )
	for row_s in f:
		L = row_s.strip().split( '\t' )
		expr_D[L[0]] = L[1:]
	f.close()

	# make the 'all data' file
	# use only those SNPs from snps_genes_D
	for snp_s in snps_genes_D:
		gene_s, ma_c, alleles_s = snps_genes_D[ snp_s ]
		try:
			expr_s = expr_D[ gene_s ] # assume that all expressions exist
		except KeyError:
			print >> sys.stderr,"Warning: missing expression for gene: %s" % gene_s
		no_s, genotypes_s = genotypes_D[ snp_s ]
		print "\t".join( [snp_s, gene_s, ma_c, alleles_s, no_s, genotypes_s, ",".join( expr_s )] )
