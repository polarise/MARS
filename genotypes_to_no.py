#!/home/paulk/software/bin/python
"""
given a genotypes file, order of genotypes and MAF_min outputs the genotypes as numbers
"""
from __future__ import division
import sys
import gzip
import operator	# will use operator.itemgetter

def genotypes_to_no_F( l_L, alleles_L, maf_f ):
	"""
	given a list of genotypes, the set of alleles and a MAF:
		- determines if MAF is met
		- converts the genotypes to a count of the minor alleles
	"""
	# ignore non-SNPs
	if "-" in alleles_L or len(alleles_L[0]) > 1 or len(alleles_L[1]) > 1:
		return '', '', 0.0, '0', False
	
	# get allele counts
	all_alleles_i = len( l_L )*2
	count_allele0_i = [ l_s.count( alleles_L[0] ) for l_s in l_L ]
	count_allele1_i = [ l_s.count( alleles_L[1] ) for l_s in l_L]
	pc_allele0_f = sum( count_allele0_i) / all_alleles_i*100
	pc_allele1_f = sum( count_allele1_i) / all_alleles_i*100
	
	# pass or fail MAF?
	if pc_allele0_f < maf_f or pc_allele1_f < maf_f:
		return '', '', 0.0, '0', False
	
	# determine the minor allele
	if pc_allele0_f < pc_allele1_f:
		minor_allele_s = alleles_L[0]
		count_allele0_s = ",".join( map( str, count_allele0_i ) )
		return count_allele0_s, minor_allele_s, pc_allele0_f, str(len(count_allele0_i)), True
	else:
		minor_allele_s = alleles_L[1]
		count_allele1_s = ",".join( map( str, count_allele1_i ) )
		return count_allele1_s, minor_allele_s, pc_allele1_f, str(len(count_allele1_i)), True	
	
try:
	genotypesfn_s = sys.argv[1]
	orderingfn_s = sys.argv[2]
	maf_f = float( sys.argv[3] )
except IndexError:
	print >> sys.stderr,"Usage: ./script.py <genotypesfn_s> <orderingfn_s> <MAF_f>"
	sys.exit(1)

# get the individuals in order
f = open( orderingfn_s )
ordering_L = [ row_s.strip().split( '\t' )[0] for row_s in f ]
f.close()

f = open( orderingfn_s )
population_L = [ '0' if row_s.strip().split( '\t' )[4] == 'Yoruban' else '1' for row in f ]
f.close()

population_s = ",".join( population_L )

# open the genotypes file
f = gzip.open( genotypesfn_s )
for row_s in f:
	L = row_s.strip().split( ' ' )
	if row_s[:3] == 'rs#':
		unordered_L = L[11:] # get the 'unordered' 
		continue		
	snp_s = L[0]
	alleles_s = L[1]
	alleles_L = L[1].split( '/' )
	genotypes_L = L[11:]
	genotypes_D = dict( zip( unordered_L,genotypes_L ) )
	ordered_genotypes_L = list()
	for o_s in ordering_L: # according to the ordering provided
		try:
			genotype_s = genotypes_D[o_s]
		except KeyError:
			print >> sys.stderr,"Missing %s" % o_s
			continue
		ordered_genotypes_L.append( genotype_s )
	ordered_no_s, minor_allele_s, pc_minor_allele_f, no_ind_s, status_b = genotypes_to_no_F( ordered_genotypes_L, alleles_L, maf_f )
	
	if status_b:
		print "\t".join( [ snp_s, minor_allele_s, alleles_s, no_ind_s, ordered_no_s, population_s ] )
	else:
		print >> sys.stderr,"%s failed to pass MAF = %.5f: genotype=%s @ %.5f" % ( snp_s, maf_f, minor_allele_s, pc_minor_allele_f )
		print >> sys.stderr,row_s
f.close()
