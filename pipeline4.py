#!/usr/bin/env python
from __future__ import division
import sys
from key_functions import *
import cPickle
import time
import os.path
import pysam

try:
	label = sys.argv[1]
except IndexError:
	print >> sys.stderr,"Error!"
	sys.exit(1)

# map of gene to transcripts
"""
f = open( "/home/paulk/RP/resources/gene_to_txs.hg19.Ens66" )
gene_to_txs = dict()
for row_s in f:
	L = row_s.strip().split( '\t' )
	if L[0] not in gene_to_txs:
		gene_to_txs[L[0]] = [L[1]]
	else:
		gene_to_txs[L[0]] += [L[1]]
f.close()

f = open( "/home/paulk/RP/resources/gene_to_txs.hg19.Ens66.pic", 'w' )
cPickle.dump( gene_to_txs, f, cPickle.HIGHEST_PROTOCOL )
f.close()
"""
start = time.time()
f = open( "/home/paulk/RP/resources/gene_to_txs.hg19.Ens66.pic" )
gene_to_txs = cPickle.load( f )
f.close()
print >> sys.stderr, "Took %.5f seconds to obtain map of genes to transcripts" % ( time.time() - start )

# map of gene coordinates
start = time.time()
f = open( os.path.expanduser( "~/RP/resources/Homo_sapiens.GRCh37.66.gene.pic" ))
gene_coords = cPickle.load( f )
f.close()
print >> sys.stderr, "Took %.5f seconds to obtain gene coordinates" % ( time.time() - start )

# map of transcripts to transcript expression
start = time.time()
if label == 'YRI':
	f = open( "/home/paulk/MARS/2_intermediate_data/2_rnaseq_data/YRI_txs_expr.pic" )
	rna_seq_data = cPickle.load( f )
	f.close()
elif label == 'CEU':
	f = open( "/home/paulk/MARS/2_intermediate_data/2_rnaseq_data/CEU_txs_expr.pic" )
	rna_seq_data = cPickle.load( f )
	f.close()	
print >> sys.stderr, "Took %.5f seconds to obtain transcript expression estimates" % ( time.time() - start )

# tabixes
probesets_tabix = pysam.Tabixfile( "/data2/paulk/RP/resources/HuEx-1_0-st-v2.na31.hg19.probeset.bed.gz")

if label == 'YRI':
	probes_tabix = pysam.Tabixfile( "/home/paulk/MARS/2_intermediate_data/0_probe_data/YRI_intensities_core_quantnorm_pmgcbg.txt.bed.gz" )
elif label == 'CEU':
	probes_tabix = pysam.Tabixfile( "/home/paulk/MARS/2_intermediate_data/0_probe_data/CEU_intensities_core_quantnorm_pmgcbg.txt.bed.gz" )

# subset genes
# setset transcripts
start = time.time()
print >> sys.stderr, "Subsetting genes to only those with transcripts...",
subset_genes = dict()
for gene in gene_to_txs:
	subset_genes[gene] = list()
	for tx in gene_to_txs[gene]:
		if tx in rna_seq_data:
			subset_genes[gene] += [tx]
print >> sys.stderr, "done!"
print >> sys.stderr, "Took %.5f seconds to subset genes" % ( time.time() - start )

genes_to_delete = list()
for gene in subset_genes:
	if len( subset_genes[gene] ) == 0:
		genes_to_delete += [gene]
#		print >> sys.stderr, subset_genes[gene]
#		i += 1
#print >> sys.stderr, i
#sys.exit(0)

for gene in genes_to_delete:
	del subset_genes[gene]


chroms = [ "chr"+x for x in map( str, range(1,23) ) + [ 'X', 'Y', 'MT' ]]

start = time.time()
data = dict()
c = 0
# for each gene with an associated transcript expression (ignore those that don't have txs expression)
for gene in subset_genes:
	if c > 100: break
	# get the coordinates
	try:
		coords = gene_coords[gene]
	except KeyError:
		print >> sys.stderr, "What's this: %s?" % gene
		continue
	
	# check to make sure that this is not one of those patch chromosomes e.g. chrHSCHR6_MHC_MCF
	chrom = coords.split( ":" )[0]
	if chrom not in chroms: continue
	
	# make the region string
	region = coords[:-2]
	
	# get the probeset
	results = probesets_tabix.fetch( region=region )
	mps_ids, ps_ids = process_bed_rows( results )
	
	all_probes = list()
	all_probes_exprs = list()
	for ps_id in ps_ids:
		region2 = "chrP:%s-%s" % ( ps_id, int(ps_id) + 1 )

		# get the probe expression intensities
		results2 = probes_tabix.fetch( region = region2 )
		probes, probes_exprs = process_bed_rows( results2 )

		all_probes += probes
		all_probes_exprs += probes_exprs

	if len( all_probes ) > 0:
		data[ gene ] = [ subset_genes[gene], ",".join( all_probes ), ",".join( all_probes_exprs )]
	c += 0
print >> sys.stderr, "Took %.5f seconds to get data for %d genes." % ( time.time() - start, c )

for d in data:
	txs, probes, probes_expr = data[d]
	no_txs = len( txs )
	no_probes = len( probes.split( ',' ))
	txs_expr = list()
	for tx in txs:
		txs_expr += rna_seq_data[tx]
	try:
		no_samples = int( len( txs_expr )/no_txs )
	except KeyError:
		print >> sys.stderr, "A gene (%s) with no transcripts! How?" % d
		continue
	print "\t".join( [ d, ",".join( txs ), str( no_txs ), str( no_samples ), str( no_probes ), probes, ",".join( map( str, txs_expr )), probes_expr ])
		

