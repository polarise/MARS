#!/usr/bin/env python
from __future__ import division
import sys
import pysam
from key_functions import *
import cPickle
import time

"""
Notes:
This script was compiled from pipeline2.py. It performs almost exactly what pipeline2.py
did with the exception that the resource files used are trimmed down substantially to 
speed up computation.

Input files:
============
1. Exon coordinates - PIC file 
	(/data2/paulk/RP/resources/exons.hg19.Ens66.1-based.pic)
2. Probeset coordinates - indexed, zipped BED file 
	(/home/paulk/RP/resources/HuEx-1_0-st-v2.na31.hg19.probeset.bed.gz)
3. Probeset-probe intensity file - indexed, zipped BED file (ad hoc) 
	(YRI: /home/paulk/MARS/2_intermediate_data/0_probe_data/YRI_intensities_core_quantnorm_pmgcbg.txt.bed.gz ;
	 CEU: /home/paulk/MARS/2_intermediate_data/0_probe_data/CEU_intensities_core_quantnorm_pmgcbg.txt.bed.gz )
4. Transcript expression - PIC file 
	(YRI:/data1/paulk/%s/isoforms.fpkm_tracking ; 
	 CEU:/data1/paulk/Montgomery_Data/output/%s/isoforms.fpkm_tracking )

Output format:
==============
tx_id | no_samples | no_probes | probe_ids | tx_expr | probe_intensities

"""

unwanted = [ '#', 'p', 't' ]

try:
	label = sys.argv[1]
except IndexError:
	print >> sys.stderr,"Error!"
	sys.exit(1)

# names for ordering
if label == 'YRI':
	samples = map(lambda(x):x.strip(), open("/home/paulk/MARS/0_prep_data/YRI_names.txt").readlines())
elif label == 'CEU':
	samples = map(lambda(x):x.strip(), open("/home/paulk/MARS/0_prep_data/CEU_names.txt").readlines())

# number of samples
no_samples = len( samples )

# get the RNA-seq data
rna_seq_data = dict()

"""
# for each of the rna-seq data files
# we will use the coordinates given to ensure that we treat each transcript even it occurs multiple times
if label == 'YRI':
	for s in samples:
		f = open("/data1/paulk/%s/isoforms.fpkm_tracking" % s)
		for row in f:
			if row[0] in unwanted: continue
			l = row.strip().split('\t')
			feat = l[0]
			if feat not in rna_seq_data:
				rna_seq_data[feat] = [float(l[9])]    # we're using the coordinates (l[6]) to make each key unique
			else:
				rna_seq_data[feat] += [float(l[9])]
		f.close()
	
	f = open( "YRI_txs_expr.pic", 'w' )
	cPickle.dump( rna_seq_data, f, cPickle.HIGHEST_PROTOCOL )
	f.close()
	
elif label == 'CEU':
	for s in samples:
		f = open("/data1/paulk/Montgomery_Data/output/%s/tophat_out/isoforms.fpkm_tracking" % s)
		for row in f:
			if row[0] in unwanted: continue
			l = row.strip().split('\t')
			feat = l[0]
			if feat not in rna_seq_data:
				rna_seq_data[feat] = [float(l[9])]    # we're using the coordinates (l[6]) to make each key unique
			else:
				rna_seq_data[feat] += [float(l[9])]
		f.close()

	f = open( "CEU_txs_expr.pic", 'w' )
	cPickle.dump( rna_seq_data, f, cPickle.HIGHEST_PROTOCOL )
	f.close()
sys.exit(0)

"""

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

start = time.time()
# read in the exon coordinates
f = open( "/data2/paulk/RP/resources/exons.hg19.Ens66.1-based.pic" )
exons = cPickle.load( f )
f.close()
print >> sys.stderr, "Took %.5f seconds to read exons file" % ( time.time() - start )

# tabixes
probesets_tabix = pysam.Tabixfile( "/data2/paulk/RP/resources/HuEx-1_0-st-v2.na31.hg19.probeset.bed.gz")
if label == 'YRI':
	probes_tabix = pysam.Tabixfile( "/home/paulk/MARS/2_intermediate_data/0_probe_data/YRI_intensities_core_quantnorm_pmgcbg.txt.bed.gz" )
elif label == 'CEU':
	probes_tabix = pysam.Tabixfile( "/home/paulk/MARS/2_intermediate_data/0_probe_data/CEU_intensities_core_quantnorm_pmgcbg.txt.bed.gz" )

chroms = [ "chr"+x for x in map( str, range(1,23) ) + [ 'X', 'Y', 'MT' ]]

data = dict()
c = 0
for exon in exons:
	if c > 10000: break
	tx,ex = exon.split( ":" )
	
	# check to make sure that this is not one of those patch chromosomes e.g. chrHSCHR6_MHC_MCF
	chrom = exons[exon].split( ":" )[0]
	if chrom not in chroms: continue
	
	# make the region string
	region = exons[exon][:-2]
	
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
		data[ exon ] = [ ",".join( all_probes ), ",".join( all_probes_exprs )]
	c += 0

transcripts = dict()
for exon in data:
	tx,ex = exon.split( ":" )
	
	# get rnaseq
	try:
		rnaseq = rna_seq_data[tx]
	except KeyError:
		print >> stderr, "Transcripts '%s' missing RNA-seq values." % tx
		continue
	
	if tx not in transcripts:
		transcripts[tx] = [ rnaseq, data[exon][0], data[exon][1] ] # add this present exon and the rnaseq
	else:
		transcripts[tx][1] += "," + data[exon][0]	# only add this present exon's probes
		transcripts[tx][2] += "," + data[exon][1]	# intensities	

for tx in transcripts:
	no_probes = len( transcripts[tx][1].split( ',' ))
	print "\t".join( [tx, str( no_samples ), str( no_probes ), transcripts[tx][1], ",".join( map( str, transcripts[tx][0] )), transcripts[tx][2] ])
