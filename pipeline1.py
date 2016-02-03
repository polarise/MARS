#!/home/paulk/software/bin/python
from __future__ import division
from sys import argv,exit,stderr
import pysam
from cPickle import load
from key_functions import *

# get all exons for Ens66
#f = open("/data2/paulk/RP/resources/exons.hg19.Ens66.1-based.pic")

# get all the genes for Ens66
f = open("/data2/paulk/RP/resources/Homo_sapiens.GRCh37.66.gene.pic")
exons = load(f)
f.close()

# get the expressions by transcript
#f = open("/data2/paulk/MARS/transcripts.expr")
"""f = open("/data2/paulk/MARS/YRI_RNA_Seq.expr")
txs_expr = load(f)
f.close()
"""

# get the expressions by gene
f = open("/data2/paulk/MARS/genes.expr")
#f = open("/data2/paulk/MARS/YRI_RNA_Seq_genes.expr")
txs_expr = load(f)
f.close()

# open a tabix object for probesets
tabixfile = pysam.Tabixfile("/data2/paulk/RP/resources/HuEx-1_0-st-v2.na31.hg19.probeset.gtf.gz")

# open a tabix object for probe intensities
tabixfile2 = pysam.Tabixfile("/data2/paulk/MARS/probe-intensities.txt.gz")
#tabixfile2 = pysam.Tabixfile("/data2/paulk/MARS/data/YRI-probe-intensities.txt.gz")

combined_data = dict()
c = 0
#f = open("/data2/paulk/RP/resources/Homo_sapiens.GRCh37.66.exons.gene")
for e in exons:	# for each gene/exon
	if c >= 20: break
#	tx,exon = e.split(':')
	tx = e
	try:
		tx_expr = txs_expr[tx]	# get the RNA-Seq values
	except:
		tx_expr = ["NA"]*8
	region = exons[e]	# get the coordinates defining the gene's span
	try:
		results = tabixfile.fetch(region=region[:-2])
	except ValueError:
		continue
	for row in results:
		feature = process_feature(row)
		ps_id = int(feature['transcript_id'])
		results2 = tabixfile2.fetch(region="chrP:%s-%s" % (ps_id,ps_id+1))
		probe_values = get_probe_values(results2)
		if tx not in combined_data:
			combined_data[tx] = probe_values	# probe_values is a dictionary
			combined_data[tx].update({'rna-seq':tx_expr})	# append a k-v pair
		else:
			combined_data[tx].update(probe_values) 
	c += 0
	
c = 0
for e in combined_data:	# for each exon in {combined_data}
	probes = [ p for p in combined_data[e].keys() if p != 'rna-seq']
	print e + "\t" + "\t".join( map( str, [ len(combined_data[e]['rna-seq']), len(combined_data[e]) - 1])) + "\t" + ",".join( [ str(p) for p in probes ] ) + "\t" + ",".join( map( str, combined_data[e]['rna-seq'] )) + "\t",
	probe_vals = ""
	for p in probes:
		probe_vals += ",".join( map( str, combined_data[e][p] )) + ","
	print probe_vals.strip(',')
