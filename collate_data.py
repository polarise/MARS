#!/usr/bin/env python
from __future__ import division
import sys
import cPickle
import argparse
import time

"""
def cnow():
	\"""
	Calendar time now.
	\"""
	return time.ctime( time.time() )
	
parser = argparse.ArgumentParser( description="Script to compile data file for \
training or testing." )
parser.add_argument( '-d', '--data-path', help="a string (delimited by \"\") ha\
ving a single asterisk indicating the folder containing quantified RNA-Seq resu\
lts (using Cufflinks)" )
parser.add_argument( '-s', '--samples-filter', help="name of file containing sa\
mple names to replace the asterisk" )
parser.add_argument( '-f', '--feature', default='gene', help="(g)ene | (t)ransc\
ripts; also possible 'txs' for transcripts [default: gene]" )
parser.add_argument( '-o', '--output', help="where to write output to [default:\
 stdout]" )
parser.add_argument( '-n', '--feature-name', help="name of file containing quan\
tified RNA-Seq data (using Cufflinks) [default: 'genes.fpkm_tracking']" )
parser.add_argument( '-l', '--library', help="name of file containing mappings \
from genes/transcripts to individual probes" )
parser.add_argument( '-p', '--probes', help="name of file containing quantile-n\
ormalised (quantnorm) and background-corrected (pm-gcbg) raw probe intensities")

args = parser.parse_args()

data_path = args.data_path
samples_filter = args.samples_filter
feature = args.feature
output = args.output
library = args.library
probes = args.probes
"""

unwanted = [ '#', 't', 'C', 'p' ]

class Region(object):
	def __init__( self, region_s ):
		self.chrom,st_sp = region_s.split( ':' )
		self.start, self.end = map( int, st_sp.split( '-' ) )
		self.length = self.end - self.start + 1
	
	def string( self ):
		return self.chrom + ":" + str(self.start) + "-" + str(self.end)

def find_max( region_list ):
	"""
	given a list of region strings returns the longest (assuming they are all on the same chromosome - unlikely to have different chromosomes)
	"""
	region_objects = map( Region, region_list )
	max_len = region_objects[0].length
	max_reg = region_list[0]
	for R in region_objects:
		if R.length > max_len:
			max_len = R.length
			max_reg = R.string()
	return max_reg	

try:
	probe_fn = sys.argv[1]
	hts_fn = sys.argv[2]
except IndexError:
	print >> sys.stderr, "Usage: ./script <probe_fn> <hts_fn>"
	sys.exit( 1 )

# genes to probesets
f = open( "/data2/paulk/MARS/0_prep_data/genes_to_probesets.bed" )
g2ps = dict()
for row in f:
	L = row.strip().split( '\t' )
	if L[7] == '-1': continue
	if L[3] not in g2ps:
		g2ps[L[3]] = [int(L[9])]
	else:
		g2ps[L[3]] += [int(L[9])]
f.close()

# probesets to probes and intensities
#f = open( "2_intermediate_data/0_probe_data/CEU_intensities_core_quantnorm_pmgcbg.txt" )
#f = open( "2_intermediate_data/0_probe_data/YRI_intensities_core_quantnorm_pmgcbg.txt" )
no_samples = 0
f = open( probe_fn )
ps2int = dict()
for row in f:
	if row[0] in unwanted: continue
	L = row.strip().split( '\t' )
	if row[0] == 'p':
		sample_names = map( lambda x: x[:7], L[7:] )
		continue
	ps = int(L[4])
	p = L[0]
	if no_samples == 0:
		no_samples = len( L[7:] )
	
	ints = ",".join( L[7:] )
	
	if ps not in ps2int:
		ps2int[ps] = dict()
	
	ps2int[ps][p] = ints

f.close()

#no_samples = len( sample_names )

rs_data = dict()
f = open( hts_fn )
for row in f:
	L = row.strip().split( '\t' )
	rs_data[ L[0] ] = ",".join( L[1:] )
f.close()	

"""
for sn in sample_names:
	f = open( "/data1/paulk/%s/genes.fpkm_tracking" % sn )
#	f = open( "/data1/paulk/Montgomery_Data/output/%s/tophat_out/genes.fpkm_tracking" % sn )
	for row in f:
		if row[0] in unwanted: continue
		L = row.strip().split( '\t' )
		feat = L[0]
		coord = L[6]
		
		if feat not in rs_data:
			rs_data[feat] = dict()
		
		if coord not in rs_data[feat]:
			rs_data[feat][coord] = L[9]
		else:
			rs_data[feat][coord] += "," + L[9]
	f.close()

# filter RNA-Seq data to have only unique genes
for feat in rs_data.keys():
	# remove all genes not expressed in all samples
	for reg in rs_data[feat].keys():
		if len( rs_data[feat][reg].split( ',' ) ) != no_samples:
			print >> sys.stderr, "Deleting region %s of gene %s due to poor representation..." % ( reg, feat )
			del rs_data[feat][reg]
	
	# remove duplicates by selecting the longest gene
	if len( rs_data[feat] ) > 1:
		max_reg = find_max( rs_data[feat].keys() )
		for reg in rs_data[feat].keys():
			if reg != max_reg:						# if the region is not the max
				print >> sys.stderr, "Deleting region %s of gene %s due to duplication..." % ( reg, feat )
				del rs_data[feat][reg]			# delete it
				
	# if the feat is empty delete it as well
	if len( rs_data[feat].keys() ) == 0:
		print >> sys.stderr, "Deleting gene %s..." % feat
		del rs_data[feat]
"""

# combine the data
c = 0
for g in g2ps:
	if c > 100: break
	
	# get the RNA-Seq values
	try:
		#rs = rs_data[g].values()[0]
		rs = rs_data[g]
	except KeyError:
		print >> sys.stderr, g
		continue
	
	# get probe intensities
	probes = ""
	intensities = ""
	for ps in g2ps[g]:
		try:
			probe_ids = ps2int[ps].keys()
		except KeyError:
			continue
		for p in probe_ids:
			if probes == '':
				probes = p
				intensities = ps2int[ps][p]
			else:
				probes += "," + p
				intensities += "," + ps2int[ps][p]

	if probes == "":
		print "\t".join( [ g, str( no_samples ), "0", "NA", rs, "NA" ] )
	else:
		print "\t".join( [ g, str( no_samples ), str( len( probes.split( ',' )) ), probes, rs, intensities ] )
	c += 0
