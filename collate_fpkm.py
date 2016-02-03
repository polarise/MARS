#!/usr/bin/env python
from __future__ import division
import sys
import cPickle
import argparse
import time

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

no_samples = 8
sample_names = map( lambda x: str( x ) + "C", range( 3, 11 ))
rs_data = dict()
for sn in sample_names:
	c = 0
	f = open( "/data2/paulk/RP/output2/%s/Genic/genes.fpkm_tracking" % sn )
	for row in f:
		if c > 100: break
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
		c += 0
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
			if reg != max_reg:					      # if the region is not the max
				print >> sys.stderr, "Deleting region %s of gene %s due to duplication..." % ( reg, feat )
				del rs_data[feat][reg]		  # delete it
				
	# if the feat is empty delete it as well
	if len( rs_data[feat].keys() ) == 0:
		print >> sys.stderr, "Deleting gene %s..." % feat
		del rs_data[feat]

for feat in rs_data:
	for coord in rs_data[feat]:
		print "\t".join( [feat] + rs_data[feat][coord].split( "," ) )
