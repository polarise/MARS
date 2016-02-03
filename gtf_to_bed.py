#!/usr/bin/env python
from __future__ import division
import sys
import subprocess
from key_functions import process_feature

"""
Take a GTF file and:
- convert it to a BED file
- sort it
- index it
"""

try:
	fn_s = sys.argv[1]
	try:
		assert fn_s[-3:].upper() == 'GTF'
	except:
		raise ValueError( "The GTF file must have a 'gtf' extension." )
		sys.exit( 2 )
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <GTF_filename>"
	sys.exit( 1 )
	
unwanted = [ '#', '-' ]

# create a BED file
f = open( fn_s )
new_fn_s = fn_s[:-4]
g = open( new_fn_s + ".bed", 'w' )
c = 0
for row_s in f:
	if c > 1000: break
	if row_s[0] in unwanted: continue
	L = row_s.strip().split( '\t' )
	features = process_feature( row_s )
	print >> g, "\t".join( [ L[0], L[3], L[4], features['gene_id'], features['transcript_id'], L[6] ])
	c += 0
f.close()
g.close()

cmd = "sort -k1,1 -k2,2n %s | bgzip > %s.gz" % (( new_fn_s + ".bed", )*2)
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
p.communicate()

cmd = "tabix -p bed %s.gz" % ( new_fn_s + ".bed" )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
p.communicate()

print >> sys.stderr, "Done!"



