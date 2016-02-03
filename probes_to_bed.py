#!/usr/bin/env python
from __future__ import division
import sys
import subprocess

"""
Take a probe intensity file and:
- convert it to a BED file
- sort it
- index it
"""

try:
	fn_s = sys.argv[1]
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <probe_intensities_filename>"
	sys.exit( 1 )
	
unwanted = [ '#', ]

# create a BED file
f = open( fn_s )
g = open( fn_s + ".bed", 'w' )
c = 0
for row_s in f:
	if c > 1000: break
	if row_s[0] in unwanted: continue
	L = row_s.strip().split( '\t' )
	probe_values_L = L[7:]
	if row_s[0] == 'p':
		no_samples_i = len( probe_values_L )
		continue
	print >> g, "\t".join( [ "chrP", L[4], L[4], L[0], ",".join( probe_values_L ) ])
	c += 0
f.close()
g.close()

cmd = "sort -k2,2n %s | bgzip > %s.gz" % (( fn_s + ".bed", )*2)
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
p.communicate()

cmd = "tabix -p bed %s.gz" % ( fn_s + ".bed" )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )

print >> sys.stderr, "Done!"
