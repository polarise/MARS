#!/usr/bin/env python
from __future__ import division
import sys

if __name__ == "__main__":
	unwatned = [ '#', 'p' ]	
	
	try:
		huge11_fn = sys.argv[1]
		huge10_fn = sys.argv[2]
		huge10_pi_fn = sys.argv[3]	# pi = probe intensities
		new_huge10_pi_fn = sys.argv[4]
	except IndexError:
		print >> sys.stderr, """\
Usage: ./script.py <huge11_fn> <huge10_fn> <huge10_pi_fn> <new_huge10_pi_fn>"""
		sys.exit( 1 )
	
	d11 = {}
	f = open( huge11_fn )
	c = 0
	for row in f:
		if c > 100: break
		if row[0] in unwatned: continue
		L = row.strip().split( "\t" )
		ps_id = L[0]
		p_id = L[4]
		seq = L[9]
		if L[0] not in d11:
			d11[ ps_id ] = {}
		d11[ ps_id ][ seq ] = p_id
		c += 0
	f.close()
			
	d10 = dict()
	f = open( huge10_fn )
	c = 0
	for row in f:
		if c > 100: break
		if row[0] in unwatned: continue
		L = row.strip().split( "\t" )
		ps_id = L[0]
		p_id = L[3]
		seq = L[8]
		if L[0] not in d10:
			d10[ ps_id ] = {}
		d10[ ps_id ][ p_id ] = seq
		c += 0
	f.close()
	
	f = open( huge10_pi_fn )
	g = open( new_huge10_pi_fn, 'w' )
	for row in f:
		if row[0] in unwatned: 
			print >> g, row.strip()
			continue
		L = row.strip().split( "\t" )
		ps_id = L[4]
		p_id = L[0]
		seq = d10[ ps_id ][ p_id ]
		new_p_id = d11[ ps_id ][ seq ]
		print >> g, "\t".join( [ new_p_id ] + L[1:] )
	f.close()
	g.close()
