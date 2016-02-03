#!/usr/bin/env python
from __future__ import division
import sys
import re
import gzip

## BUGGY!!!!

class Probe( object ):
	def __init__( self, probe_id, _type, gc_count, probe_length, interrogation_position, probe_sequence ):
		self.probe_id = probe_id
		self.type = _type
		self.gc_count = gc_count
		self.probe_length = probe_length
		self.interrogation_position = interrogation_position
		self.probe_sequence = probe_sequence
	
	def __repr__( self ):
		probe_str = "Probe ID: %s\nProbe type: %s\nGC count: %s\nProbe length: %s\nInter. Pos.: %s\nSequence: %s" % ( self.probe_id, self.type, self.gc_count, self.probe_length, self.interrogation_position, self.probe_sequence )
		return probe_str

class Probeset( object ):
	def __init__( self, probeset_id, _type ):
		self.probeset_id = probeset_id
		self.type = _type
		self.probeset_name = None
		
		self.probes = dict()
	
	def add_probe( self, atom_id, Probe ):
		if atom_id not in self.probes:
			self.probes[ atom_id ] = [ Probe ]
		else:
			self.probes[ atom_id ] += [ Probe ]
		
	def __repr__( self ):
		return "Probeset: %s; No. of probes: %s" % ( self.probeset_id, len( self.probes ))

class PGF( object ):
	def __init__( self, pgf_fn ):
		self.name = pgf_fn
		self.chip_type = list()
		self.lib_set_name = None
		self.lib_set_version = None
		self.create_date = None
		self.guid = None
		self.pgf_format_version = None
		self.rows = None
		self.cols = None
		self.probesets = None
		self.datalines = None
		self.sequential = None
		self.order = None
		
		self.Probesets = dict()
		
	def read_header( self ):
		def get_item( row ):
			return row.strip().split( "=" )[1]
		
		with open( self.name ) as f:
			for row in f:
				if re.search( "chip_type", row.strip() ): self.chip_type += [ get_item( row ) ]
				elif re.search( "lib_set_name", row.strip() ): self.lib_set_name = get_item( row )
				elif re.search( "lib_set_version", row.strip() ): self.lib_set_version = get_item( row )
				elif re.search( "create_date", row.strip() ): self.create_date = get_item( row )
				elif re.search( "guid", row.strip() ): self.guid = get_item( row )
				elif re.search( "pgf_format_version", row.strip() ): self.pgf_format_version = get_item( row )
				elif re.search( "rows", row.strip() ): self.rows = get_item( row )
				elif re.search( "cols", row.strip() ): self.cols = get_item( row )
				elif re.search( "probesets", row.strip() ): self.probesets = get_item( row )
				elif re.search( "datalines", row.strip() ): self.datalines = get_item( row )
				elif re.search( "sequential", row.strip() ): self.sequential = get_item( row )
				elif re.search( "order", row.strip() ): self.order = get_item( row )
				elif row[0] != "#": break
		
		return			
	
	def get_data( self ):
		numbers = map( str, range( 1, 10 ))
		
		with open( self.name ) as f:
			c = 0
			for row in f:
				if c > 100: break
				if row[0] == "#": continue
				
				if row[0] in numbers:
					L = row.strip().split( "\t" )
					probeset = Probeset( L[0], L[1] )
					
					row = f.next()

					while row[0] not in numbers:
						L = row.split( "\t" )
						if L[1] != "":
							atom_id = L[1]
						elif L[1] == "":
							probe = Probe( L[2], L[3], L[4], L[5], L[6], L[7].strip() )
							probeset.add_probe( atom_id, probe )
						try:
							row = f.next()
						except StopIteration:
							break
					
					self.Probesets[ probeset.probeset_id ] = probeset
				c += 0
				
		return self
	
	def overlap( self, other_pgf ):
		c = 0
		for P in self.Probesets:
			if c > 1000: break
			print self.Probesets[ P ]
			print "Probes in %s..." % self.name
			for p in self.Probesets[ P ].probes:
				for r in self.Probesets[ P ].probes[ p ]:
					print r
			print
			print "Probes in %s..." % other_pgf.name
			for q in other_pgf.Probesets[ P ].probes:
				for s in other_pgf.Probesets[ P ].probes[ q ]:
					print s
			print
			c += 1
		
	def tabbed_output( self, out_fn="" ):
		if out_fn == "":
			f = sys.stdout
		else:
			f = open( out_fn, 'w' )
		
		for pss in self.Probesets:
			for ps in self.Probesets[ pss ].probes:
				for p in self.Probesets[ pss ].probes[ ps ]:
					print >> f, "\t".join([ self.Probesets[ pss ].probeset_id, p.probe_id, p.type, p.gc_count, p.probe_length, p.interrogation_position, p.probe_sequence ])
			
		f.close()
	
	def __repr__( self ):
		pgf_str = "Filename: %s\nChip type: %s\nLibrary set name: %s\nLibrary version: %s\nCreated: %s\nGUID: %s\nPGF version: %s\nNo. of rows: %s\nNo. of cols: %s\nNo. of probesets: %s\nNo. of datalines: %s\nSequential: %s\nOrder: %s" % ( self.name, ",".join( self.chip_type ), self.lib_set_name, self.lib_set_version, self.create_date, self.guid, self.pgf_format_version, self.rows, self.cols, self.probesets, self.datalines, self.sequential, self.order )
		return pgf_str
	
if __name__ == "__main__":
	try:
		pgf_fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "Usage: ./script.py <pgf_fn1>"
		sys.exit( 0 )
	
	pgf = PGF( pgf_fn )
	
	pgf.read_header()

	pgf.get_data()
	
	pgf.tabbed_output( pgf_fn + ".txt" )
	
	
	
	


