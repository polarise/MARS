#!/usr/bin/env python

def f( in_file, out_file ):
	the_file = open( in_file )
	other_file = open( out_file, 'w' )
	
	for row in the_file:
		print >> other_file, row.strip()
	
	the_file.close()
	other_file.close()
	
	return

f( "infile.txt", "outfile3.txt" )
f( "infile.txt", "outfile4.txt" )
