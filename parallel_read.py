#!/usr/bin/env python
import multiprocessing

def f( in_file, out_file ):
	the_file = open( in_file )
	other_file = open( out_file, 'w' )
	
	for row in the_file:
		print >> other_file, row.strip()
	
	the_file.close()
	other_file.close()
	
	return

p1 = multiprocessing.Process( target=f, args=( "infile.txt", "outfile1.txt", ) )
p1.start()

p2 = multiprocessing.Process( target=f, args=( "infile.txt", "outfile2.txt", ) )
p2.start()

p1.join()
p2.join()
