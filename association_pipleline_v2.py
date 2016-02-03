#!/usr/bin/env python
from __future__ import division
import sys
import os
import time
import subprocess

# global variables
HOME_s = os.environ['HOME']
PWD_s = os.environ['PWD']

genotypes_dir_s = HOME_s + "/RP/resources/variation"
gene_coordsfn_s = HOME_s + "/RP/resources/Homo_sapiens.GRCh37.66.full_chr.bed"

# get the input files
try:
	chrom_s = sys.argv[1]
	orderfn_s = sys.argv[2]
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <chrom> <expressionfn_s>"
	sys.exit(1)

try:
	maf_f = float( sys.argv[3] )
except IndexError:
	maf_f = 0.05

f = open( orderfn_s )
if f.readline().strip().split( '\t' )[4] == 'Yoruban':
	pop_s = 'YRI'
else:
	pop_s = 'CEU'
f.close()

genotypesfn_s = genotypes_dir_s + "/%s/genotypes_%s_%s_r28_nr.b36_fwd.txt.gz" % ( pop_s, chrom_s, pop_s )

print >> sys.stderr, "[%s]\tRun variables:" % time.ctime( time.time() )
MY_PYTHON_PATH_s = HOME_s + "/MARS/python_scripts"
MY_R_PATH_s = HOME_s + "/MARS/R_scripts"
print >> sys.stderr, "Genotypes file:   %s" %( genotypesfn_s )
print >> sys.stderr, "Order file:       %s" % ( orderfn_s )
print >> sys.stderr, "Gene coordinates: %s" % ( gene_coordsfn_s )
print >> sys.stderr
print >> sys.stderr, "Current directory is %s." % ( os.environ['PWD'] )
print >> sys.stderr, "Python scripts path is %s." % ( MY_PYTHON_PATH_s )
print >> sys.stderr, "R scripts path is %s." % ( MY_R_PATH_s )

# creating a folder for this chromosome
print >> sys.stderr, "[%s]\tCreating a new directory for '%s' data... " % ( time.ctime(time.time()), chrom_s ),
try:
	os.mkdir( chrom_s )
except OSError:
	print >> sys.stderr
	print >> sys.stderr, "[%s]\tDirectory already exists..." % ( time.ctime( time.time() )),
print >> sys.stderr, "done."

print >> sys.stderr, "[%s]\tChanging to newly created directory..." % time.ctime(time.time()),
os.chdir( chrom_s )
print >> sys.stderr, "done."

# convert genotypes file to BED
print >> sys.stderr, "[%s]\tConverting genotypes file to BED... " % time.ctime(time.time()),
cmd = "python %s/snp_coord_to_bed.py %s > genotypes_%s.bed" % ( MY_PYTHON_PATH_s, genotypesfn_s, pop_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# convert BED to indexed tabix
print >> sys.stderr, "[%s]\tConverting BED to indexed tabix... " % time.ctime(time.time()),
cmd = "sort -k1,1 -k2,2n genotypes_%s.bed | bgzip > genotypes_%s.bed.gz" % ( pop_s, pop_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()

cmd = "tabix -p bed genotypes_%s.bed.gz" % ( pop_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# convert genotypes into numbers
print >> sys.stderr, "[%s]\tConverting genotypes into integers... " % time.ctime(time.time()),
cmd = "python %s/genotypes_to_no.py %s %s %s > converted_genotypes_%s.txt 2> converted_genotypes_%s.miss" % ( MY_PYTHON_PATH_s, genotypesfn_s, orderfn_s, maf_f, pop_s, pop_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

print >> sys.stderr, "[%s]\tReturning to previous directory... " % time.ctime(time.time()),
os.chdir( ".." )
print >> sys.stderr, "done."
print >> sys.stderr

