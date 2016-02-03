#!/usr/bin/env python
from __future__ import division
import sys
import os
import time
import subprocess

# global variables
HOME_s = os.environ['HOME']
PWD_s = os.environ['PWD']

genotypes_dir_s = HOME_s + "/RP/resources/variation/CEU"
#gene_coordsfn_s = HOME_s + "/RP/resources/Entrez_genes.txt" # **
gene_coordsfn_s = HOME_s + "/RP/resources/Homo_sapiens.GRCh37.66.full_chr.bed"
#orderfn_s = HOME_s + "/MARS/0_prep_data/combined_reduced_CEU.pauls.txt" # **
orderfn_s = HOME_s + "/MARS/0_prep_data/combined_reduced_CEU.txt"

# get the input files
try:
	chrom_s = sys.argv[1]
	expressionfn_s = sys.argv[2]
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <chrom> <expressionfn_s>"
	sys.exit(1)

genotypesfn_s = genotypes_dir_s + "/genotypes_%s_CEU_r28_nr.b36_fwd.txt.gz" % ( chrom_s )

try:
	aug_i = int( sys.argv[3] )
except IndexError:
	aug_i = 10000
	
try:
	maf_f = float( sys.argv[4] )
except IndexError:
	maf_f = 0.05

print >> sys.stderr, "[%s]\tRun variables:" % time.ctime( time.time() )
MY_PYTHON_PATH_s = HOME_s + "/MARS/python_scripts"
MY_R_PATH_s = HOME_s + "/MARS/R_scripts"
print >> sys.stderr, "Genotypes file: %s" %( genotypesfn_s )
print >> sys.stderr, "Order file: %s" % ( orderfn_s )
print >> sys.stderr, "Gene coordinates: %s" % ( gene_coordsfn_s )
print >> sys.stderr, "Gene expression: %s" % ( expressionfn_s )
print >> sys.stderr, "Searching for SNPs within %s bases" % ( aug_i )
print >> sys.stderr
print >> sys.stderr, "Current directory is %s." % ( os.environ['PWD'] )
print >> sys.stderr, "Python scripts path is %s." % ( MY_PYTHON_PATH_s )
print >> sys.stderr, "R scripts path is %s." % ( MY_R_PATH_s )

# creating a folder for this chromosome
print >> sys.stderr, "[%s]\tCreating a new directory for '%s' data... " % ( time.ctime(time.time()), chrom_s ),
try:
	os.mkdir( chrom_s )
except OSError:
	print >> sys.stderr, "directory already exists!"
	sys.exit(1)
print >> sys.stderr, "done."

print >> sys.stderr, "[%s]\tChanging to newly created directory..." % time.ctime(time.time()),
os.chdir( chrom_s )
print >> sys.stderr, "done."

# convert genotypes file to BED
print >> sys.stderr, "[%s]\tConverting genotypes file to BED... " % time.ctime(time.time()),
cmd = "python %s/snp_coord_to_bed.py %s > genotypes.bed" % ( MY_PYTHON_PATH_s, genotypesfn_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# convert BED to indexed tabix
print >> sys.stderr, "[%s]\tConverting BED to indexed tabix... " % time.ctime(time.time()),
cmd = "sort -k1,1 -k2,2n genotypes.bed | bgzip > genotypes.bed.gz"
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()

cmd = "tabix -p bed genotypes.bed.gz"
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# convert genotypes into numbers
print >> sys.stderr, "[%s]\tConverting genotypes into integers... " % time.ctime(time.time()),
cmd = "python %s/genotypes_to_no.py %s %s %s > converted_genotypes.txt 2> converted_genotypes.miss" % ( MY_PYTHON_PATH_s, genotypesfn_s, orderfn_s, maf_f )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# created augmented convert gene coordinates to BED
print >> sys.stderr, "[%s]\tConverting gene coordinates to BED... " % time.ctime(time.time()),
cmd = "python %s/augment_gene_coords_to_bed.py %s %s > aug_gene_coords.bed" % ( MY_PYTHON_PATH_s, gene_coordsfn_s, aug_i )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# convert BED to indexed tabix
print >> sys.stderr, "[%s]\tConverting BED indexed tabix... " % time.ctime(time.time()),
cmd = "sort -k1,1 -k2,2n aug_gene_coords.bed | bgzip > aug_gene_coords.bed.gz"
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()

cmd = "tabix -p bed aug_gene_coords.bed.gz"
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# get overlapping SNP
print >> sys.stderr, "[%s]\tObtaining overlapping SNPs... " % time.ctime(time.time()),
cmd = "python %s/overlapping_snps.py genotypes.bed.gz aug_gene_coords.bed.gz > overlapping.txt 2> overlapping.miss" % ( MY_PYTHON_PATH_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# combine all data
print >> sys.stderr, "[%s]\tCombining all data together... " % time.ctime(time.time()),
cmd = "python %s/combine_all_data.py overlapping.txt converted_genotypes.txt %s > combined_data.txt 2> combined_data.miss" % ( MY_PYTHON_PATH_s, expressionfn_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# test for associations
print >> sys.stderr, "[%s]\tPerforming simple association tests... " % time.ctime(time.time()),
cmd = "/cm/shared/apps/R/2.15.0/bin/Rscript %s/association_test.r combined_data.txt" % ( MY_R_PATH_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# converting results from Rdata to text
print >> sys.stderr, "[%s]\tConverting results from Rdata to text... " % time.ctime(time.time()),
cmd = "/cm/shared/apps/R/2.15.0/bin/Rscript %s/Rdata_to_text.r association_results.Rdata association_results.txt" % ( MY_R_PATH_s )
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

print >> sys.stderr, "[%s]\tReturning to previous directory... " % time.ctime(time.time()),
os.chdir( ".." )
print >> sys.stderr, "done."
print >> sys.stderr

