#!/usr/bin/env python
from __future__ import division
import sys
import os
import time
import subprocess

try:
	chrom_s = sys.argv[1]
	expressionfn_s = sys.argv[2]
except IndexError:
	print >> sys.stderr, "Usage: ./script.py <chrom> <expressionfn_s>"
	sys.exit(1)

try:
	aug_i = int( sys.argv[3] )
except IndexError:
	aug_i = 10000
	
# global variables
HOME_s = os.environ['HOME']
PWD_s = os.environ['PWD']

genotypes_dir_s = HOME_s + "/RP/resources/variation"
gene_coordsfn_s = HOME_s + "/RP/resources/Homo_sapiens.GRCh37.66.full_chr.bed"

MY_PYTHON_PATH_s = HOME_s + "/MARS/python_scripts"
MY_R_PATH_s = HOME_s + "/MARS/R_scripts"
print >> sys.stderr, "Gene expression:  %s" % ( expressionfn_s )
print >> sys.stderr, "Searching for SNPs within %s bases" % ( aug_i )
print >> sys.stderr

print >> sys.stderr, "[%s]\tChanging to created directory '%s'..." % ( time.ctime(time.time()), chrom_s ),
os.chdir( chrom_s )
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

# combine genotypes
print >> sys.stderr, "[%s]\tCombining genotypes... " % time.ctime(time.time()),
cmd = "bedtools intersect -a genotypes_YRI.bed -b genotypes_CEU.bed | sort -k1,1 -k2,2n | bgzip > genotypes.bed.gz"
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

cmd = "tabix -p bed genotypes.bed.gz"
p = subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True )
output = p.communicate()
print >> sys.stderr, "done."

# combining converted genotypes file
print >> sys.stderr, "[%s]\tCombining converted genotypes... " % time.ctime(time.time()),
cmd = "python %s/combine_genotypes.py converted_genotypes_YRI.txt converted_genotypes_CEU.txt > converted_genotypes.txt 2> converted_genotypes.miss" % ( MY_PYTHON_PATH_s )
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
cmd = "/cm/shared/apps/R/2.15.0/bin/Rscript %s/association_test_MIX.r combined_data.txt" % ( MY_R_PATH_s )
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
