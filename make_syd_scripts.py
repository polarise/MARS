names = [i.strip() for i in open("../0_prep_data/CEU_names.txt")]
for ID in names:
	f = open("cufflinks_%s.sh" % ID,'w')
	print >> f,"""\
#!/bin/bash

#$ -N cufflinks
#$ -q all.q
#$ -cwd
#$ -v PATH
#$ -v LD_LIBRARY_PATH
"""

	print >> f,"/home/paulk/software/bin/cufflinks -o /data1/paulk/Montgomery_Data/output/%s/tophat_out -p 10 -G /data2/paulk/RP/resources/Homo_sapiens.GRCh37.66.full_chr.gtf -b /data2/paulk/RP/resources/refs/hg19/Homo_sapiens.GRCh37.66.dna.chr.fa -u /data1/paulk/Montgomery_Data/output/%s/tophat_out/accepted_hits.bam" % (ID,ID)
	f.close()

