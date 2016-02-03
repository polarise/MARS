#!/home/paulk/software/bin/python
import sys
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument('-c','--command',help="the command to be run (include multiple lines)")
parser.add_argument('-t','--template',default="/home/paulk/RP/bash_general/empty_subjob.sh",help="the template to be used")
parser.add_argument('-g','--generate-template',action='store_true',help="print to stdout a blank template for use")
parser.add_argument('-j','--job-name',default="my_job_%s" % date,help="the name of the job")
parser.add_argument('-n','--script-name',default="my_script_%s.sh" % date,help="the name of the script files") # check this
parser.add_argument('-s','--submit',default=False,action='store_true',help="submit the job if you're sure it's OK")
parser.add_argument('-p','--path',default=".",help="default path to use; assumes that all data and destination files will be at the same path")
parser.add_argument('-b','--bin-path',default=".",help="path to executible; assumes only one executible")

args = parser.parse()
command = args.command
template = args.template
generate_template = args.generate_template
job_name = args.job_name

if 




# the template
f = open("/home/paulk/RP/bash_general/empty_subjob.sh")
job = "".join(f.readlines())
f.close()

print job
#cmd = "/home/paulk/tophat-1.4.1.Linux_x86_64/tophat -o /data1/paulk/%s -p 10 /home/paulk/bowtie-0.12.7/scripts/hg19 %s"
path = "/data2/paul/raw_hapmap_sequences/pickrell/fastq/%s"
cmd = "/home/paulk/software/bin/cufflinks -o /data1/paulk/%s -p 10 -G /data2/paulk/RP/resources/Homo_sapiens.GRCh37.66.full_chr.gtf -b /data2/paulk/RP/resources/refs/hg19/Homo_sapiens.GRCh37.66.dna.chr.fa -u /data1/paulk/%s/accepted_hits.bam"

f = open(sys.argv[1])
for row in f:
	l = row.strip().split('\t')
	#g = open("cufflinks_%s.sh" % l[0],'w')
	print >> g,job % (cmd % (l[0],l[0]))
	#print >> g,job % (cmd % (l[0],",".join(map(lambda(x):path % x,l[-1].split(',')))))
	#g.close()
f.close()
