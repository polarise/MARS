#!/home/paulk/software/bin/python
from __future__ import division
from sys import argv,exit,stderr

try:
    pgf_file = argv[1]
    out_file = argv[2]
except IndexError:
    print >> stderr,"""\
Script to collate probes by probeset from a PGF file.
Usage: ./<script.py> <PGF_file> <output_file>"""
    exit(1)

f = open(pgf_file)
c = 0
probes = dict()
current_probeset = ''
count_probesets = 0
count_probes = 0
for line in f:
    if line[0] == '#': continue
    if c > 20: break
    l = line.strip().split('\t')
    #print len(l),":",l
    length = len(l)
    if length == 2:
        current_probe = l[0]
        probes[current_probe] = [] # start simple; just get the probeset to probe mappings; can be modified to get the GC, the sequence etc.
        count_probesets += 1
    elif length == 1:
        pass
    elif length == 6:
        probes[current_probe] += [l[0]]
        count_probes += 1
    c += 0
f.close()

f = open(argv[2],'w')
for p in probes:
    print >> f,p+"\t"+",".join(probes[p])
f.close()

print >> stderr,"Written %s probesets with a total of %s probes to %s"%(count_probesets,count_probes,out_file)
