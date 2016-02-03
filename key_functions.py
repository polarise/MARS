from __future__ import division
from sys import stderr
from scipy import *
from scipy import stats
import pysam
#import rpy2.robjects as R
#options = R.r['options']
#options(warn=-1)

"""
TODO: Document this file
"""

def motif_abundance(seq,motif_length):
	"""
	given a sequence and a sequence length return the N-length motifs and the number of times
	they occur and the total number of possible motifs
	"""
	motifs = dict()
	total_possible_motifs = 0
	for i in xrange(len(seq)-motif_length+1):
		motif = seq[i:i+motif_length]
		if motif not in motifs:
			motifs[motif] = 1
		else:
			motifs[motif] += 1
		total_possible_motifs += 1
	return motifs,total_possible_motifs

def n_ary_sum(N1,N2,digits,base):
	"""
	given two numbers in base returns a string of their sum in that base
	when overflow_digits == False will return a string of length as max(N1,N2)
	"""
	N1_list = map(int,list(str(N1).zfill(digits)))
	N2_list = map(int,list(str(N2).zfill(digits)))
	
	# sanity checks
	try:
		assert len(N1_list) == len(N2_list)
	except:
		raise ValueError("Input digits of unequal length")
		sys.exit(1)
	
	if base in N1_list or base in N2_list:
		raise ValueError("Base digit cannot appear in numbers")
		sys.exit(1)

	carry = [0] * (digits + 1)
	result = [0] * (digits + 1)

	indices = range(digits)[::-1]

	for i in indices:
		result[i + 1] = (N1_list[i] + N2_list[i] + carry[i + 1]) % base
		carry[i] = (N1_list[i] + N2_list[i] + carry[i + 1]) // base
	result[0] = carry[0]

	return "".join(map(str,result)).zfill(digits + 1)

def motif_generator(alpha,size):
	base = len(alpha)
	
	instances = list()
	result = 0
	for i in xrange(base**size):
		result = n_ary_sum(result,1,size,base)[1:]
		instance = ""
		for r in result:
			instance += alpha[int(r)]
		instances.append(instance)
	return instances	

def get_probe_values(results):
	values = dict()
	for row in results:
		l = row.split('\t')
		values[int(l[3])] = map(float,l[4:])
	return values

def intron_length(region1,region2,pos1,pos2):
	"""
	compute the distance given two region strings and where to start 
	pos1 is the coordinate on the first region
	pos2 is ditto
	length is returned as a string
	"""
	c1,st1_sp1,sd1 = region1.split(':')
	c2,st2_sp2,sd2 = region2.split(':')
	st1,sp1 = map(int,st1_sp1.split('-'))
	st2,sp2 = map(int,st2_sp2.split('-'))
	if c1 != c2: raise ValueError("Conflict in chromosome names")
	if sd1 != sd2: raise ValueError("Conflict in strands")
	if sd1 == '+': return str(abs((st1+pos1) - (sp2-pos2))+1)
	elif sd1 == '-': return str(abs((st2+pos2) - (sp1-pos1))+1)
	
def simple_intron_length(region1,region2):
	c1,st1_sp1,sd1 = region1.split(':')
	c2,st2_sp2,sd2 = region2.split(':')
	st1,sp1 = map(int,st1_sp1.split('-'))
	st2,sp2 = map(int,st2_sp2.split('-'))
	if c1 != c2: raise ValueError("Conflict in chromosome names")
	if sd1 != sd2: raise ValueError("Conflict in strands")
	if sd1 == '+': return str(abs(sp1-st2)-1)
	if sd2 == '-': return str(abs(sp2-st1)-1)

def PrintStatic(line):
	stderr.write("\r%s"%line.ljust(50))
	stderr.flush()

def cigar2end(cigar,start):
	"""
	given simple CIGAR string (with the M and N only) returns the position of the last nucleotide
	modification: if it denotes a junction will return the exon junctions
	"""
	pass1 = cigar.split('N')
	pass2 = [i.split('M') for i in pass1]
	pass3 = []
	for i in pass2:
		for j in i:
			if j != '':
				pass3 += [int(j)]
	length = sum(pass3)
	# for the case of one jxn there are two additional co-ordinates
	exon1_sp = pass3[0]
	exon2_st = pass3[2]
	end = start+length-1
	return end,start+exon1_sp-1,end-exon2_st

def process_bed_rows( rows ):
	"""
	takes a set of bed rows; returns hidden values
	"""
	col3 = list()
	col4 = list()
	for row in rows:
		L = row.split( '\t' )
		col3 += [ L[3] ]
		col4 += [ L[4] ]
	return col3, col4


def process_feature(feature):
	"""
	given a row from a GTF file returns a dictionary of the items available together with a key called 'keys' of all the keys found
	can be made to be more elegant
	"""
	feature_dict = dict()
	line = feature.strip().split('\t')
	Annot = [k.strip().strip(';').strip('"') for k in line[8].strip().split(' ')]
	if 'gene_id' in Annot: feature_dict['gene_id'] = Annot[Annot.index('gene_id')+1]
	if 'transcript_id' in Annot: feature_dict['transcript_id'] = Annot[Annot.index('transcript_id')+1]
	if 'exon_number' in Annot: feature_dict['exon_number'] = Annot[Annot.index('exon_number')+1]
	if 'gene_name' in Annot: feature_dict['gene_name'] = Annot[Annot.index('gene_name')+1]
	if 'tx_name' in Annot: feature_dict['tx_name'] = Annot[Annot.index('tx_name')+1]
	if 'locus_id' in Annot: feature_dict['locus_id'] = Annot[Annot.index('locus_id')+1]
	if 'reads' in Annot: feature_dict['reads'] = Annot[Annot.index('reads')+1]
	if 'length' in Annot: feature_dict['length'] = Annot[Annot.index('length')+1]
	if 'RPKM' in Annot: feature_dict['RPKM'] = Annot[Annot.index('RPKM')+1]
	feature_dict['keys'] = feature_dict.keys()
	return feature_dict

def region_to_GTF(name,region,transcript=False):
	"""
	name is 'ENSGXXXXXXX:<exon_number>'
	regions are defined by chr6:56399919-56400089:-
	"""
	feat,exno = name.split(':')
	chrom,st_sp,sd = region.split(':')
	st,sp = map(lambda n:n+1,map(int,st_sp.split('-')))
	if not transcript:
		gtf_line = "%s\tprotein_coding\texon\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; transcript_id \"-\"; exon_number \"%s\"; gene_name \"-\"; gene_biotype \"protein_coding\"; transcript_name \"-\";" % (chrom,st,sp,sd,feat,exno)
	elif transcript:
		gtf_line = "%s\tprotein_coding\texon\t%s\t%s\t.\t%s\t.\tgene_id \"-\"; transcript_id \"%s\"; exon_number \"%s\"; gene_name \"-\"; gene_biotype \"protein_coding\"; transcript_name \"-\";" % (chrom,st,sp,sd,feat,exno)
	return gtf_line	

def terminal_exon(exon,gex_coords):
	"""
	check whether the given exon is either the first or the last
	if the first or the last return (True,None,None)
	else return (False,<prev_exon>,<next_exon>)
	"""
	gene,the_exon = exon.split(':')
	exno = int(the_exon)
	if exno == 1: return True,None,None
	try:
		next_exon = gex_coords[gene+':'+str(exno+1)] # check if it's the last exon
	except KeyError:
		next_exon = None				# it is the last exon
	if next_exon == None:
		return True,None,None
	else:
		return False,gene+':'+str(exno-1),gene+':'+str(exno+1)

def parse_lines(lines,strand,get_transcripts=False):
	"""
	given a batch of lines separated by newlines return either gene-exons or true-exons together with the number of lines processed
	"""
	no_lines = 0
	gene_exons = set()
	for l in lines:
		line = l.strip().split('\t')
		if strand == '*':
			pass
		elif line[6] != strand:
			continue
		try:
			Annot = [k.strip().strip(';').strip('"') for k in line[8].strip().split(' ')]
		except IndexError:
			raise ValueError("offending line %s" % l.strip())
		if 'gene_id' in Annot: gene_id = Annot[Annot.index('gene_id')+1]
		if 'transcript_id' in Annot: tx_id = Annot[Annot.index('transcript_id')+1]
		if 'exon_number' in Annot: exon_no = Annot[Annot.index('exon_number')+1]
		if 'gene_name' in Annot: gene_name = Annot[Annot.index('gene_name')+1]
		if 'tx_name' in Annot: tx_name = Annot[Annot.index('tx_name')+1]
		if not get_transcripts: gene_exons.add(gene_id+":"+exon_no)
		else: gene_exons.add(tx_id+":"+exon_no)
		no_lines += 1
	return gene_exons,no_lines

def complementDNA(dna):
	compDNA = ''
	for n in dna:
		if n == 'A':
			compDNA += 'T'
		elif n == 'C':
			compDNA += 'G'
		elif n == 'G':
			compDNA += 'C'
		elif n == 'T':
			compDNA += 'A'
	return compDNA

def augment_region(region,extend):
	chrom,st_sp,sd = region.split(':')
	st,sp = st_sp.split('-')
	nst = str(int(st)-extend)
	nsp = str(int(sp)+extend)
	return chrom+":"+nst+"-"+nsp+":"+sd

def GTFrow_to_5p3pcoords(GTFrow,offset=(25,25,25,25),full_chrom_names=True):
	"""
	given a GTFrow and an offsets about the splice junctions returns 5' and 3' splice site boundaries coordidates
	0-based
	"""
	l = GTFrow.strip().split('\t')
	chrom = l[0]
	if full_chrom_names:
		if chrom[:3] != "chr":
			chrom = "chr"+chrom
	else:
		if chrom[:3] == "chr":
			chrom = chrom[3:]
	st = int(l[3])
	sp = int(l[4])
	sd = l[6]
	if sd == '+':
		five_st = str(sp-offset[0]+1)
		five_sp = str(sp+offset[1])
		three_st = str(st-offset[2])
		three_sp = str(st+offset[3]-1)
	elif sd == '-':
		five_st = str(st-offset[1])
		five_sp = str(st+offset[0]-1)
		three_st = str(sp-offset[3]+1)
		three_sp = str(sp+offset[2])
	#annot = l[8]
	five = chrom+":"+five_st+"-"+five_sp+":"+sd
	three = chrom+":"+three_st+"-"+three_sp+":"+sd
	return five,three

def normalise(sm_feat_data,lg_feat_data,sm_lg_map,kind):
	normalised_data = dict()
	missing_lg = list()
	for sm in sm_feat_data:
		try:
			lg = sm_lg_map[sm]
			nd = list()
			for i in xrange(len(sm_feat_data[sm])):
				try:
					nd += [(sm_feat_data[sm][i]-lg_feat_data[lg][i]) if kind == 'difference' else (sm_feat_data[sm][i]/lg_feat_data[lg][i])]
				except ZeroDivisionError:
					nd += ['NA']
			normalised_data[lg,sm] = nd
			#normalised_data[lg,sm] = array(sm_feat_data[sm])/array(lg_feat_data[lg])
		except KeyError:
			missing_lg.append(sm)
	
	return normalised_data,missing_lg

def normalise_rnaseq(sm_feat_data,lg_feat_data):
	normalised_data = dict()
	missing_lg = list()
	for sm in sm_feat_data:
		try:
			lg = sm.split(':')[0]
			nd = list()
			for i in xrange(len(sm_feat_data[sm])):
				try:
					nd += [sm_feat_data[sm][i]/lg_feat_data[lg][i]]
				except ZeroDivisionError:
					nd += ['NA']
			if nd.count(1) == 8: continue
			normalised_data[sm] = nd
			"""
			try:
				normalised_data[lg,sm] = array(sm_feat_data[sm])/array(lg_feat_data[lg])
				if sum(normalised_data[lg,sm]) == 8.0: del normalised_data[lg,sm]
			"""
		except KeyError:
			missing_lg.append(lg)
			continue
		
	return normalised_data,missing_lg

def correlations():
	return

def statistical_test(normalised_data,cascon,test_type,var_equal):
	no_all = len(cascon)
	cas = cascon[:no_all//2]
	con = cascon[no_all//2:]

	if test_type == "unpaired t-test":
		test = R.r['t.test']
		paired = True
		#test = stats.ttest_ind
	elif test_type == "unpaired wilcoxon":
		test = R.r['wilcox.test']
		paired = True
	
	test_results = dict()
	for n in normalised_data:
		if normalised_data[n].count(0) <= no_all/2 and normalised_data[n].count('NA') <= 0 and sum([normalised_data[n][i-3] for i in cas]+[normalised_data[n][j-3] for j in con]) != len(cascon):
			result = test(R.FloatVector([normalised_data[n][i-3] for i in cas]),R.FloatVector([normalised_data[n][j-3] for j in con]),paired=paired,var_equal=var_equal)
			test_results[n] = result[0][0],result[2][0]
	return test_results

def dict_to_file(mydict,outfile,t):
	f = open(outfile,'w')
	for i in mydict:
		#print i
		if t == 'p': print >> f,"\t".join([str(i[0]),str(i[1])]+map(str,mydict[i]))
		elif t == 'e': print >> f,"\t".join([i]+map(str,mydict[i]))
	f.close()
