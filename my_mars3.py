#!/home/paulk/software/bin/python
from __future__ import division
import sys
import scipy
import time
import itertools
import Orange
import multiprocessing

"""
serial version
"""

def current_line(f,g):
	while f and g:
		yield f.readline(),g.readline()
		
def line_batch(f,g,M=10):
	while f and g:
		yield ((f.readline(),g.readline()) for i in xrange(M))

def Seekize(row,q=None):
	f_line,g_line = row
	g,no_samps,no_probes,probe_IDs,fpkm,probes = f_line.strip().split('\t')

	no_samps = int(no_samps)
	no_probes = int(no_probes)
	fpkm = map(float,fpkm.split(','))
#	if fpkm.count(0.0) > 0.95*len(fpkm):
#		continue #print >> sys.stderr,"Warning: too many zeros!\n"; continue
	probes = map(float,probes.split(','))

	# required Orange constructs
	predictors = [Orange.feature.Continuous(name) for name in J[:no_probes]]
	response = Orange.feature.Continuous('RS')
	domain = Orange.data.Domain(predictors + [response])

	probes_matrix = [[scipy.log(probes[no_samps*j + i]) for j in xrange(no_probes)] + [scipy.log(fpkm[i]+1e-20)] for i in xrange(no_samps)]
	data = Orange.data.Table(domain,probes_matrix)
	
	g,no_samps,no_probes,probe_IDs,fpkm,probes = g_line.strip().split('\t')
	no_samps = int(no_samps)
	no_probes = int(no_probes)
	fpkm = map(float,fpkm.split(','))
	probes = map(float,probes.split(','))	

	try:
		M = Orange.regression.earth.EarthLearner(data,degree=1,terms=None)
	except ValueError:
		M = None
	
	if M == None:
		print >> sys.stderr,"[%s]\tInput data error: variance in response in zero" % time.ctime(time.time())
		return g,fpkm,["NA"]*no_samps,"NA","NA"
	else:
		probes_matrix = [[scipy.log(probes[no_samps*j + i]) for j in xrange(no_probes)] for i in xrange(no_samps)]
		domain = Orange.data.Domain(predictors)
		instances = [Orange.data.Instance(domain,p) for p in probes_matrix]
		try:
			M_predicted = map(lambda x:round(scipy.exp(M.predict(x)[0]),5),instances)
		except RuntimeWarning:
			print >> sys.stderr,"[%s]\tWarning: Failed to make a prediction because of something" % time.ctime(time.time())
			M_predicted = None
		if M_predicted != None:
			R,p = scipy.stats.mstats.spearmanr(fpkm,M_predicted,use_ties=True)
			return g,fpkm,M_predicted,R,p
#			q.put((g,fpkm,M_predicted,R,p))
		else:
			return g,fpkm,"NA"*no_samps,"NA","NA"
#			q.put((g,fpkm,"NA"*no_samps,"NA","NA"))

try:
	f_fn = sys.argv[1]
	g_fn = sys.argv[2]
except IndexError:
	print >> sys.stderr,"Usage:./script.py <train> <test>"
	sys.exit(1)

# prepare the predictor names
I = itertools.product('abcdefghijklmnopqrstuvwxyzABCDEF',repeat=2)
J = ["".join(i) for i in I]
	
F = open(f_fn)
G = open(g_fn)

#c = 0
#for f_line,g_line in current_line(F,G):
#	if c > 20: break
#	g,fpkm,M_predicted,R,p = Seekize(f_line,g_line)
#	print "\t".join(map(str,[g,",".join(map(str,fpkm)),",".join(map(str,M_predicted)),R,p]))
#	c += 1

def printer(q,child):
	while child.recv() != 'started':
		pass
	print 'Child recvd \'started\' msg'
	done = False
	while not q.empty() or not done:
		g,fpkm,M_predicted,R,p = q.get()
		print "\t".join(map(str,[g,",".join(map(str,fpkm)),",".join(map(str,M_predicted)),R,p]))
		if not done:
			try:
				signal = child.recv()
			except IOError:
				signal = ''
			if signal == 'done':
				print 'Child recvd \'done\' msg'
				done = True


#parent,child = multiprocessing.Pipe()
#queue = multiprocessing.Queue()
#q_pro = multiprocessing.Process(target=printer,args=(queue,child,))
#q_pro.start()

pro_list = list()
c = 0
#parent.send('started')
for row in itertools.izip(F,G):
	if c > 10: break
	g,fpkm,M_predicted,R,p = Seekize(row)
	print "\t".join(map(str,[g,",".join(map(str,fpkm)),",".join(map(str,M_predicted)),R,p]))
	c += 0
#parent.send('done')
#sys.exit(0)

F.close()
G.close()

