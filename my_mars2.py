#!/home/paulk/software/bin/python
from __future__ import division
import sys
import scipy
import time
import itertools
import Orange
import multiprocessing
import Queue

"""
parallel version
"""
def remove_uncorrelated(p_matrix,min_cor=0.5):
	p_matrix = scipy.array(p_matrix)	# make an array
	fpkm = p_matrix[:,-1]							# get the last column of fpkm values
	base_matrix = p_matrix[:,:-1]			# get the rest of matrix
	new_p_matrix = []									# the list of accepted columns
	selected = []											# record of what is accepted (1) /rejected (0)
	for i in xrange(scipy.shape(base_matrix)[1]):
		current_column = base_matrix[:,i]
		cor = scipy.stats.mstats.spearmanr(fpkm,current_column)[0]
		if cor >= min_cor:
			new_p_matrix.append(list(current_column))
			selected.append(1)
		else:
			selected.append(0)

	if len(new_p_matrix) == 0:
		return [],selected
	else:
		new_p_matrix.append(list(fpkm))	# void function (duh!)
		new_p_matrix = list(map(list,scipy.array(new_p_matrix).transpose()))
	
		return new_p_matrix,selected
	
def subset_to_selected(p_matrix,selected):
	p_matrix = scipy.array(p_matrix)
	new_p_matrix = []
	for i in xrange(scipy.shape(p_matrix)[1]):
		if selected[i]:
			new_p_matrix.append(list(p_matrix[:,i]))
	new_p_matrix = list(map(list,scipy.array(new_p_matrix).transpose()))
		
	return new_p_matrix

def Seekize(row,deg,min_cor,q=None):
	f_line,g_line = row
	
	# training
	g_tr,no_samps_tr,no_probes_tr,probe_IDs_tr,fpkm_tr,probes_tr = f_line.strip().split('\t')

	no_samps_tr = int(no_samps_tr)
	no_probes_tr = int(no_probes_tr)
	fpkm_tr = map(float,fpkm_tr.split(','))
	probes_tr = map(float,probes_tr.split(','))


	probes_matrix_tr = [[probes_tr[no_samps_tr*j + i] for j in xrange(no_probes_tr)] + [fpkm_tr[i]] for i in xrange(no_samps_tr)]
	probes_matrix_tr,selected = remove_uncorrelated(probes_matrix_tr,min_cor)
	if probes_matrix_tr == []:
#		print >> sys.stderr,"[%s]\tError: Failed to find correlated predictors" % time.ctime(time.time())
		q.put((g_tr,fpkm_tr,["NA"]*no_samps_tr,"NA","NA","NA","NA"))
	else:
		# required Orange constructs :: generic (used for training and testing)
		predictors = [Orange.feature.Continuous(name) for name in J[:sum(selected)]]
		response = Orange.feature.Continuous('RS')
	
		# training domain object
		domain_tr = Orange.data.Domain(predictors + [response])
		
		data = Orange.data.Table(domain_tr,probes_matrix_tr)
	
		# testing
	
		g_te,no_samps_te,no_probes_te,probe_IDs_te,fpkm_te,probes_te = g_line.strip().split('\t')
		no_samps_te = int(no_samps_te)
		no_probes_te = int(no_probes_te)
		fpkm_te = map(float,fpkm_te.split(','))
		probes_te = map(float,probes_te.split(','))	

		# now build the model object
		try:
			M = Orange.regression.earth.EarthLearner(data,degree=deg,terms=None)
		except ValueError:
			M = None
	
		if M == None:
#			print >> sys.stderr,"[%s]\tInput data error: variance in response is zero" % time.ctime(time.time())
			q.put((g_te,fpkm_te,["NA"]*no_samps_te,"NA","NA","NA","NA"))
		else:
			probes_matrix_te = [[probes_te[no_samps_te*j + i] for j in xrange(no_probes_te)] for i in xrange(no_samps_te)]
			probes_matrix_te = subset_to_selected(probes_matrix_te,selected)
			domain_te = Orange.data.Domain(predictors)
		
			# build instance objects for each set of predictor values
			instances = [Orange.data.Instance(domain_te,p) for p in probes_matrix_te]
			M_predicted = map(lambda x:round(M.predict(x)[0],5),instances)
			#M_predicted = map(lambda x: 0 if x < 0 else x,M_predicted)
	#		R,pval = scipy.stats.mstats.spearmanr(fpkm_te,M_predicted,use_ties=True)
			R,Rpval = scipy.stats.mstats.pearsonr(fpkm_te,M_predicted)
			T,Tpval = scipy.stats.ttest_ind(fpkm_te,M_predicted)
	#		print >> sys.stderr,M_predicted,R,pval
			q.put((g_te,fpkm_te,M_predicted,R,Rpval,T,Tpval))

def printer(q,child):
	# wait until data starts arriving
	time.sleep(2)
	msg = ''
	msg = child.recv()
	while msg != 'started':
		try:
			msg = child.recv()
		except IOError:
			continue
	
	# now I know that data has started arriving
#	print >> sys.stderr,"[%s]\tChild recvd \'started\' msg" % time.ctime(time.time())
	
	# assmue we're not done yet and wait
	done = False
	wait = True
	while not q.empty() or not done or not wait:
#		print >> sys.stderr,q.empty(),done,wait
		g,fpkm,M_predicted,R,Rp,T,Tp = q.get()
		print "\t".join(map(str,[g,",".join(map(str,fpkm)),",".join(map(str,M_predicted)),R,Rp,T,Tp]))
		if not done:
			try:
				signal = child.recv()
			except IOError:
				signal = ''
			if signal == 'done':
#				print >> sys.stderr,"[%s]\tChild recvd \'done\' msg" % time.ctime(time.time())
				done = True
		if q.empty():
			wait = False

try:
	deg = sys.argv[1]
	min_cor = sys.argv[2]
	f_fn = sys.argv[3]
	g_fn = sys.argv[4]
except IndexError:
	print >> sys.stderr,"Usage:./script.py <train> <test>"
	sys.exit(1)
	
deg = int(deg)
min_cor = float(min_cor)

# prepare the predictor names
I = itertools.product('abcdefghijklmnopqrstuvwxyzABCDEF',repeat=2)
J = ["".join(i) for i in I]

# the files
F = open(f_fn)
G = open(g_fn)

parent,child = multiprocessing.Pipe()
queue = multiprocessing.Queue()
q_pro = multiprocessing.Process(target=printer,args=(queue,child,))
q_pro.start()
#q_pro.join()

unwanted = [' ']

pro_list = list()
c = 0
parent.send('started')
for row in itertools.izip(F,G):
	if len(row[0].strip().split('\t')) != 6: continue # lines that are not fully formatted
	if c > 10: break
	while len(pro_list) > 15:
		for p in pro_list:
			if not p.is_alive(): pro_list.remove(p)
	pro = multiprocessing.Process(target=Seekize,args=(row,deg,min_cor,queue,))
	pro_list.append(pro)
	pro.start()
	c += 0
parent.send('done')

queue.join_thread()

F.close()
G.close()


	"""def EARTH(self,row1,row2):
		"""
#		the EARTH function as implemented by <SOMEONE>
		"""
		# get the training data ('*_tr': training data; '*_te': testing data)
		l1 = row1.strip().split('\t')
		gene_name = l1[0]
		no_samples_tr = int(l1[1])
		no_probes_tr = int(l1[2])
		fpkm_tr = map(float,l1[4].split(','))

		# get the testing data
		l2 = row2.strip().split('\t')
		no_samples_te = int(l2[1])
		no_probes_te = int(l2[2])
		fpkm_te = map(float,l2[4].split(','))
	
		try:	# check whether you can do any prediction
			K = len(l1[5].split(','))
		except IndexError: # we can't; push zeros
			return "\t".join(map(str,[gene_name,no_probes_te,'NA'] + [l2[4]] + [",".join(['0.0']*no_samples_te)] + ['NA']*6))		# OUTPUT
			
		probes_tr = map(float,l1[5].split(','))

		probes_te = map(float,l2[5].split(','))

		probes_matrix_tr = [[probes_tr[no_samps_tr*j + i] for j in xrange(no_probes_tr)] + [fpkm_tr[i]] for i in xrange(no_samps_tr)]
		probes_matrix_tr,selected = remove_uncorrelated(probes_matrix_tr,min_cor)
		if probes_matrix_tr == []:
	#		print >> sys.stderr,"[%s]\tError: Failed to find correlated predictors" % time.ctime(time.time())
			q.put((g_tr,fpkm_tr,["NA"]*no_samps_tr,"NA","NA","NA","NA"))
		else:
			# required Orange constructs :: generic (used for training and testing)
			predictors = [Orange.feature.Continuous(name) for name in J[:sum(selected)]]
			response = Orange.feature.Continuous('RS')
	
			# training domain object
			domain_tr = Orange.data.Domain(predictors + [response])
		
			data = Orange.data.Table(domain_tr,probes_matrix_tr)
	
			# testing
	
			g_te,no_samps_te,no_probes_te,probe_IDs_te,fpkm_te,probes_te = g_line.strip().split('\t')
			no_samps_te = int(no_samps_te)
			no_probes_te = int(no_probes_te)
			fpkm_te = map(float,fpkm_te.split(','))
			probes_te = map(float,probes_te.split(','))	

			# now build the model object
			try:
				M = Orange.regression.earth.EarthLearner(data)
			except ValueError:
				M = None
	
			if M == None:
	#			print >> sys.stderr,"[%s]\tInput data error: variance in response is zero" % time.ctime(time.time())
				q.put((g_te,fpkm_te,["NA"]*no_samps_te,"NA","NA","NA","NA"))
			else:
				probes_matrix_te = [[probes_te[no_samps_te*j + i] for j in xrange(no_probes_te)] for i in xrange(no_samps_te)]
				probes_matrix_te = subset_to_selected(probes_matrix_te,selected)
				domain_te = Orange.data.Domain(predictors)
		
				# build instance objects for each set of predictor values
				instances = [Orange.data.Instance(domain_te,p) for p in probes_matrix_te]
				M_predicted = map(lambda x:round(M.predict(x)[0],5),instances)
				#M_predicted = map(lambda x: 0 if x < 0 else x,M_predicted)
		#		R,pval = scipy.stats.mstats.spearmanr(fpkm_te,M_predicted,use_ties=True)
				R,Rpval = scipy.stats.mstats.pearsonr(fpkm_te,M_predicted)
				T,Tpval = scipy.stats.ttest_ind(fpkm_te,M_predicted)
		#		print >> sys.stderr,M_predicted,R,pval
				q.put((g_te,fpkm_te,M_predicted,R,Rpval,T,Tpval))"""
