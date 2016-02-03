#!/home/paulk/software/bin/python
from __future__ import division
import sys
import rpy2.robjects as R
from rpy2.robjects.packages import importr
import scipy
import scipy.stats
import itertools
import multiprocessing
import time
import operator

# R-related imports and functions
importr('mda')
importr('utils')
importr('rpart')
mars = R.r['mars']
predict = R.r['predict']
matrix = R.r['matrix']
cor_test = R.r['cor.test']
cor = R.r['cor']
t_test = R.r['t.test']
options = R.r['options']
Rsum = R.r['sum']
Rlog = R.r['log']
Rsqrt = R.r['sqrt']
rpart = R.r['rpart']
t = R.r['t']
dim = R.r['dim']

options(warn=-1)

class Worker(multiprocessing.Process):
	"""
	The Working Class ;-)
	"""
	def __init__(self,i,q1,q2,min_cor):
		multiprocessing.Process.__init__(self)
		self.proc_no = i
		self.q1 = q1
		self.q2 = q2
		self.min_cor = min_cor
		
	def return_medians(self,probes_tr,no_samples_tr,no_probes_tr):
#		grouped_probes_tr = [probes_tr[i*no_samples_tr : (i + 1)*no_samples_tr] for i in xrange(no_probes_tr)]		# group the probes by probe
		grouped_probes_tr = [[probes_tr[no_samples_tr*j + i] for j in xrange(no_probes_tr)] for i in xrange(no_samples_tr)]
		medians_tr = map(scipy.median,grouped_probes_tr)
		return medians_tr,[1]
	
	def my_mode(self,values):
		R = values
		r_min = min(R)
		r_max = max(R)
		w = 10.
		k = (r_max - r_min)/w
		rint = []
		i = 0
		while i < k:
			rint += [(r_min + i*w,r_min + (i+1)*w)]
			i += 1
			
		bins = {ri:0 for ri in rint}
		
		for r in R:
			done = False
			for ri in bins:
				if ri[0] <= r < ri[1]:
					bins[ri] += 1
					done = True
				if done:
					break
		
		rep = max(bins.iteritems(),key=operator.itemgetter(1))
		return scipy.mean(rep[0])

	def return_pseudo_mode(self,probes_tr,no_samples_tr,no_probes_tr):
		grouped_probes_tr = [[probes_tr[no_samples_tr*j + i] for j in xrange(no_probes_tr)] for i in xrange(no_samples_tr)]		# group the probes by probe
		modes_tr = map(self.my_mode,grouped_probes_tr)
		return modes_tr,[1]
	
	def remove_uncorrelated(self,probes_tr,fpkm_tr,no_samples_tr,no_probes_tr,ignore):
		"""
		remove all predictors that have Spearman correlation less than min_cor with the response
		"""
		if ignore or sum(fpkm_tr) == 0:
			return probes_tr,[1]*no_probes_tr
		# make an array
		grouped_probes_tr = [probes_tr[i*no_samples_tr : (i + 1)*no_samples_tr] for i in xrange(no_probes_tr)]		# group the probes by probe
		correlations = [scipy.stats.mstats.spearmanr(fpkm_tr,g,use_ties=True)[0] for g in grouped_probes_tr]
		selected_probes_tr = [grouped_probes_tr[i] for i in xrange(len(grouped_probes_tr)) if correlations[i] >= self.min_cor]
		selected = [1 if correlations[i] >= self.min_cor else 0 for i in xrange(len(grouped_probes_tr))]
		filt_probes_tr = [item for sublist in selected_probes_tr for item in sublist]
		return filt_probes_tr,selected
				
	def subset_to_selected(self,probes_te,selected,no_samples_te,no_probes_te,ignore):
		if ignore:
			return probes_te
		grouped_probes_te = [probes_te[i*no_samples_te : (i + 1)*no_samples_te] for i in xrange(no_probes_te)]		# group the probes by probe
		selected_probes_te = [grouped_probes_te[i] for i in xrange(len(grouped_probes_te)) if selected[i]]
		filt_probes_te = [item for sublist in selected_probes_te for item in sublist]
		return filt_probes_te
	
	def MARS(self,row1,row2):
		"""
		the MARS function as implemented by Hastie and Tibshirani
		"""
		# get the training data ('*_tr': training data; '*_te': testing data)
		l1 = row1.strip().split('\t')
		gene_name = l1[0]
		no_samples_tr = int(l1[1])
		no_probes_tr = int(l1[2])
		fpkm_tr = R.FloatVector(map(lambda x:x,map(float,l1[4].split(','))))

		# get the testing data
		l2 = row2.strip().split('\t')
		no_samples_te = int(l2[1])
		no_probes_te = int(l2[2])
		fpkm_te = R.FloatVector(map(float,l2[4].split(',')))
		
		try:	# check whether you can do any prediction
			K = len(l1[5].split(','))
		except IndexError: # we can't; push zeros
			return "\t".join(map(str,[gene_name,no_probes_te,'NA'] + [l2[4]] + [",".join(['0.0']*no_samples_te)] + ['NA']*6))		# OUTPUT
			
		if l1[5] == 'NA':
			return "\t".join(map(str,[gene_name,no_probes_te,'NA'] + [l2[4]] + [",".join(['0.0']*no_samples_te)] + ['NA']*6))		# OUTPUT
		
		# we can! make the probe matrix for prediction
		# filter probes by correlation
		ignore = True			# whether to ignore the min_cor value or not
		filt_probes_tr,selected = self.remove_uncorrelated(map(float,l1[5].split(',')),map(float,l1[4].split(',')),no_samples_tr,no_probes_tr,ignore)
#		filt_probes_tr,selected = self.return_pseudo_mode(map(float,l1[5].split(',')),no_samples_te,no_probes_te)					# using median of each
		no_filt_probes_tr = sum(selected)

		if len(selected) == 0 or not no_filt_probes_tr:
			return "\t".join(map(str,[gene_name,no_probes_te,0] + [l2[4]] + [",".join(['0.0']*no_samples_te)] + ['NA']*6))		# OUTPUT
		
		probe_matrix_tr = matrix(R.FloatVector(filt_probes_tr),nrow=no_samples_tr,ncol=no_filt_probes_tr)

		# build the model
		model = mars(probe_matrix_tr,fpkm_tr,nk=600,prune=False)

		# construct the predictor matrix
		filt_probes_te = self.subset_to_selected(map(float,l2[5].split(',')),selected,no_samples_te,no_probes_te,ignore)
#		filt_probes_te,selected = self.return_pseudo_mode(map(float,l2[5].split(',')),no_samples_te,no_probes_te)					# using median of each
 		
 		probe_matrix_te = matrix(R.FloatVector(filt_probes_te),nrow=no_samples_te,ncol=no_filt_probes_tr)
		
		# make a prediction
		predicted_fpkm = predict(model,probe_matrix_te)
		
		predicted_fpkm = map(lambda x:x if x > 0.0 else 0.0,predicted_fpkm)		# negative->zero
#		predicted_fpkm = map(lambda x:x,predicted_fpkm)													# leave negative as-is
		rho,rho_pval = scipy.stats.mstats.pearsonr(map(float,l2[4].split(',')),map(lambda x:x,predicted_fpkm))
		t_stat,pval = scipy.stats.ttest_ind(map(float,l2[4].split(',')),map(lambda x:x,predicted_fpkm))
		mean_te = scipy.mean(map(float,l2[4].split(',')))
		var_pr = scipy.var(predicted_fpkm)
		
		if scipy.isnan(rho):
			rho,rho_pval = 'NA','NA'
		
		if scipy.isnan(t_stat):
			t_stat,pval = 'NA','NA'
		
		# return the results	OUTPUT MODIFIED!!!
		return "\t".join(map(str,[gene_name,no_probes_te,no_probes_tr] + [l2[4]] + [",".join(map(lambda x:str(x),predicted_fpkm))] + [rho,rho_pval,t_stat,pval,mean_te,var_pr]))	# don't suppress negatives
	
	def run(self):
		print >> sys.stderr,"[%s]\tProcess %s starting." % (time.ctime(time.time()),self.proc_no)
		while True:
			row1,row2 = self.q1.get()
			
			if row1 == None and row2 == None: # is this the poison pill?
				print >> sys.stderr,"[%s]\tProcess %s exiting normally." % (time.ctime(time.time()),self.proc_no)
				self.q2.put('END')	# send a poison pill down the output queue as well
				
				self.q1.task_done()
				self.q2.task_done()
				break
			
			# change the prediction method as required
			result = self.MARS(row1,row2)
			self.q2.put(result)
			
			self.q1.task_done()			# this was the culprit!!! remember to signal "task_done" if you have to use 'continue'
			self.q2.task_done()
		return

#-------------------------------------------------------------------------------
# program begins here...
#-------------------------------------------------------------------------------
if __name__ == '__main__':
	try:
		train_fn = sys.argv[1]
		test_fn = sys.argv[2]
		min_cor = float(sys.argv[3])
	except IndexError:
		print >> sys.stderr,"Usage: script.py <train_fn> <test_fn> <min_cor>"
		sys.exit(0)
		
	try:
		assert -1.0 <= min_cor <= 1.0
	except:
		print >> sys.stderr,"Warning: invalid 'min_cor'; will use min_cor = -1"
		min_cor = -1.0
	
	num_of_processes = multiprocessing.cpu_count()

	# the queues...
	q1 = multiprocessing.JoinableQueue()	# read queue
	q2 = multiprocessing.JoinableQueue()	# write queue

	# read in the data
	f = open(train_fn)
	g = open(test_fn)
	c = 0
	for row1,row2 in itertools.izip(f,g):
		if c > 1000: break
		q1.put((row1,row2))
		c += 1
	f.close()
	g.close()
	
	# poison pill
	for i in xrange(num_of_processes):
		q1.put((None,None))

	# start the worker processes
	procs = [Worker(i,q1,q2,min_cor) for i in xrange(num_of_processes)]
	for p in procs:
		p.daemon = True
		p.start()
	
	# wait until the read queue is done with
	q1.join()
	print >> sys.stderr,"[%s]\tInput queue joined successfully." % time.ctime(time.time())

	# get the output
	end = list()
	while len(end) < num_of_processes:
		output = q2.get()
		if output != 'END':
			print output
		else:
			end.append(output)
	print >> sys.stderr,"[%s]\tOutput queue emptied successfully." % time.ctime(time.time())
	
	# join the write queue
	q2.join()
	print >> sys.stderr,"[%s]\tOutput queue joined successfully." % time.ctime(time.time())
