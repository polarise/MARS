# launch N = 10 procesess
# how does each process know that the work hasn't begun?
# how does each process obtain data?
# how does each process process data?
# how does each process know when the job is done?


# go through the file and hand M = 1000 lines to each process in turn

import multiprocessing
import time
import random

# construct big generators
def make_gen(f_obj,M):
	while f_obj:
		yield (f_obj.readline() for i in xrange(M))

# works!
#f = open("together.alla")
#c = 0
#for X in make_gen(f,1000):
#	if c > 5: break
#	for row in X:
#		print row.strip()
#	print
#	c += 1
	
# define the compute function
def fxn(conn):
	"""
	data is a generator of M lines or a single line
	"""
	data = None
	while data is None:
		data = conn.recv()
	result = map(lambda x:x**2,data)
	conn.close()

# parallel version
t = time.time()
pro_list = list()
par1,child1 = multiprocessing.Pipe()
# generate random data 100 times
for i in xrange(1000):
	data = [random.random() for i in xrange(10)]
	
	# use at most 10 processes
	while len(pro_list) > 10:
		for p in pro_list:
			if not p.is_alive(): pro_list.remove(p)
	p = multiprocessing.Process(target=fxn,args=(child1,))
	pro_list.append(p)
	p.start()
	par1.send(data)
par1.close()
print "Parallel: %s" % (time.time() - t)
		
t = time.time()
# serial version
for i in xrange(1000):
	data = [random.random() for i in xrange(10)]
	result = map(lambda x:x**2,data)
print "Serial: %s" % (time.time() - t)



# create a new write process 
# open a pipe for this process 

# open the file
# pass data to the compute process through the pipe


# pass the results from the compute process to the write process through the pipe



