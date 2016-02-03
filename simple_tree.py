#!/home/paulk/software/bin/python
import scipy
import scipy.stats
import random

def RSS(y_L,y_R):
	y_L_star = scipy.mean(y_L)
	y_R_star = scipy.mean(y_R)
	RSS_L = 0
	for i in y_L:
		RSS_L += (i - y_L_star)**2
	RSS_R = 0
	for i in y_R:
		RSS_R += (i - y_R_star)**2
	
	return RSS_L + RSS_R

# find the best x
def best_x(x_y): # x_y is a list of tuples ((x_1,y_1),...,(x_n,y_n))
	res = list()
	for x,y in x_y:
		y_L = [j[1] for j in x_y if j[0] <= x]
		y_R = [j[1] for j in x_y if j[0] > x]
		rss = RSS(y_L,y_R)
		res.append((x,rss))
		print "RSS at (%.5f,%.5f)\t=\t%.5f" % (x,y,rss)
	
	# get the x at which the minimum RSS occurs
	x_B = scipy.Inf
	rss_B = scipy.Inf
	for x,rss in res:
		if rss < rss_B:
			x_B = x
			rss_B = rss
	return x_B,rss_B

def split_data(x_y,x_B):
	x_y_L = [j for j in x_y if j[0] <= x_B]
	x_y_R = [j for j in x_y if j[0] > x_B]
	return x_y_L,x_y_R
	

def build_tree(x_y):
	x_B,rss_B = best_x(x_y)
	print "The best value of x is %.5f (RSS = %.5f)." % (x_B,rss_B)
	x_y_L,x_y_R = split_data(x_y,x_B)
	print "Left and right lengths: %d,%d" % (len(x_y_L),len(x_y_R))
	print
	
	x_BL,rss_BL = best_x(x_y_L)
	print "The best value of x_L is %.5f (RSS = %.5f)." % (x_BL,rss_BL)
	x_y_LL,x_y_LR = split_data(x_y_L,x_BL)
	print "Left and right lengths: %d,%d" % (len(x_y_LL),len(x_y_LR))
	print
	
	x_BR,rss_BR = best_x(x_y_R)
	print "The best value of x_R is %.5f (RSS = %.5f)." % (x_BR,rss_BR)
	x_y_RL,x_y_RR = split_data(x_y_R,x_BR)
	print "Left and right lengths: %d,%d" % (len(x_y_RL),len(x_y_RR))
	print

# generate x and y
x = scipy.stats.poisson.rvs(10,size=100)
x.sort()
x_y = [(i,scipy.stats.poisson.rvs(5)) for i in x]
#x_B,rss_B = best_x(x_y)

build_tree(x_y)

#print "The best value of x is %.5f (RSS = %.5f)" % (x_B,rss_B)
