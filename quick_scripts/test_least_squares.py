#!/home/paulk/software/bin/python
from __future__ import division
from sys import argv,exit,stderr
import scipy
import scipy.linalg as SLA

# generate random values for X and Y
X = scipy.transpose(scipy.mat(scipy.random.random((5,5))))
Y = scipy.transpose(scipy.mat(scipy.random.random(5)))

print "X = \n",X
print
print "Y = \n",Y
print
print "-----------------------------------------"
# use least squares to get B, the regression parameters
A = SLA.lstsq(X,Y)

print "scipy.lstsq: A = \n",A[0]
print
print "hermetian: A = \n",SLA.inv(scipy.conjugate(scipy.transpose(X))*X)*scipy.conjugate(scipy.transpose(X))*Y
print
print "non-hermetian: A = \n",SLA.inv(scipy.transpose(X)*X)*scipy.transpose(X)*Y
print
