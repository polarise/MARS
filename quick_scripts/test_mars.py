#!/home/paulk/software/bin/python
from __future__ import division
from sys import argv,exit,stderr
import rpy2.robjects as R
from rpy2.robjects.packages import importr

# import the mda package that contains the MARS implementation
mda = importr("mda")

"""
Script to read in the RNA-Seq expressions and the probes per gene
Columns in both files are samples
Columns are arranged similarly in both files
"""

# obtain the RNA-Seq data
f = open("../data/test_rnaseq.expr")
response = dict()
for line in f:
    l = line.strip().split('\t')
    response[l[0]] = {'y':R.FloatVector(map(float,l[1:]))}
f.close()

# read in the probe data; there are multiple probes per gene/transcript
f = open("../data/test_exonarray.expr")
predictors = dict()
for line in f:
    l = line.strip().split('\t')
    probes = [p.split(',') for p in l[1:]]
    for i in xrange(len(probes[0])): # from the lists of predictors pick the first, second, and so on from all the samples 
        if l[0] not in predictors: # do this so that making the dfs will be easy 
            predictors[l[0]] = {'x'+str(i):R.FloatVector([float(probes[j][i]) for j in xrange(len(l[1:]))])}
        else:
            predictors[l[0]]['x'+str(i)] = R.FloatVector([float(probes[j][i]) for j in xrange(len(l[1:]))])
f.close()

mars = R.r['mars']
predict = R.r['predict']

# print out the data for verification
for g in response: # for each gene
    df_response = R.DataFrame(response[g]) # create a data frame for the response
    df_predictors = R.DataFrame(predictors[g]) # create a df for the predictors
    rows = R.IntVector(range(1,7)) # indexes to be used to subset the df from rows 1 to 7
    fit = mars(df_predictors.rx(rows,True),df_response.rx(rows,True)) # build the mars  model
    print g,df_response.rx(7,True)[0],predict(fit,df_predictors.rx(7,True))[0] # predict
    print g,df_response.rx(8,True)[0],predict(fit,df_predictors.rx(8,True))[0] # predict
    #print df_response,fit[11]
    print
    #break
