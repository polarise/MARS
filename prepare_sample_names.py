#!/usr/bin/env python
from __future__ import division
import sys
import random

def pick_samples( names ):
  train = list()
  test = list()
  
  while len( train ) < 87: ## WATCH OUT!!!
    pick = random.choice( names )
    train.append( pick )
    names.remove( pick )
  
  test = names
  
  return train, test  
  
def make_row( name, train=True, test_present=True, mapped=False ):
	if not mapped:
		if train and test_present:
		  print "%s\t%s.CEL" % ( name, name )
		elif not train and test_present:
		  print "*%s\t%s.CEL" % ( name, name )
		elif not train and not test_present:
		  print "*NA\t%s.CEL" % ( name )
		else:
		  raise ValueError( "Hmmm... I'm not designed for this. What should I do, my lord?" )
		  return
	elif mapped:
		if train and test_present:
		  print "%s\t%s" % ( name, mapped_names[ name ] )
		elif not train and test_present:
		  print "*%s\t%s" % ( name, mapped_names[ name ] )
		elif not train and not test_present:
		  print "*NA\t%s" % ( mapped_names[ name ] )
		else:
		  raise ValueError( "Hmmm... I'm not designed for this. What should I do, my lord?" )		
  		return

f = open( "/home/paulk/MARS/0_prep_data/combined_reduced_ALL.txt.ordered" )
#f = open( sys.argv[1] )
names = [ row.strip().split( '\t' )[0] for row in f ]
f.close()

#f = open( sys.argv[2] )
#mapped_names = {}
#for row in f:
#	L = row.strip().split( '\t' )
#	mapped_names[L[0]] = L[1]
#f.close()

train, test = pick_samples( names )

print "hts\tma"
for s in train:
  make_row( s )#, mapped=True )
for s in test:
  make_row( s, train=False )#, mapped=True )
