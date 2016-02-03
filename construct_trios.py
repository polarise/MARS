from __future__ import division
import sys
import gzip
import scipy
import scipy.stats
import rpy2.robjects as R
from rpy2.robjects.packages import importr
summary = R.r['summary']

stats = importr('stats')

class Gene( object ):
	"""
	gene class
	"""
	def __init__( self, ID ):
		"""
		gene constructor
		"""
		self.ID = ID
		self.mid_parent = list()
		self.offspring = list()
		self.h_sqrd = None
	
	def __repr__( self ):
		"""
		gene representation
		"""
		return "%s: %s" % ( self.ID, self.exprs )
		
	def add_expr( self, exprs, trios ):
		"""
		add the expressions only for trios
		"""
		for t in trios:
			T = trios[t] # a family
			if T.is_normal():
				self.mid_parent.append( ( exprs[T.parents[0]] + exprs[T.parents[0]] )/2 )
				self.offspring.append( exprs[T.offspring[0]] )
		return
		
	def heritability( self ):
		"""
		compute the regression coefficient (need to get the p-value)
		"""
		x = R.FloatVector( self.mid_parent )
		y = R.FloatVector( self.offspring )
		formula = R.Formula( 'y ~ x' )
		env = formula.environment
		env['x'] = x
		env['y'] = y
		results = stats.lm( formula )
		print summary( results )[3][1], summary( results )[3][7]
		r = results[0][1]
		self.h_sqrd = 2*r**2

class Trio( object ):
	"""
	trio class
	"""
	def __init__( self, name ):
		"""
		trio constructor
		"""
		self.name = name
		self.parents = list()
		self.offspring = list()
	
	def __repr__( self ):
		"""
		trio representation
		"""
		return "Parents: (%s); Child(ren): (%s)\nParent1 Expr: %s\nParent2 Expr: %s\nChild Expr: %s\n" % ( ", ".join( self.parents ), ", ".join( self.offspring ), self.parent1_expr, self.parent2_expr, self.child_expr )
	
	def add_parent( self, parent_id ):
		"""
		add parent id
		"""
		self.parents.append( parent_id )
		
	def add_child( self, child_id ):
		"""
		add child id
		"""
		self.offspring.append( child_id )
		
	def is_normal( self ):
		"""
		is this a literal trio i.e. three members only?
		"""
		return len( self.parents ) == 2 and len( self.offspring ) == 1
		
def process_ped( row ):
	"""
	function to process each row in the ped file
	"""
	L = row.strip().split( '\t' )
	name = L[5].split( ':' )[4].split( '.' )[0]
	ID = L[6].split( ':' )[4]
	return name, ID

def process_expr( row ):
	"""
	function to process each row in the expr file
	"""
	L = row.strip().split( '\t' )
	return L[0], map( float, L[1:] )

def is_parent( ID ): return ID in parent_IDs
def is_child( ID ): return ID in child_IDs


# make list of parent IDs
parent_IDs = [ row.strip() for row in open( "/home/paulk/MARS/A_offspring/parent_IDs.txt" ) ]

# make list of child IDs
child_IDs = [ row.strip() for row in open( "/home/paulk/MARS/A_offspring/offspring_IDs.txt" ) ]

# read the pedinfo file (cat the YRI and CEU data)
f = open( sys.argv[1] ) # pedinfo file
trios = dict()
for row in f:
	name, ID = process_ped( row )
	if name not in trios:
		trios[name] = Trio( name )
	if is_parent( ID ): trios[name].add_parent( ID ) # assumes a lot
	elif is_child( ID ): trios[name].add_child( ID )
f.close()
			
f = gzip.open( sys.argv[2] ) # expr
c = 0
sample_names = list()
genes = dict()
for row in f:
	if c > 10: break
	if row[0] == 'g':
		sample_names = row.strip().split( '\t' )[1:]
		continue
	gene_id, exprs = process_expr( row )
	genes[gene_id] = Gene( gene_id )
	genes[gene_id].add_expr( dict( zip( sample_names, exprs )), trios )
	genes[gene_id].heritability()
	print genes[gene_id].h_sqrd
	c += 1
f.close()

