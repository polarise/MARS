#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import key_functions

#with open( "/home/paulk/MARS/9A_additional_analyses/1_processed_data/HTS_samplenames.txt" ) as f:
#	fns = [ fn.strip() for fn in	f ]

fns = map( lambda x: str( x ) + "C", range( 3, 11 ))

txs_rpkm = dict()
for fn in fns:
#	f = open( "/home/paulk/MARS/9A_additional_analyses/1_processed_data/HTS/%s/flux_capacitor_output/Brain_gene_transcript_rpkm.gtf" % fn )
	f = open( "/home/paulk/MARS/F_GTEx/GTEx_RP_Data/HTS/%s/flux_capacitor_output/RP_gene_transcript_rpkm.gtf" % fn )
	
	c = 0
	for row in f:
		if c > 10: break
		pr = key_functions.process_feature( row )
		tx_id = pr['transcript_id'].split( "." )[0]
		gene_id = pr['gene_id'].split( "." )[0]
		tx_rpkm = pr['RPKM']
		
		if tx_id not in txs_rpkm:
			txs_rpkm[ tx_id ] = { "gene": gene_id, "rpkm": [ tx_rpkm ] }
		else:
			txs_rpkm[ tx_id ][ "rpkm" ] += [ tx_rpkm ]		
		c += 0

	f.close()
	
print( "tx_id", "gene_id", *fns, sep="\t" )
for tx in txs_rpkm:
	print( tx, txs_rpkm[tx]["gene"], *txs_rpkm[tx]["rpkm"], sep="\t" )



