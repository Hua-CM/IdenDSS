#!/bin/bash

WDIR=`pwd`

echo ${FILE}
python ${WDIR}\..\src\identification_dss.py database \
	-i $WDIR/example.fasta \
	-l 40 \
	-o $WDIR/example.fasta_db.fasta
python ${WDIR}\..\src\identification_dss.py iden \
	-l 40 \
	-m $WDIR/example_input.tsv \
	-d $WDIR/example.fasta_db.fasta \
	-@ 8 \
	-o $WDIR
