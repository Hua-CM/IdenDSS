#!/bin/bash

WDIR=`pwd`

echo ${FILE}
python ${WDIR}/../IDSS.py database \
	-i ${WDIR}/example.fasta \
	-l 40 \
	-o ${WDIR}/example.fasta_db.fasta
python ${WDIR}/../IDSS.py iden \
	-l 40 \
	-m ${WDIR}/example_input.tsv \
	-d ${WDIR}/example.fasta_db.fasta \
	-@ 8 \
	-o $WDIR

python ${WDIR}/../IDSS.py plugin \
  -c \
  -i ${WDIR}/example.lst \
  -d ${WDIR}/example.fasta_db.fasta \
  -o ${WDIR}/RFLP
