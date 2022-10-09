#!/bin/bash

IdenDSS index \
	-i example.fasta \
	-o database/example_db.fasta \
	-l 40

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
