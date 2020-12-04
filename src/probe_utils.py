# -*- coding: utf-8 -*-
# @Time : 2020/12/4 13:01
# @Author : Zhongyi Hua
# @FileName: probe_utils.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from functools import reduce
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline


def probe_generate_c(_seq_obj, _length=20):
    """

    :param _seq_obj:
    :param _length:
    :return:
    """
    sequence = _seq_obj.seq + _seq_obj.seq[0:_length]
    _seq_obj.id = _seq_obj.id
    assembly_dict = {}
    for index in range(len(sequence) - (_length - 1)):
        assembly_dict.setdefault(str(sequence[index: index + _length]), _seq_obj.id + '_' + str(index + 1))
    return assembly_dict


def probe_generate_l(_seq_obj, _length=20):
    """

    :param _seq_obj:
    :param _length:
    :return:
    """
    _seq_obj.id = _seq_obj.id
    assembly_dict = {}
    for index in range(len(_seq_obj.seq) - (_length - 1)):
        assembly_dict.setdefault(str(_seq_obj.seq[index: index + _length]), _seq_obj.id + '_' + str(index + 1))
    return assembly_dict


class BLASTParse:
    def __init__(self, _seq_path, _asm_list, _database_path, _tmp_dir_path, _length, _thread, _exec=''):
        self.seq_path = _seq_path
        self.asm_list = _asm_list
        self.db_path = _database_path
        self.tmp_dir = _tmp_dir_path
        self.length = _length
        self.thread = _thread
        self.exec = _exec

    def blast_cmd(self):
        _group = os.path.splitext(os.path.split(self.seq_path)[1])[0]
        blastn_cmd = NcbiblastnCommandline(
            cmd=os.path.join(self.exec, 'blastn'),
            query=self.seq_path,
            db=self.db_path,
            task='blastn-short',
            dust='no',
            word_size=min(round(self.length/2)-1, 11),
            outfmt='\"6 qacc sacc length pident evalue\"',
            num_threads=self.thread,
            evalue=10,
            out=os.path.join(self.tmp_dir, _group + '.txt'),
            max_hsps=1,
            max_target_seqs=len(self.asm_list)+1)
        blastn_cmd()

    def blast_parse(self):
        _putative_kmer = []
        _group = os.path.splitext(os.path.split(self.seq_path)[1])[0]
        _blast_df = pd.read_table(os.path.join(self.tmp_dir, _group + '.txt'),
                                  names=['query', 'subject', 'length', 'identity', 'evalue'])
        _blast_df = _blast_df[(_blast_df['length'] == self.length) & (_blast_df['identity'] == 100)]
        # Find DSS **in** other species
        _blast_df1 = _blast_df.loc[~(_blast_df['subject'].isin(self.asm_list))].drop_duplicates(subset=['query'])
        # Find DSS conserved in species
        _blast_df2 = _blast_df[_blast_df['subject'].isin(self.asm_list)]
        _blast_df2 = _blast_df2.groupby('query')['subject'].count().reset_index()
        _blast_df2 = _blast_df2[_blast_df2['subject'] == len(self.asm_list)]
        # Find intersection as species DSS
        _dss_list = list(set(_blast_df2['query'].to_list()) - set(_blast_df1['query'].to_list()))
        # DSS set empty
        if not _dss_list:
            _result_tb = pd.DataFrame.from_dict({0: {'seq':  '', 'position': '', 'GC': ''}}, 'index')
        else:
            _kmer_dict = {_.id: _.seq for _ in
                          SeqIO.parse(self.seq_path, 'fasta')
                          }
            _tmp_list = []
            for _dss in _dss_list:
                _tmp_list.append(
                                 {'assembly': '_'.join('NCABA_2_1'.split('_')[:-1]),
                                  'seq': _kmer_dict[_dss],
                                  'position': _dss.split('_')[-1] +
                                              '-' +
                                              str(int(_dss.split('_')[-1])+self.length-1),
                                  'GC': GC(_kmer_dict[_dss])}
                                 )
            _result_tb = pd.DataFrame(_tmp_list)
        return _result_tb


def iden_main(args):
    if args.circular:
        probe_generate = probe_generate_c
    else:
        probe_generate = probe_generate_l
    seq_dict = {_.id: _ for _ in SeqIO.parse(args.database, 'fasta')}
    _table = pd.read_table(args.meta, names=['group', 'sample', 'assembly'], dtype=str)
    # by group for generate probe
    _group_table = _table.groupby(['group'])['assembly'].apply(list).reset_index(name='assembly')
    # by sample for pre-BLAST
    _group_sample_table = _table.groupby(['group', 'sample'])['assembly'].apply(list).reset_index(name='assembly')
    for _idx, _row in _group_table.iterrows():
        # pre-probe
        _super_probe = defaultdict(str)
        for _asm in _row['assembly']:
            for _k, _v in probe_generate(seq_dict[_asm], args.length):
                _super_probe[_k] = _v
        _probe_lines = ['>' + seq_id + '\n' + seq for seq, seq_id in _super_probe.items()]
        with open(os.path.join(args.tmp, 'pre_probe.fasta'), 'w') as f:
            f.write('\n'.join(_probe_lines))
        del _super_probe, _probe_lines
        # pre-blast
        ## make blast database
        _db_lines = ['>' + seq_id + '\n' + seq_dict[seq_id].seq for seq_id in _row['assembly']]
        with open(os.path.join(args.tmp, 'pre_blast_db.fasta'), 'w') as f:
            f.write('\n'.join(_db_lines))
        del _db_lines
        _database_cmd = NcbimakeblastdbCommandline(
            cmd=os.path.join(args.blast, 'makeblastdb'),
            dbtype='nucl',
            input_file=os.path.join(args.tmp, 'pre_blast_db.fasta'))
        _database_cmd()
        ## blast
        _blastn_cmd = NcbiblastnCommandline(
            cmd=os.path.join(args.blast, 'blastn'),
            query=os.path.join(args.tmp, 'pre_probe.fasta'),
            db=os.path.join(args.tmp, 'pre_blast_db.fasta'),
            task='blastn-short',
            dust='no',
            word_size=min(round(args.length/2)-1, 11),
            outfmt='\"6 qacc sacc length pident evalue\"',
            num_threads=args.threads,
            evalue=10,
            out=os.path.join(args.tmp, 'pre_blast.txt'),
            max_hsps=1,
            max_target_seqs=len(_row['assembly']))
        _blastn_cmd()
        ## parse result
        _blast_df = pd.read_table(os.path.join(args.tmp, 'pre_blast.txt'),
                                  names=['query', 'subject', 'length', 'identity', 'evalue'])
        _blast_df = _blast_df[(_blast_df['length'] == args.length) & (_blast_df['identity'] == 100)]
        _check_list = _group_sample_table.loc[_group_sample_table['group'] == _row['group'], 'assembly'].to_list()
        _probe_list = []
        for _sample in _check_list:
            _probe_list.append(set(_blast_df.loc[~(_blast_df['subject'].isin(_sample))]['query'].to_list()))
        _probe_result = list(reduce(lambda x,y: x & y, _probe_list))
        del _probe_list
        _probe_list = [_ for _ in SeqIO.parse(os.path.join(args.tmp, 'pre_probe.fasta'), 'fasta') if _.id in
                       _probe_result]
        SeqIO.write(_probe_list, os.path.join(args.tmp, 'probe.fasta'),'fasta')
        del _probe_list
        ## BLAST
        ap = BLASTParse(os.path.join(args.tmp, 'probe.fasta'),
                        _row['assembly'],
                        args.database,
                        args.tmp,
                        args.length,
                        args.threads,
                        args.blast)
        ap.blast_cmd()
        _result_tb = ap.blast_parse()
        _result_tb['group'] = _row['group']
        _result_tb[['group', 'assembly', 'seq', 'position', 'GC']]. \
            to_csv(os.path.join(args.output, _row['group'] + '.txt'),
                   sep='\t',
                   index=False)
