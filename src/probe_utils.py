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
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from functools import reduce
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastnCommandline


def del_dir(_dir):
    for r, d, f in os.walk(_dir):
        for files in f:
            os.remove(os.path.join(r, files))
        os.removedirs(r)


def probe_generate(_seq_obj, _length=20):
    """

    :param _seq_obj: a Bio.Seq object
    :param _length:
    :return:
    """
    sequence = _seq_obj.seq + _seq_obj.seq[0:_length]
    _seq_obj.id = _seq_obj.id
    assembly_dict = {}
    for index in range(len(sequence) - (_length - 1)):
        assembly_dict.setdefault(str(sequence[index: index + _length]), _seq_obj.id + '_' + str(index + 1))
    return assembly_dict


class BLASTParse:
    def __init__(self, _group, _seq_path, _asm_list, _database_path, _tmp_dir_path, _length, _thread, _exec=''):
        self.group = _group
        self.seq_path = _seq_path
        self.asm_list = _asm_list
        self.db_path = _database_path
        self.tmp_dir = _tmp_dir_path
        self.length = _length
        self.thread = _thread
        self.exec = _exec

    def blast_cmd(self):
        blastn_cmd = NcbiblastnCommandline(
            cmd=os.path.join(self.exec, 'blastn'),
            query=self.seq_path,
            db=self.db_path,
            task='blastn-short' if self.length < 30 else 'blastn',
            dust='no',
            word_size=round(self.length/2)-1,
            outfmt='\"6 qacc sacc length pident evalue\"',
            num_threads=self.thread,
            evalue=10,
            out=os.path.join(self.tmp_dir, self.group + '.txt'),
            max_hsps=1,
            max_target_seqs=len(self.asm_list)+1)
        blastn_cmd()

    def blast_parse(self):
        _blast_df = pd.read_table(os.path.join(self.tmp_dir, self.group + '.txt'),
                                  names=['query', 'subject', 'length', 'identity', 'evalue'])
        keep_set = set(_blast_df['query'].to_list())
        _blast_df = _blast_df[(_blast_df['length'].values == self.length) & (_blast_df['identity'].values == 100)]
        # Find DSS **in** other species
        drop_set = set(_blast_df.loc[_blast_df['subject'].map(lambda x: x not in self.asm_list)]['query'].to_list())
        _dss_list = list(keep_set - drop_set)
        # DSS set empty
        if not _dss_list:
            _result_tb = pd.DataFrame.from_dict(
                {0: {'assembly': '', 'seq':  '', 'position': '', 'GC': ''}},
                'index')
        else:
            _kmer_dict = {_.id: _.seq for _ in
                          SeqIO.parse(self.seq_path, 'fasta')
                          }
            _tmp_list = []
            for _dss in _dss_list:
                _tmp_list.append(
                                 {'assembly': '_'.join(_dss.split('_')[:-1]),
                                  'seq': _kmer_dict[_dss],
                                  'position': _dss.split('_')[-1] +
                                              '-' +
                                              str(int(_dss.split('_')[-1])+self.length-1),
                                  'GC': GC(_kmer_dict[_dss])}
                                 )
            _result_tb = pd.DataFrame(_tmp_list)
        return _result_tb


def iden_main(args):
    seq_dict = {_.id: _ for _ in SeqIO.parse(args.database, 'fasta')}
    _table = pd.read_table(args.meta, names=['group', 'sample', 'assembly'], dtype=str)
    # by group for generate probe
    _group_table = _table.groupby(['group'])['assembly'].apply(list).reset_index(name='assembly')
    # by sample for pre-BLAST
    _group_sample_table = _table.groupby(['group', 'sample'])['assembly'].apply(list).reset_index(name='assembly')
    for _idx, _row in _group_table.iterrows():
        os.mkdir(args.tmp)
        # pre-probe
        _tmp_list = _group_sample_table.\
            loc[_group_sample_table['group'].values == _row['group']].\
            to_dict(orient='records')
        sample_assembly = {_sample['sample']: _sample['assembly'] for _sample in _tmp_list}
        del _tmp_list
        fix_sample = list(sample_assembly.keys())[0]
        _sample_probe = defaultdict()
        for _sample, _assembly in sample_assembly.items():
            _probe = defaultdict()
            if _sample == fix_sample:
                for _asm in _assembly:
                    for _k, _v in probe_generate(seq_dict[_asm], args.length).items():
                        _probe[_k] = _v
                _sample_probe[_sample] = _probe
            else:
                for _asm in _assembly:
                    for _k, _v in probe_generate(seq_dict[_asm], args.length).items():
                        _probe[_k] = _v
                        _probe[Seq(_k).reverse_complement().__str__()] = _v
                _sample_probe[_sample] = _probe
        try:
            _probe_result = set(_sample_probe[fix_sample].keys()) & \
                reduce(lambda x, y: x & y,
                       [set(_sample_probe[_sample].keys()) for _sample in list(sample_assembly.keys())[1:]])
        except TypeError:
            _probe_result = set(_sample_probe[fix_sample].keys())
        _probe_result = {_probe: _sample_probe[fix_sample][_probe] for _probe in list(_probe_result)}
        _probe_lines = ['>' + seq_id + '\n' + seq for seq, seq_id in _probe_result.items()]
        with open(os.path.join(args.tmp, 'probe.fasta'), 'w') as f:
            f.write('\n'.join(_probe_lines))
        del _probe_result, _probe_lines
        # BLAST
        ap = BLASTParse(_row['group'],
                        os.path.join(args.tmp, 'probe.fasta'),
                        _row['assembly'],
                        args.database,
                        args.tmp,
                        args.length,
                        args.threads,
                        args.blast)
        ap.blast_cmd()
        _result_tb = ap.blast_parse()
        _result_tb['group'] = _row['group']
        _result_tb['assembly'] = _row['assembly'][0]
        _result_tb[['group', 'assembly', 'seq', 'position', 'GC']]. \
            to_csv(os.path.join(args.output, _row['group'] + '.txt'),
                   sep='\t',
                   index=False)
        del_dir(args.tmp)
