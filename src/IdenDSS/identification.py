# -*- coding: utf-8 -*-
# @Time : 2020/12/4 13:01
# @Author : Zhongyi Hua
# @FileName: probe_utils.py
# @Usage: The core module for identifying DSS
# @Note:
# @E-mail: njbxhzy@hotmail.com

from functools import reduce
from collections import defaultdict
from typing import Dict
import logging

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Blast.Applications import NcbiblastnCommandline


def probe_generate(_seq_obj, _length=20):
    """
    Generate probes for DSS identification
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


class IdentifyOneGroup:
    """
    Identify DSS for one group
    """
    def __init__(self, data_info, set_info, _group: str, sam2asm_dict: Dict, db_dict: Dict):
        """
        : return data_info    : DataInfo    : Two attributes: length, database
            : db     : Path : Path to the BLAST database
            : length : int  : DSS length
        : return set_info     : SettingInfo :
            : tmp     : Path : Temporary directory Path
            : threads : int  : The threads number for parallel software.
            : blast   : Path : BLAST executor path. If it in the PATH, this could be empty
        : return group        : str  : Group name
        : return sam2asm_dict : Dict : sample and assembly dict.
            {'sample id 1': [assembly1, assembly2, ...], 'sample id 2': [assembly1, assembly2, ...]}
        : return db_dict      : Dict : database sequences dict.
            {'assembly_id1': SeqRecord obj 1, 'assembly_id2': SeqRecord obj 2, ...}
        : return _asm_lst     : list : cache for asm_lst
        """
        self.datainfo = data_info
        self.setinfo = set_info
        self.group = _group
        self.sam2asm_dict = sam2asm_dict
        self.db_dict = db_dict
        self._asm_lst = None

    @property
    def asm_lst(self):
        """
        Generate assmbly list for this group from sample2assembly dict
        """
        if not self._asm_lst:
            self._asm_lst = []
            for _asmlst in self.sam2asm_dict.values():
                self._asm_lst += _asmlst
        return self._asm_lst

    def convserved_kmer(self) -> None:
        """
        Generate k-mers converved in intraspecies assemblies.
        FUTURE PROSCPECT: Using KMC3 accelerate this process
        """
        fix_sample = list(self.sam2asm_dict.keys())[0]
        _sample_probe = defaultdict()
        for _sample, _assembly in self.sam2asm_dict.items():
            _probe = defaultdict()
            if _sample == fix_sample:
                for _asm in _assembly:
                    for _k, _v in probe_generate(self.db_dict[_asm], self.datainfo.length).items():
                        _probe[_k] = _v
                _sample_probe[_sample] = _probe
            else:
                for _asm in _assembly:
                    for _k, _v in probe_generate(self.db_dict[_asm], self.datainfo.length).items():
                        _probe[_k] = _v
                        _probe[str(Seq(_k).reverse_complement())] = _v
                _sample_probe[_sample] = _probe
        try:
            _probe_result = set(_sample_probe[fix_sample].keys()) & \
                reduce(lambda x, y: x & y,
                       [set(_sample_probe[_sample].keys()) for _sample in list(self.sam2asm_dict.keys())[1:]])
        except TypeError:
            _probe_result = set(_sample_probe[fix_sample].keys())
        _probe_result = {_probe: _sample_probe[fix_sample][_probe] for _probe in list(_probe_result)}
        _probe_lines = ['>' + seq_id + '\n' + seq for seq, seq_id in _probe_result.items()]
        if len(_probe_lines) == 0:
            print(f'No conserved k-mers in {self.group} samples')
        (self.setinfo.tmp / 'probe.fasta').write_text('\n'.join(_probe_lines))

    def blast_cmd(self) -> None:
        """
        Generate second run BLAST and run it.
        """
        blastn_cmd = NcbiblastnCommandline(
            cmd=str(self.setinfo.blastn),
            query=self.setinfo.tmp / 'probe.fasta',
            db=self.datainfo.db,
            task='blastn-short' if self.datainfo.length < 30 else 'blastn',
            dust='no',
            word_size=round(self.datainfo.length/2)-1,
            outfmt='\"6 qacc sacc length pident evalue\"',
            num_threads=self.setinfo.threads,
            evalue=10,
            out=self.setinfo.tmp / (self.group + '.txt'),
            max_hsps=1,
            max_target_seqs=len(self.asm_lst)+1)
        blastn_cmd()

    def blast_parse(self) -> pd.DataFrame:
        """
        parse sconde run result and return a DataFrame
        """
        _blast_df = pd.read_table(self.setinfo.tmp / (self.group + '.txt'),
                                  names=['query', 'subject', 'length', 'identity', 'evalue'])
        keep_set = set(_blast_df['query'].to_list())
        _blast_df = _blast_df[(_blast_df['length'].values == self.datainfo.length) & (_blast_df['identity'].values == 100)]
        # Find DSS **in** other species
        drop_set = set(_blast_df.loc[_blast_df['subject'].map(lambda x: x not in self.asm_lst)]['query'].to_list())
        _dss_list = list(keep_set - drop_set)
        # DSS set empty
        if not _dss_list:
            _result_tb = pd.DataFrame.from_dict(
                {0: {'assembly': '', 'seq':  '', 'position': '', 'GC': ''}},
                'index')
        else:
            _kmer_dict = {_.id: _.seq for _ in
                          SeqIO.parse(self.setinfo.tmp / 'probe.fasta', 'fasta')
                          }
            _tmp_list = []
            for _dss in _dss_list:
                _tmp_list.append(
                                 {'assembly': '_'.join(_dss.split('_')[:-1]),
                                  'seq': _kmer_dict[_dss],
                                  'position': _dss.split('_')[-1] +
                                              '-' +
                                              str(int(_dss.split('_')[-1]) + self.datainfo.length - 1),
                                  'GC': GC(_kmer_dict[_dss])}
                                 )
            _result_tb = pd.DataFrame(_tmp_list)
        return _result_tb


def iden_main(data_info, set_info) -> None:
    """
    The main interface for DSS identification
    """
    seq_dict = {_.id: _ for _ in SeqIO.parse(data_info.db, 'fasta')}
    _table = pd.read_table(data_info.input, names=['group', 'sample', 'assembly'], dtype=str)
    # by group for generate probe
    _group_table = _table.groupby(['group'])['assembly'].apply(list).reset_index(name='assembly')
    # by sample for pre-BLAST
    set_info.logger.setLevel(level=logging.INFO)
    _group_sample_table = _table.\
        groupby(['group', 'sample'])['assembly'].\
        apply(list).\
        reset_index(name='assembly')
    for _idx, _row in _group_table.iterrows():
        _group = _row['group']
        set_info.logger.info(f'Begin to identify {_group} putative DSS, please wait for a moment!')
        set_info.initiate()
        _tmp_list = _group_sample_table.\
            loc[_group_sample_table['group'].values == _group].\
            to_dict(orient='records')
        sam2asm_dict = {_sample['sample']: _sample['assembly'] for _sample in _tmp_list}
        # BLAST
        one_group = IdentifyOneGroup(data_info, set_info, _row['group'], sam2asm_dict, seq_dict)
        one_group.convserved_kmer()
        set_info.logger.info('Intraspeceis conserved k-mers, done.')
        one_group.blast_cmd()
        set_info.logger.info('Putative DSS identification, done')
        _result_tb = one_group.blast_parse()
        _result_tb['group'] = _group
        _result_tb[['group', 'assembly', 'seq', 'position', 'GC']]. \
            to_csv(data_info.output / (_group + '.txt'),
                   sep='\t',
                   index=False)
        set_info.logger.info('Parse results, done')
        set_info.autoclean()
    set_info.logger.info('All groups done!')
