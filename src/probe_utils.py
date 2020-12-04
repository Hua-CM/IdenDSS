# -*- coding: utf-8 -*-
# @Time : 2020/12/4 13:01
# @Author : Zhongyi Hua
# @FileName: probe_utils.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import tempfile
import pandas as pd
from Bio import SeqIO

def parse_meta_info(_table_path):
    _table = pd.read_table(_table, names=['group', 'sample', 'assembly'], ,dtype=str)
    _table.groupby('group')['assembly'].apply(list).reset_index(name='assembly')

class GenerateProbe:
    def __init__(self, _table, _seq):
    """
    :param _table: meta info path, a str object
    :param _seq: database seq path, a str object
    """
    self.table = pd.read_table(_table, names=['group', 'sample', 'assembly'], ,dtype=str)
    self.seq = SeqIO.parse(,'fasta')
    self.assembly_dict = {}

    def handle_sample(self):
        

    def probe_generate(self, _length=20):
        sequence = self.seq.seq + self.seq.seq[0:_length]
        self.seq.id = self.seq.id.split('.')[0]
        _tmp_list = []
        for index in range(len(sequence) - (_length-1)):
            _tmp_list.append((str(sequence[index: index + _length]),
                              self.seq.id + '_' + str(index + 1)))
        self.assembly_dict = {_probe: _index for _probe, _index in _tmp_list}
        return self.assembly_dict

    def probe_save(self, save_path):
        probe_lines = ['>' + seq_id + '\n' + seq for seq, seq_id in self.assembly_dict.items()]
        with open(save_path, 'w') as f:
            f.write('\n'.join(probe_lines))
    
    def pre_blast():
        