# -*- coding: utf-8 -*-
# @Time : 2021/11/13 10:31
# @Author : Zhongyi Hua
# @FileName: plugins.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Restriction import Restriction as Rst
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict
from functools import reduce
from collections import defaultdict


def get_seq(_seq, _start, _end, _flanking_length):
    if _start <= _flanking_length:
        return _seq.seq[len(_seq.seq)-_flanking_length-1+_start:] + _seq.seq[0: _end+_flanking_length]
    elif _end >= len(_seq.seq)-_flanking_length:
        return _seq.seq[_start-_flanking_length-1:] + _seq.seq[:_flanking_length-(len(_seq.seq)-_end)]
    else:
        return _seq.seq[_start-_flanking_length-1: _end+_flanking_length]


class Database:
    def __init__(self, fa_path, db_path, _exec):
        """
        Since chloroplast was a circular, we need cut the base at the beginning which equal to the length of the DSS
        to the end
        :param fa_path: the input fasta path
        :param db_path: the output prepared database path
        :param exec: the BLAST executor path. If it in the PATH, this could be empty
        """
        self.in_path = fa_path
        self.out_path = db_path
        self.exec = _exec

    def database_generate(self, _length=20):
        fasta_in = SeqIO.parse(self.in_path, 'fasta')
        fasta_out = []
        for assembly in fasta_in:
            assembly.seq = assembly.seq + assembly.seq[0:_length-1]
            fasta_out.append(assembly)
        SeqIO.write(fasta_out, self.out_path, 'fasta')

    def copy_file(self):
        open(self.out_path, 'wb').write(open(self.in_path, 'rb').read())

    def database_blast(self):
        database_cmd = NcbimakeblastdbCommandline(
            cmd=os.path.join(self.exec, 'makeblastdb'),
            dbtype='nucl',
            input_file=self.out_path)
        database_cmd()
        print('IdenDSS database created success!')


class RFLP:
    def __init__(self, _seqs, _start, _end, _enzymes):
        """
        :param _seqs: a list of Bio.SeqRecord.SeqRecord
        :param _start:
        :param _end:
        :param _enzymes: a list of character contained Restriction Enzymes' names
        """
        self.enzymes = _enzymes
        self.seqs = _seqs
        self.start = _start
        self.end = _end

    @staticmethod
    def create_enzyme(name):
        e_types = [x for t, (x, y) in typedict.items() if name in y][0]
        enzyme_types = tuple(getattr(Rst, x) for x in e_types)

        return Rst.RestrictionType(name, enzyme_types, rest_dict[name])

    def get_site_interval(self, name):
        enzyme = self.create_enzyme(name)
        rco_pos, *_, rco_site = enzyme.characteristic()
        return {name: np.array([rco_pos, len(rco_site) - rco_pos])}

    def analysis_pos(self, _seq, _rb, _rbi):
        """
        :param _seq: Bio.Seq.Seq
        :param _rb: Bio.Restriction.RestrictionBatch
        :param _rbi: Dict used to correct the whole recognized site (Bio.Restriction just return the cut position)
        :return:  list: enzymes' names
        """
        rb_res = _rb.search(_seq)
        res = {}
        res_enzymes = []
        for _enzyme, _pos_list in rb_res.items():
            _enzyme_list = []
            if len(_pos_list) == 1:
                _enzyme_list.append(_pos_list[0] + _rbi[str(_enzyme)])
            res[_enzyme] = _enzyme_list
        for _enzyme, _pos in res.items():
            if len(_pos) > 0:
                if _pos[0][1] > self.start and self.end > _pos[0][0]:
                    res_enzymes.append(str(_enzyme))
        return res_enzymes

    def rflp(self):
        enzyme_list = self.enzymes
        rb = reduce(lambda x, y: x + y, map(self.create_enzyme, enzyme_list))
        rbi = defaultdict()
        for _ in map(self.get_site_interval, enzyme_list):
            rbi.update(_)
        res_list = []
        for _seq in self.seqs:
            enzymes = self.analysis_pos(_seq.seq, rb, rbi)
            if enzymes:
                _tmp_dict = dict(zip(('assembly', 'start', 'end'), _seq.id.split('-')))
                _tmp_dict['enzyme'] = ','.join(enzymes)
                res_list.append(_tmp_dict)
            else:
                continue
        return pd.DataFrame(res_list)


def search_rflp(args):
    # Restriction Enzymes
    with open(os.path.abspath(os.path.join(__file__, '../../template/RestrictionEnzymes'))) as _f_in:
        enzymes_list = _f_in.read().strip().split('\n')
    # import DSS result
    with open(args.input, 'r') as f_in:
        file_list = f_in.read().strip().split('\n')
    # import database
    _meta_fasta = SeqIO.to_dict(SeqIO.parse(args.database, 'fasta'))
    # do RFLP search
    for _file in file_list:
        _dss_tb = pd.read_table(_file)
        if sum(_dss_tb['seq'].isna()) == 1:
            continue
        _dss_tb[['start', 'end']] = _dss_tb.apply((lambda x: x['position'].split('-')), axis=1, result_type="expand")
        in_list = []
        _write_list = []
        for _idx, _item in _dss_tb.iterrows():
            if not args.circular:
                if int(_item['start']) < 250 or int(_item['end']) > len(_meta_fasta[_item['assembly']].seq) - 250:
                    continue
            in_list.append(
                SeqRecord(seq=get_seq(_meta_fasta[_item['assembly']], int(_item['start']), int(_item['end']), 250),
                          id=_item['assembly']+'-'+_item['start']+'-'+_item['end'])
            )
        rflp_ins = RFLP(in_list, 251, 290, enzymes_list)
        _tmp_df = rflp_ins.rflp()
        if _tmp_df.empty:
            continue
        else:
            _tmp_df = _tmp_df.merge(_dss_tb, how='left')
            _tmp_df[['group', 'assembly', 'seq', 'position', 'GC', 'enzyme']].to_csv(
                os.path.join(args.output, os.path.splitext(os.path.split(_file)[1])[0] + '_rflp.txt'),
                sep='\t',
                index=False
            )
