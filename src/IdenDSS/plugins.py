# -*- coding: utf-8 -*-
# @Time : 2021/11/13 10:31
# @Author : Zhongyi Hua
# @FileName: plugins.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com

import re
from functools import reduce
from collections import defaultdict
from typing import List
from pathlib import Path
from random import sample

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Restriction import Restriction as Rst
from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict
from Bio.SeqUtils import GC


def get_flanking(_seq, _start, _end, _flanking_length):
    """
    Get DSS with flanking sequence from thw whole sequences. 1-based index.
    """
    if _start <= _flanking_length:
        return _seq.seq[len(_seq.seq)-_flanking_length-1+_start:] + _seq.seq[0: _end+_flanking_length]
    if _end >= len(_seq.seq)-_flanking_length:
        return _seq.seq[_start-_flanking_length-1:] + _seq.seq[:_flanking_length-(len(_seq.seq)-_end)]
    return _seq.seq[_start-_flanking_length-1: _end+_flanking_length]


def get_seq(_seq, _start, _end):
    """
    Get sequence from thw whole sequences. 1-based index.
    """
    return _seq.seq[_start-1: _end]



class IndexDb:
    """
    Index reference database for DSS identification
    """
    def __init__(self, data_info, set_info):
        """
        : param data_info : DataInfo    : containing: input, output, circular, length attributes.
            : input  : the input fasta path
            : output : the output prepared database path
        : param set_info  : SettingInfo : containing blast attribute.
            : blast  : BLAST executor path. If it in the PATH, this could be empty
        """
        self.datainfo = data_info
        self.setinfo = set_info

    def _database_generate(self) -> None:
        """
        If the sequences are circular (e.g. chloroplast), then copy a part of 
        the sequence from the head to the end.
        """
        fasta_in = SeqIO.parse(self.datainfo.input, 'fasta')
        fasta_out = []
        for assembly in fasta_in:
            assembly.seq = assembly.seq + assembly.seq[0: self.datainfo.length-1]
            fasta_out.append(assembly)
        SeqIO.write(fasta_out, self.datainfo.output, 'fasta')

    def _copy_file(self) -> None:
        """
        If the sequences are not circular, just copy the file.
        """
        with self.datainfo.output.open(mode='wb') as fid:
            fid.write(self.datainfo.input.read_bytes())

    def _index_blast(self) -> None:
        """
        make blast db
        """
        database_cmd = NcbimakeblastdbCommandline(
            cmd= self.setinfo.makeblastdb,
            dbtype='nucl',
            input_file=self.datainfo.output)
        database_cmd()
        self.setinfo.logger.info('IdenDSS database created success!')

    def index(self):
        """
        The main index interface
        """
        try:
            if self.datainfo.circular:
                self._database_generate()
            else:
                self._copy_file()
            self._index_blast()
        except Exception as e:
            self.setinfo.logger.error(e)

class RFLP:
    """
    Scan RFLP sites on a fixed region (*start* to *end*) on a batch of sequences.
    """
    def __init__(self, _seqs: List[SeqRecord], _start: int, _end: int, _enzymes: List[str]):
        """
        : param _seqs    : a list of Bio.SeqRecord.SeqRecord
        : param _start   : The start positon of the scanned region. 1-based.
        : param _end     : The end position of the scanned region. 1-based.
        : param _enzymes : a list of character contained Restriction Enzymes' names
        """
        self.seqs = _seqs
        self.start = _start
        self.end = _end
        self.enzymes = _enzymes

    @staticmethod
    def create_enzyme(name):
        """
        Create a restriction endonuclease enzyme object
        """
        e_types = [x for t, (x, y) in typedict.items() if name in y][0]
        enzyme_types = tuple(getattr(Rst, x) for x in e_types)
        return Rst.RestrictionType(name, enzyme_types, rest_dict[name])

    def get_site_interval(self, name):
        """
        Get bidirectional enzyme recognization position
        : param name : str : The restriction endonuclease enzyme name
        """
        enzyme = self.create_enzyme(name)
        rco_pos, *_, rco_site = enzyme.characteristic()
        return {name: np.array([rco_pos, len(rco_site) - rco_pos])}

    def analysis_pos(self, _seq, _rb, _rbi):
        """
        : param _seq : Bio.Seq.Seq : The seq need to be analysed
        : param _rb  : Bio.Restriction.RestrictionBatch : The object used to search sequence recognization site.
        : param _rbi : Dict : The dict used to correct the whole recognized site (Bio.Restriction just return the cut position)
        :return      : list : The names of enzymes that could digest this seq.
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
        """
        main function to perform the search
        """
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


def search_rflp(data_info, set_info, enzymes_list):
    """
    The main interface for searching RFLP sites on DSSs
    """
    # Restriction Enzymes
    # import DSS result
    with open(data_info.input, 'r', encoding='utf-8') as f_in:
        file_list = f_in.read().strip().split('\n')
    # import database
    _meta_fasta = SeqIO.to_dict(SeqIO.parse(data_info.db, 'fasta'))
    # do RFLP search
    for _file in file_list:
        try:
            group_name = Path(_file).stem
            set_info.logger.info(f'Begin to search {group_name} RFLP sites')
            _dss_tb = pd.read_table(_file)
            if sum(_dss_tb['seq'].isna()) == 1:
                continue
            _dss_tb[['start', 'end']] = _dss_tb.apply((lambda x: x['position'].split('-')), axis=1, result_type="expand")
            in_list = []
            for _idx, _item in _dss_tb.iterrows():
                if not data_info.circular:
                    if int(_item['start']) < 250 or int(_item['end']) > len(_meta_fasta[_item['assembly']].seq) - 250:
                        continue
                in_list.append(
                    SeqRecord(seq=get_flanking(_meta_fasta[_item['assembly']], int(_item['start']), int(_item['end']), 250),
                            id=_item['assembly']+'-'+_item['start']+'-'+_item['end'])
                )
            rflp_ins = RFLP(in_list, 251, 290, enzymes_list)
            _tmp_df = rflp_ins.rflp()
            if _tmp_df.empty:
                continue
            else:
                _tmp_df = _tmp_df.merge(_dss_tb, how='left')
                _tmp_df[['group', 'assembly', 'seq', 'position', 'GC', 'enzyme']].to_csv(
                    data_info.output/ (group_name + '_rflp.txt'),
                    sep='\t',
                    index=False
                )
            set_info.logger.info(f'{group_name} RFLP sites searches done')
        except Exception as e:
            set_info.logger.error(e)
    set_info.logger.info(f'All RFLP sites searches done')

def _combine_per_dss(_dss_tb):
    """Generate combined DSS from one DSS result file
    Args:
        _dss_tb (pd.DataFrame): The input result dataframe
    """
    _dss_tb[['start', 'end']] = _dss_tb.position.str.split('-', expand=True)
    _dss_tb[['start', 'end']] = _dss_tb[['start', 'end']].astype(int)
    _dss_tb.sort_values(by=["assembly", "start"], inplace=True)
    _combined_list = []
    _seq = ''
    _asm = _dss_tb.iloc[0, 1]
    _start_pos = None
    _end_pos = None
    pointer = -1
    for _idx, _row in _dss_tb.iterrows():
        if (_row['start'] != pointer + 1) or (_row['assembly'] != _asm):
            _combined_list.append({
                'assembly': _asm,
                'seq': _seq,
                'start': _start_pos,
                'end': _end_pos,
                'GC': GC(_seq)})
            _asm = _row['assembly']
            _seq = _row['seq']
            _start_pos = _row['start']
            _end_pos = _row['end']
            pointer = _row['start']
        else:
            _end_pos = _row['end']
            _seq += _row['seq'][-1]
            pointer += 1
    # add the last one!
    _combined_list.append({'assembly': _asm,'seq': _seq, 'start': _start_pos, 'end': _end_pos, 'GC': GC(_seq)})
    _combined_list = _combined_list[1:]
    combined_res = pd.DataFrame(_combined_list)
    combined_res['group'] = _dss_tb.iloc[0, 0]
    combined_res['position'] = combined_res.apply(lambda x: str(x['start']) + '-' + str(x['end']), axis=1)
    return combined_res[['group', 'assembly', 'seq', 'position', 'GC']]

    
def combine_dss(data_info, set_info):
    """
    Generate combined DSS from DSS results
    """
    file_lst = data_info.input.read_text().strip().split('\n')
    for _file in file_lst:
        try:
            _dss_tb = pd.read_table(Path(_file))
            if sum(_dss_tb['seq'].isna()) == 1:
                continue
            group_name = Path(_file).stem
            set_info.logger.info(f'Begin to combine {group_name} DSSs')
            _combined_dss_tb = _combine_per_dss(_dss_tb)
            _combined_dss_tb.to_csv(data_info.output / (group_name + '_combined.txt'), sep='\t', index=False)
            set_info.logger.info(f'Combining {group_name} DSSs done')
        except Exception as e:
            set_info.logger.error(e)
    set_info.logger.info('Combining all groups DSSs done')


def summary_dss(data_info, set_info):
    """
    Summary DSS number for a list of DSS results.
    """
    file_lst = data_info.input.read_text().strip().split('\n')
    tmp_num_lst = []
    set_info.logger.info(f'Begin to summary all groups results')
    for _file in file_lst:
        try:
            _dss_tb = pd.read_table(_file)
            group_name = _dss_tb.iloc[0,0]
            if sum(_dss_tb['seq'].isna()) == 1:
                dss_num = 0
            else:
                dss_num = len(_dss_tb)
            tmp_num_lst.append({'group': group_name, 'dss_num': dss_num})
        except Exception as e:
            set_info.logger.error(e)
    pd.DataFrame(tmp_num_lst).to_csv(data_info.output / 'summary.tsv', sep='\t', index=False)
    set_info.logger.info('Summary all groups results done')


def flanking(data_info, set_info):
    """Generate flanking sequence, but not design primer.
    """
    _meta_fasta = SeqIO.to_dict(SeqIO.parse(data_info.db, 'fasta'))
    file_lst = data_info.input.read_text().strip().split('\n')
    f_len = data_info.length
    for _file in file_lst:
        try:
            group_name = Path(_file).stem
            set_info.logger.info(f'Begin to generate {group_name} flanking sequences')
            _dss_tb = pd.read_table(_file)
            if sum(_dss_tb['seq'].isna()) == 1:
                continue
            _dss_tb[['start', 'end']] = _dss_tb.apply((lambda x: x['position'].split('-')), axis=1, result_type="expand")
            in_list = []
            for _idx, _item in _dss_tb.iterrows():
                _ref_asm = _meta_fasta[_item['assembly']]
                _start = int(_item['start']) - 1
                _end = int(_item['end']) + 1
                if not data_info.circular:
                    if _start < f_len or _end > len(_ref_asm) - f_len:
                        continue
                if _start <= f_len:
                    left_flanking = get_seq(_ref_asm, len(_ref_asm) - f_len - 1, len(_ref_asm)) + get_seq(_ref_asm, 1, _start)
                    right_flanking = get_seq(_ref_asm, _end, _end + f_len)
                if _end >= len(_ref_asm)- f_len:
                    left_flanking = get_seq(_ref_asm, _start - f_len, _start)
                    right_flanking = get_seq(_ref_asm, _end, len(_ref_asm)) + get_seq(_ref_asm,1, f_len-(len(_ref_asm)-_end))
                else:
                    left_flanking = get_seq(_ref_asm, _start - f_len, _start)
                    right_flanking = get_seq(_ref_asm, _end, _end + f_len)
                in_list.append({'left_flanking' : left_flanking,
                                'right_flanking': right_flanking,
                                'start'         : _item['start'],
                                'end'           : _item['end']})
            _tmp_df = pd.DataFrame(in_list)
            _tmp_df = _tmp_df.merge(_dss_tb, how='left')
            _tmp_df[['group', 'assembly', 'seq', 'left_flanking', 'right_flanking', 'position', 'GC']].\
            to_csv(
                data_info.output/ (group_name + '_flanking.txt'),
                sep='\t',
                index=False
            )
            set_info.logger.info(f'{group_name} flanking sequence extract done')
        except Exception as e:
            set_info.logger.error(e)


def convert(data_info, set_info):
    """Convert DSS on old reference to the new reference
    Note: Notice the DSS across the gap in circular sequence
    """
    _meta_fasta = SeqIO.to_dict(SeqIO.parse(data_info.db, 'fasta'))
    line_lst = data_info.input.read_text().strip().split('\n')
    for _line in line_lst:
        try:
            _file, _new_ref = _line.split('\t')
            group_name = Path(_file).stem
            set_info.logger.info(f'Begin to convert {group_name} reference sequences')
            _new_ref_str = str(_meta_fasta[_new_ref].seq)
            _dss_tb = pd.read_table(_file)
            group_name = _dss_tb.iloc[0,0]
            if sum(_dss_tb['seq'].isna()) == 1:
                continue
            out_lst = []
            for _idx, _item in _dss_tb.iterrows():
                _seq = _item['seq']
                _gc = _item['GC']
                _match_info = re.search(_item['seq'], _new_ref_str)
                if not _match_info:
                    _seq = str(Seq(_seq).reverse_complement())
                    _match_info = re.search(_seq, _new_ref_str )
                    _gc = '{:.1f}'.format(100 - float(_item['GC']))
                if not _match_info:
                    # some reverse complement DSS across the gap
                    # may cause this error. Just skip them.
                    continue
                _new_start, _new_end = _match_info.span()
                _new_start += 1 # for 1-based index
                out_lst.append({
                    'assembly': _new_ref,
                    'seq': _seq,
                    'position': str(_new_start) + '-' + str(_new_end),
                    'GC': _gc
                })
            converted_tb = pd.DataFrame(out_lst)
            converted_tb['group'] = group_name
            converted_tb[['group', 'assembly', 'seq', 'position', 'GC']].\
                    to_csv(data_info.output / (group_name + '_converted.txt'), sep='\t', index=False)
            set_info.logger.info(f'{group_name} has been converted to new reference sequence')
        except Exception as item_exception:
            set_info.logger.error(item_exception)

def generate_dict(combined_dss_tb: pd.DataFrame, length: int):
    """Generate a start position dict for each combined DSS.

    Args:
        combined_dss_tb (pd.DataFrame): The combined DSS table
        length (int): The DSS length

    Returns:
        pos_start_dict (dict[List]): {0: [55873,55874,...], 1:[66112,66113,...]}
    """
    cdss_start_dict = defaultdict(list)
    combined_dss_tb[['start', 'end']] = combined_dss_tb.position.str.split('-', expand=True)
    combined_dss_tb[['start', 'end']] = combined_dss_tb[['start', 'end']].astype(int)
    for _idx, _row in combined_dss_tb.iterrows():
        cdss_start_dict[_idx] = [_row['assembly'] + str(_) for _ in range(_row['start'], _row['end'] - length + 2)]
    return cdss_start_dict

def resample(data_info, set_info):
    """_summary_

    Args:
        data_info (_type_): _description_
        set_info (_type_): _description_
    """
    file_lst = data_info.input.read_text().strip().split('\n')
    sample_num = data_info.length
    for _file in file_lst:
        try:
            group_name = Path(_file).stem
            set_info.logger.info(f'Begin to sample {sample_num} records from {group_name}')
            _dss_tb = pd.read_table(_file)
            _dss_tb[['start', 'end']] = _dss_tb.position.str.split('-', expand=True)
            _dss_tb[['start', 'end']] = _dss_tb[['start', 'end']].astype(int)
            _dss_tb['tmpid'] = _dss_tb.apply(lambda x: x['assembly'] + str(x['start']), axis=1)
            dss_length = len(_dss_tb['seq'][0])
            if len(_dss_tb) < sample_num:
                sample_tb = _dss_tb
            else: 
                sample_start_lst = []
                _combine_dss_tb = _combine_per_dss(_dss_tb)
                pos_start_dict = generate_dict(_combine_dss_tb, dss_length)
                # give priority to different combined DSS
                while len(sample_start_lst) < sample_num:
                    sample_start_lst += [_.pop() for _ in pos_start_dict.values() if _]
                # reduce to the sample number
                sample_start_lst = sample(sample_start_lst, sample_num)
                sample_tb = _dss_tb[_dss_tb['tmpid'].isin(sample_start_lst)].copy()
            sample_tb.drop(['start', 'end', 'tmpid'], axis=1, inplace=True)
            sample_tb.to_csv(
                    data_info.output/ (group_name + '_sampled.txt'),
                    sep='\t',
                    index=False
                )
            set_info.logger.info(f'{group_name} has been sampled to {len(sample_tb)} records')
        except Exception as item_exception:
            set_info.logger.error(item_exception)
