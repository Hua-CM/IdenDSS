# -*- coding: utf-8 -*-
# @Time : 2021/6/26 17:48
# @Author : Zhongyi Hua
# @FileName: primer_utils.py
# @Usage:
# @Note:

# @E-mail: njbxhzy@hotmail.com

import subprocess
from collections import deque
from pathlib import Path

from Bio import SeqIO
import pandas as pd

from .plugins import get_seq


def pre_cfg(_meta_tb, _meta_fasta, _outpath, _circular):
    """
    Prepare Primer3 input file
    """
    _write_list = []
    for _idx, _item in _meta_tb.iterrows():
        if not _circular:
            if int(_item['start']) < 400 or int(_item['end']) > len(_meta_fasta[_item['assembly']].seq)-400:
                continue
        _sequence_id = "SEQUENCE_ID=" + \
                       str(_meta_fasta[_item['assembly']].id) + \
                       "_" + str(_item['start']) + \
                       "_" + str(_item['end'])
        _sequence_template = "SEQUENCE_TEMPLATE=" + str(
            get_seq(_meta_fasta[_item['assembly']], int(_item['start']), int(_item['end']), 400)
        )
        _write_list.append(_sequence_id)
        _write_list.append(_sequence_template)
        _write_list.append('SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,380,' +
                           str(381 + int(_item['end']) - int(_item['start'])) +
                           ',400')
        _write_list.append('=')
        _outpath.write_text('=\n')
        _outpath.write_text('\n'.join(_write_list))
        _outpath.write_text("\n")


def parse_output(_p3out, _output) -> None:
    """
    Parse Primer3 result to readable tabular text
    : param  _p3out : Path : Primer3 result path
    """
    _file = deque(_p3out.read_text().strip().split('\n'))
    _module = []
    _record_list = []
    for _line in _file:
        _module.append(_line)
        if _line == '=':
            _model_dict = dict(_.split('=') for _ in _module)
            if 'SEQUENCE_ID' in _model_dict:
                if not _model_dict['PRIMER_PAIR_NUM_RETURNED'] == '0':
                    _record = {'Assembly': '_'.join(_model_dict['SEQUENCE_ID'].split('_')[:-2]),
                               'GenomePosition': '-'.join(_model_dict['SEQUENCE_ID'].split('_')[-2:]),
                               'Product Size': _model_dict['PRIMER_PAIR_0_PRODUCT_SIZE'],
                               'Forward':  _model_dict['PRIMER_LEFT_0_SEQUENCE'],
                               'Forward_position': _model_dict['PRIMER_LEFT_0'].split(',')[0],
                               'Forward_length': _model_dict['PRIMER_LEFT_0'].split(',')[1],
                               'Forward_TM': _model_dict['PRIMER_LEFT_0_TM'],
                               'Reverse': _model_dict['PRIMER_RIGHT_0_SEQUENCE'],
                               'Reverse_position': _model_dict['PRIMER_RIGHT_0'].split(',')[0],
                               'Reverse_length': _model_dict['PRIMER_RIGHT_0'].split(',')[1],
                               'Reverse_TM': _model_dict['PRIMER_RIGHT_0_TM'],
                               'Seq': _model_dict['SEQUENCE_TEMPLATE']
                               }
                    _record_list.append(_record)
                else:
                    _record = {'Assembly': '_'.join(_model_dict['SEQUENCE_ID'].split('_')[:-2]),
                               'GenomePosition': '-'.join(_model_dict['SEQUENCE_ID'].split('_')[-2:]),
                               'Seq': _model_dict['SEQUENCE_TEMPLATE']
                               }
                    _record_list.append(_record)
                _module = []
            else:
                continue
    pd.DataFrame(_record_list).to_csv(_output, sep='\t', index=False)


def primer_main(data_info, set_info):
    """
    The main interface for Primer design.
    """
    file_list = data_info.input.read_text().strip().split('\n')
    _meta_fasta = SeqIO.to_dict(SeqIO.parse(data_info.database, 'fasta'))
    for _file in file_list:
        _meta_tb = pd.read_table(_file)
        _prefix = Path(_file).stem
        if sum(_meta_tb['seq'].isna()) == 1:
            continue
        _meta_tb[['start', 'end']] = _meta_tb.apply((lambda x: x['position'].split('-')), axis=1, result_type="expand")
        set_info.tmp.initiate()
        pre_cfg(_meta_tb[['assembly', 'start', 'end']],
                _meta_fasta,
                set_info.tmp/ (_prefix + '.p3in'),
                data_info.circular)
        
        setting_file = (Path(__file__).parent / 'template' / 'DSS_settings.txt').resolve()
        subprocess.run(
            ['primer3_core',
             '--p3_settings_file=' + setting_file,
             '--output=' + (set_info.tmp / (_prefix + '.p3out')),
             (set_info.tmp / _prefix + '.p3in')],
             check=True
        )
        parse_output(set_info.tmp / (_prefix + '.p3out'),
                     data_info.output / (_prefix + '_primer.txt'))
        set_info.tmp.autoclean()
