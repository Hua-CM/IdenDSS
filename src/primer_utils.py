# -*- coding: utf-8 -*-
# @Time : 2021/6/26 17:48
# @Author : Zhongyi Hua
# @FileName: primer_utils.py
# @Usage:
# @Note:

# @E-mail: njbxhzy@hotmail.com


import pandas as pd
import os
from Bio import SeqIO
from collections import deque
from probe_utils import del_dir
from plugins import get_seq


def pre_cfg(_meta_tb, _meta_fasta, _outpath, _circular):
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
    with open(_outpath, 'w') as f2:
        f2.write('=\n')
        f2.write('\n'.join(_write_list))
        f2.write("\n")


def parse_output(_p3out, _output):
    _a = open(_p3out, 'r')
    _file = deque(_a.read().splitlines())
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


def primer_main(args):
    with open(args.input) as _f_in:
        file_list = _f_in.read().strip().split('\n')
    _meta_fasta = SeqIO.to_dict(SeqIO.parse(args.database, 'fasta'))
    for _file in file_list:
        _meta_tb = pd.read_table(_file)
        if sum(_meta_tb['seq'].isna()) == 1:
            continue
        _meta_tb[['start', 'end']] = _meta_tb.apply((lambda x: x['position'].split('-')), axis=1, result_type="expand")
        os.mkdir(args.tmp)
        pre_cfg(_meta_tb[['assembly', 'start', 'end']],
                _meta_fasta,
                os.path.join(args.tmp, _file + '.p3in'),
                args.circular)
        setting_file = os.path.abspath(os.path.join(__file__, '../../template/DSS_settings.txt'))
        os.system(' '.join(['primer3_core',
                            '--p3_settings_file=' + setting_file,
                            '--output=' + os.path.join(args.tmp, _file + '.p3out'),
                            os.path.join(args.tmp, _file + '.p3in')]
                           )
                  )
        parse_output(os.path.join(args.tmp, _file + '.p3out'),
                     os.path.join(args.output, _file + '_primer.txt'))
        del_dir(args.tmp)
