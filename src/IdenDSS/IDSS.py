# -*- coding: utf-8 -*-
# @Time : 2020/11/1 19:18
# @Author : Zhongyi Hua
# @FileName: IDSS.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import argparse
import logging
from pathlib import Path
from tempfile import gettempdir
from sys import platform

from .identification import iden_main
from .primer_utils import primer_main
from .plugins import IndexDb, search_rflp, combine_dss, summary_dss, flanking, convert, resample
from .utils import enzymes_list, SettingInfo, DataInfo
from .validation import v_main


class CustomFormatter(argparse.HelpFormatter):
    """Reduce the redundant metavar
    """
    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

def get_args():
    """
    parse commend line arguments
    """
    parser = argparse.ArgumentParser(
        prog='IdenDSS',
        description='This script was for identifying DNA signature sequences(DSS)',
        formatter_class=CustomFormatter)
    sub_parser = parser.add_subparsers(title='Available', dest='database/iden/plugin/validate')
    sub_parser.required = True

    database_parser = sub_parser.add_parser(
        'index', help='Generate database for DSS identification', formatter_class=CustomFormatter)
    database_parser.add_argument('-i', '--input', type=Path, required=True, metavar='<FILE>',
                                 help='The genome fasta file')
    database_parser.add_argument('-l', '--length', type=int, default=40, metavar='<INT>',
                                 help='DSS length <default=40>')
    database_parser.add_argument('-c', '--circular',  action='store_true',
                                 help='the sequences are circular or not')
    database_parser.add_argument('-o', '--output', type=Path, required=True, metavar='<FILE>',
                                 help='The genome fasta output file (a BLAST database would '
                                      'be constructed)')
    database_parser.add_argument('--blast', type=Path, default=Path(''), dest='bin_dir', metavar='<DIR>',
                                 help='BLAST exec directory <If your BLAST software not in PATH>')
    database_parser.set_defaults(subcmd='index')

    probe_parser = sub_parser.add_parser(
        'identify', help='Identification DSS based on database', formatter_class=CustomFormatter)
    probe_parser.add_argument('-m', '--meta', type=Path, required=True, metavar='<FILE>',
                              help='The meta file')
    probe_parser.add_argument('-d', '--database', type=Path, required=True, metavar='<FILE>',
                              help='Database fasta <Corresponding BLAST database file must in the same '
                                   'directory>')
    probe_parser.add_argument('-l', '--length', type=int, default=40, metavar='<INT>',
                              help='DSS length <default=40>')
    probe_parser.add_argument('-c', '--circular', action='store_true',
                              help='The sequences are circular')
    probe_parser.add_argument('-o', '--output', type=Path, required=True, metavar='<DIR>',
                              help='The result directory')
    probe_parser.add_argument('-@', '--threads', type=int, default=4, metavar='<INT>',
                              help='Threads used in BLAST <Default 4>')
    probe_parser.add_argument('-t', '--tmp', type=Path, default=Path(gettempdir()), metavar='<DIR>',
                              help='Temporary directory path (Default system temporary directory)')
    probe_parser.add_argument('--blast', type=Path, default=Path(''), dest='bin_dir', metavar='<DIR>',
                              help='BLAST exec directory <If your BLAST software not in PATH>')
    probe_parser.set_defaults(subcmd='identify')

    plugin_parser = sub_parser.add_parser(
        'plugin', help='Some plugin for DSS results', formatter_class=CustomFormatter)
    plugin_parser.add_argument('-i', '--input', type=Path, required=True, metavar='<FILE>',
                               help='A meta file. One DSS result file path per line. \
                                     Some plugins may need an extra column.')
    plugin_parser.add_argument('-d', '--database', type=Path, required=True, metavar='<FILE>',
                               help='Database fasta (The database used to identify DSS)')
    plugin_parser.add_argument('-c', '--circular', action='store_true',
                               help='The sequences are circular')
    plugin_parser.add_argument('-o', '--output', type=Path, required=True, metavar='<DIR>',
                               help='The result directory')
    plugin_parser.add_argument('-t', '--tmp', type=Path, default=Path(gettempdir()), metavar='<DIR>',
                               help='Temporary directory path <Default system temporary directory>')
    plugin_parser.add_argument('--primer', action='store_true',
                               help='Design Primer (Need primer3_core in your PATH)')
    plugin_parser.add_argument('--rflp', action='store_true',
                               help='Identify restriction enzyme sits on DSS for putative RFLP method  (combined result is OK)')
    plugin_parser.add_argument('--combine', action='store_true',
                               help='Generate the corresponding combined DSS file')
    plugin_parser.add_argument('--statistic', action='store_true',
                               help='Count the DSS number')
    plugin_parser.add_argument('--flank', type=int, default=0, metavar='<INT>',
                               help='Generate the flanking sequence for each. e.g., --flank 200')
    plugin_parser.add_argument('--convert', action='store_true',
                               help='Convert the DSS to the another reference assembly \
                                     (Need a scenond column with the new reference name in meta file).')
    plugin_parser.add_argument('--sample', type=int, default=0, metavar='<INT>',
                                help='Sample N records from raw result file. This plugin give priority to \
                                      sample DSS from different combined DSS.')
    plugin_parser.add_argument('--bin', type=Path, default=Path(''), dest='bin_dir', metavar='<DIR>',
                               help='Primer3 exec directory <If your Primer3 software not in PATH>')
    plugin_parser.set_defaults(subcmd='plugin')

    validate_parser = sub_parser.add_parser(
        'validate', help='Validate DSS using HTS file. MUST BE DSS, NOT COMBINED DSS', formatter_class=CustomFormatter)
    validate_parser.add_argument('-i', '--input', type=Path, required=True, metavar='<FILE>',
                                 help='The DSS result file path')
    validate_parser.add_argument('--sp', required=True, type=str, dest='sp', metavar='<FILE>',
                                 help='The intraspecies HTS FASTQ file \
                                 (could be a list, seperate by comma)')
    validate_parser.add_argument('--bg', required=True, type=str, dest='bg', metavar='<FILE>',
                                 help='The  background species HTS \
                                 FASTQ file (could be a list, seperate by comma)')
    validate_parser.add_argument('-o', '--output', type=Path, required=True, metavar='<FILE>',
                                 help='The validated DSS result path')
    validate_parser.add_argument('-t', '--tmp', type=Path, default=Path(gettempdir()), metavar='<DIR>',
                                 help='Temporary directory path \
                                       <Default system temporary directory>')
    validate_parser.add_argument('--bin', type=Path, default=Path(''), metavar='<DIR>',
                                 help='KMC3 exec directory, containing kmc and kmc_dump \
                                       <If your KMC3 software not in PATH>')
    validate_parser.set_defaults(subcmd='validate')

    args = parser.parse_args()
    return args


def main():
    """
    The main interface
    """
    args = get_args()

    logger = logging.getLogger('IdenDSS')
    logger.setLevel(level=logging.INFO)
    logger.addHandler(logging.StreamHandler())

    if args.subcmd == 'index':
        set_info = SettingInfo(
            bin_dir=args.bin_dir,
            logger=logger)
        data_info = DataInfo(in_put=args.input,
                             output=args.output,
                             circular=args.circular,
                             length=args.length)
        index_db = IndexDb(data_info, set_info)
        index_db.index()

    if args.subcmd == 'identify':
        set_info = SettingInfo(tmp=args.tmp,
                               bin_dir=args.bin_dir,
                               threads=args.threads,
                               logger=logger)
        data_info = DataInfo(in_put=args.meta,
                             output=args.output,
                             db=args.database,
                             length=args.length,
                             circular=args.circular)
        iden_main(data_info, set_info)

    if args.subcmd == 'plugin':
        set_info = SettingInfo(tmp=args.tmp,
                               logger=logger)
        if platform == 'win32':
            set_info.clean_file(args.input)
        data_info = DataInfo(in_put=args.input,
                             output=args.output,
                             db=args.database,
                             circular=args.circular)
        if args.primer:
            primer_main(data_info, set_info)
        if args.rflp:
            search_rflp(data_info, set_info, enzymes_list)
        if args.combine:
            combine_dss(data_info, set_info)
        if args.statistic:
            summary_dss(data_info, set_info)
        if args.flank:
            data_info.length = args.flank
            flanking(data_info, set_info)
        if args.convert:
            convert(data_info, set_info)
        if args.sample:
            data_info.length = args.sample
            resample(data_info, set_info)

    if args.subcmd == 'validate':
        set_info = SettingInfo(tmp=args.tmp,
                               bin_dir=args.bin,
                               logger=logger)
        data_info = DataInfo(in_put=args.input,
                             output=args.output,
                             db=[Path(_) for _ in args.bg.split(',')],
                             db2=[Path(_) for _ in args.sp.split(',')])
        v_main(data_info, set_info)

if __name__ == '__main__':
    main()
