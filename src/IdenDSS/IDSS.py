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

from .identification import iden_main
from .primer_utils import primer_main
from .plugins import IndexDb, search_rflp, combine_dss, summary_dss
from .utils import enzymes_list, SettingInfo, DataInfo
from .validation import v_main

def get_args():
    """
    parse commend line arguments
    """
    parser = argparse.ArgumentParser(
        prog='IdenDSS', description='This script was for identifying DNA signature sequences(DSS)')
    sub_parser = parser.add_subparsers(title='Available', dest='database/iden/plugin/validate')
    sub_parser.required = True

    database_parser = sub_parser.add_parser(
        'index', help='Generate database for DSS identification')
    database_parser.add_argument('-i', '--input', type=Path, required=True,
                                 help='<File path>  The genome fasta file')
    database_parser.add_argument('-l', '--length', type=int, default=40,
                                 help='<Int> DSS length <default=40>')
    database_parser.add_argument('-c', '--circular',  action='store_true',
                                 help='the sequences are circular or not')
    database_parser.add_argument('-o', '--output', type=Path, required=True,
                                 help='<File path>  The genome fasta output file (a BLAST database would '
                                      'be constructed)')
    database_parser.add_argument('--blast', type=Path, default=Path(''), dest='bin_dir',
                                 help='<Directory path> BLAST exec directory <If your BLAST software not in PATH>')                                  
    database_parser.set_defaults(subcmd='index')

    probe_parser = sub_parser.add_parser(
        'identify', help='Identification DSS based on database')
    probe_parser.add_argument('-m', '--meta', type=Path, required=True,
                              help='<File path> The meta file')
    probe_parser.add_argument('-d', '--database', type=Path, required=True,
                              help='<File path> Database fasta <Corresponding BLAST database file must in the same '
                                   'directory>')
    probe_parser.add_argument('-l', '--length', type=int, default=40,
                              help='<Int> DSS length <default=40>')
    probe_parser.add_argument('-c', '--circular', action='store_true',
                              help='the sequences are circular or not')
    probe_parser.add_argument('-o', '--output', type=Path, required=True,
                              help='<Firectory path> result directory')
    probe_parser.add_argument('-@', '--threads', type=int, default=4,
                              help='<Int> Threads used in BLAST <Default 4>')
    probe_parser.add_argument('-t', '--tmp', type=Path, default=Path(gettempdir()),
                              help='<Directory path> Temporary directory path \
                                    (Default system temporary directory)')
    probe_parser.add_argument('--blast', type=Path, default=Path(''), dest='bin_dir',
                              help='<Directory path> BLAST exec directory <If your BLAST software not in PATH>')
    probe_parser.set_defaults(subcmd='identify')

    plugin_parser = sub_parser.add_parser(
        'plugin', help='Some plugin for DSS results')
    plugin_parser.add_argument('-i', '--input', type=Path, required=True,
                               help='<File path>  A meta file. One DSS result file path per line (combined result is OK)')
    plugin_parser.add_argument('-d', '--database', type=Path, required=True,
                               help='<File path> Database fasta (The database used to identify DSS)')
    plugin_parser.add_argument('-c', '--circular', action='store_true',
                               help='the sequences are circular or not')
    plugin_parser.add_argument('-o', '--output', type=Path, required=True,
                               help='<Directory path> result directory')
    plugin_parser.add_argument('-t', '--tmp', type=Path, default=Path(gettempdir()),
                               help='<Dir path> Temporary directory path <Default system temporary directory>')
    plugin_parser.add_argument('--primer', action='store_true',
                               help='Design Primer (Need primer3_core in your PATH)')
    plugin_parser.add_argument('--rflp', action='store_true',
                               help='Identify restriction enzyme sits on DSS for putative RFLP method')
    plugin_parser.add_argument('--combine', action='store_true',
                               help='Generate the corresponding combined DSS file')
    plugin_parser.add_argument('--statistic', action='store_true',
                               help='Count the DSS number')
    plugin_parser.add_argument('--bin', type=Path, default=Path(''), dest='bin_dir',
                               help='<Directory path> Primer3 exec directory <If your Primer3 software not in PATH>')
    plugin_parser.set_defaults(subcmd='plugin')

    validate_parser = sub_parser.add_parser(
        'validate', help='Validate DSS using HTS file. MUST BE DSS, NOT COMBINED DSS')
    validate_parser.add_argument('-i', '--input', type=Path, required=True,
                                 help='<File path> ')
    validate_parser.add_argument('--sp', required=True, type=str, dest='sp',
                                 help='<File path> The intraspecies HTS FASTQ file \
                                 (could be a list, seperate by comma)')
    validate_parser.add_argument('--bg', required=True, type=str, dest='bg',
                                 help='<File path> The  background species HTS \
                                 FASTQ file (could be a list, seperate by comma)')
    validate_parser.add_argument('--min', type=int, default=2,
                                 help='<Int> The minium k-mer occurance')
    validate_parser.add_argument('-o', '--output', type=Path, required=True,
                                 help='<File path> validated DSS result path')
    validate_parser.add_argument('-t', '--tmp', type=Path, default=Path(gettempdir()),
                                 help='<Directory path> Temporary directory path \
                                       <Default system temporary directory>')
    validate_parser.add_argument('--bin', type=Path, default=Path(''),
                                 help='<Directory path> KMC3 exec directory, containing kmc and kmc_dump \
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
