# -*- coding: utf-8 -*-
# @Time : 2020/11/1 19:18
# @Author : Zhongyi Hua
# @FileName: IDSS.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com


import tempfile, os
from src.probe_utils import iden_main
from src.primer_utils import primer_main
from src.plugins import Database, search_rflp, combine
import argparse


def getArgs():
    parser = argparse.ArgumentParser(
        prog='DSS identification', description='This script was for identifying DNA signature sequences(DSS)')
    sub_parser = parser.add_subparsers(title='Available', dest='database/iden')
    sub_parser.required = True

    database_parser = sub_parser.add_parser(
        'database', help='Generate database for DSS identification')
    database_parser.add_argument('-i', '--input', required=True,
                                 help='<file_path>  The genome fasta file')
    database_parser.add_argument('-l', '--length', type=int, default=20,
                                 help='<int> DSS length <default=20>')
    database_parser.add_argument('-c', '--circular',  action='store_true',
                                 help='the sequences are circular or not')                           
    database_parser.add_argument('-o', '--output', required=True,
                                 help='<file_path>  The genome fasta output file (a BLAST database would '
                                      'be constructed)')
    database_parser.add_argument('-b', '--blast', default='',
                                 help='<directory path> BLAST exec directory <If your BLAST software not in PATH>')                                  
    database_parser.set_defaults(subcmd="database")

    probe_parser = sub_parser.add_parser(
        'iden', help='Identification DSS based on database')
    probe_parser.add_argument('-m', '--meta', required=True,
                              help='<file path> The meta file')
    probe_parser.add_argument('-d', '--database', required=True,
                              help='<file path> Database fasta <Corresponding BLAST database file must in the same '
                                   'directory>')
    probe_parser.add_argument('-l', '--length', type=int, default=20,
                              help='<int> DSS length <default=20>')
    probe_parser.add_argument('-c', '--circular',  action='store_true',
                              help='the sequences are circular or not')
    probe_parser.add_argument('-@', '--threads', type=int, default=4,
                              help='<int> Threads used in BLAST <Default 4>')
    probe_parser.add_argument('-t', '--tmp', default=None,
                              help='<dir path> Temporary directory path <Default system temporary directory>')
    probe_parser.add_argument('-o', '--output', required=True,
                              help='<directory path> result directory')
    probe_parser.add_argument('-b', '--blast', default='',
                              help='<directory path> BLAST exec directory <If your BLAST software not in PATH>')
    probe_parser.set_defaults(subcmd="iden")

    plugin_parser = sub_parser.add_parser(
        'plugin', help='Some plugin for DSS results')
    plugin_parser.add_argument('-c', '--circular',  action='store_true',
                               help='the sequences are circular or not')
    plugin_parser.add_argument('-i', '--input', required=True,
                               help='<file_path>  A meta file. One DSS result file path per line (combined result is OK)')
    plugin_parser.add_argument('-o', '--output', required=True,
                               help='<directory path> result directory')
    plugin_parser.add_argument('-d', '--database', required=True,
                               help='<file path> Database fasta (The database used to identify DSS)')
    plugin_parser.add_argument('-t', '--tmp', default=None,
                               help='<dir path> Temporary directory path <Default system temporary directory>')
    plugin_parser.add_argument('--primer', action='store_true',
                               help='Design Primer (Need primer3_core in your PATH)')
    plugin_parser.add_argument('--rflp', action='store_true',
                               help='Identify restriction enzyme sits on DSS for putative RFLP method')
    plugin_parser.add_argument('--combine', action='store_true',
                               help='Generate the corresponding combined DSS file')
    plugin_parser.add_argument('--statistic', action='store_true',
                               help='Count the DSS number')
    plugin_parser.set_defaults(subcmd="plugin")

    args = parser.parse_args()
    return args


def main():
    args = getArgs()

    if args.subcmd == "database":
        db = Database(args.input, args.output, args.blast)
        if args.circular:
            db.database_generate(args.length)
        else:
            db.copy_file()
        db.database_blast()

    else:
        if not os.path.exists(args.output):
            os.mkdir(args.output)
            #args.tmp = tempfile.mktemp(dir=args.tmp)

        if args.subcmd == "iden":
            # set temporary directory
            args.tmp = tempfile.mktemp(dir=args.tmp)
            iden_main(args)

        elif args.subcmd == "plugin":
            args.tmp = tempfile.mktemp(dir=args.tmp)
            if args.primer:
                primer_main(args)
            if args.rflp:
                search_rflp(args)
            if args.combine:
                combine(args)


if __name__ == '__main__':
    main()
