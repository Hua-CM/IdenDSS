# -*- coding: utf-8 -*-
# @Time : 2020/11/1 19:18
# @Author : Zhongyi Hua
# @FileName: identification_dss.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import tempfile
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from probe_utils import iden_main
import argparse


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


def getArgs():
    parser = argparse.ArgumentParser(
        prog='DSS identification', description='This script was for identifying DNA signature sequences(DSS)')
    sub_parser = parser.add_subparsers(title='Available', dest='database/iden')
    sub_parser.required = True

    database_parser = sub_parser.add_parser(
        'database', help='Generate database for DSS identification')
    database_parser.add_argument('-i', '--input_fasta', required=True,
                                 help='<file_path>  The genome fasta file')
    database_parser.add_argument('-l', '--length', type=int, default=20,
                                 help='<int> DSS length <default=20>')
    database_parser.add_argument('-c', '--circular',  action='store_true',
                                 help='the sequences are circular or not')                           
    database_parser.add_argument('-o', '--output_fasta', required=True,
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
    args = parser.parse_args()
    return args


def main():
    args = getArgs()
    if args.subcmd == "database":
        db = Database(args.input_fasta, args.output_fasta, args.blast)
        if args.circular:
            db.database_generate(args.length)
        else:
            db.copy_file()
        db.database_blast()

    elif args.subcmd == "iden":
        # set temporary directory
        args.tmp = tempfile.mktemp(dir=args.tmp)
        iden_main(args)


if __name__ == '__main__':
    main()
