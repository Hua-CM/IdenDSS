# -*- coding: utf-8 -*-
# @Time : 2020/11/1 19:18
# @Author : Zhongyi Hua
# @FileName: identification_dss.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import tempfile
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline


class GenerateProbe:
    # !!!! notice the key and value in assembly dict
    def __init__(self, _seq_obj):
        """
        :param _seq_obj: a Bio.Seq object
        """
        self.seq = _seq_obj
        self.assembly_dict = {}

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


class BLASTParse:
    def __init__(self, _seq_path, _asm_list, _database_path, _tmp_dir_path, _length, _thread, _exec=''):
        self.seq_path = _seq_path
        self.asm_list = _asm_list
        self.db_path = _database_path
        self.tmp_dir = _tmp_dir_path
        self.length = _length
        self.thread = _thread
        self.exec = _exec

    def blast_cmd(self):
        _group = os.path.splitext(os.path.split(self.seq_path)[1])[0]
        blastn_cmd = NcbiblastnCommandline(
            cmd=os.path.join(self.exec, 'blastn'),
            query=self.seq_path,
            db=self.db_path,
            task='blastn-short',
            dust='no',
            word_size=min(round(self.length/2)-1, 11),
            outfmt='\"6 qacc sacc length pident evalue\"',
            num_threads=self.thread,
            evalue=10,
            out=os.path.join(self.tmp_dir, _group + '.txt'),
            max_hsps=1,
            max_target_seqs=len(self.asm_list)+1)
        blastn_cmd()

    def blast_parse(self):
        _putative_kmer = []
        _group = os.path.splitext(os.path.split(self.seq_path)[1])[0]
        _blast_df = pd.read_table(os.path.join(self.tmp_dir, _group + '.txt'),
                                  names=['query', 'subject', 'length', 'identity', 'evalue'])
        _blast_df = _blast_df[(_blast_df['length'] == self.length) & (_blast_df['identity'] == 100)]
        # Find DSS **in** other species
        _blast_df1 = _blast_df.loc[~(_blast_df['subject'].isin(self.asm_list))].drop_duplicates(subset=['query'])
        # Find DSS conserved in species
        _blast_df2 = _blast_df[_blast_df['subject'].isin(self.asm_list)]
        _blast_df2 = _blast_df2.groupby('query')['subject'].count().reset_index()
        _blast_df2 = _blast_df2[_blast_df2['subject'] == len(self.asm_list)]
        # Find intersection as species DSS
        _dss_list = list(set(_blast_df2['query'].to_list()) - set(_blast_df1['query'].to_list()))
        # DSS set empty
        if not _dss_list:
            _result_tb = pd.DataFrame.from_dict({0: {'seq':  '', 'position': '', 'GC': ''}}, 'index')
        else:
            _kmer_dict = {_.id: _.seq for _ in
                          SeqIO.parse(self.seq_path, 'fasta')
                          }
            _tmp_list = []
            for _dss in _dss_list:
                _tmp_list.append(
                                 {'seq': _kmer_dict[_dss],
                                  'position': _dss.split('_')[-1] +
                                              '-' +
                                              str(int(_dss.split('_')[-1])+self.length-1),
                                  'GC': GC(_kmer_dict[_dss])}
                                 )
            _result_tb = pd.DataFrame(_tmp_list)
        return _result_tb


class Database:
    def __init__(self, fa_path, db_path, exec):
        """
        Since chloroplast was a circular, we need cut the first XX bp to the end
        :param fa_path:
        :param db_path:
        """
        self.in_path = fa_path
        self.out_path = db_path
        self.exec = exec

    def database_generate(self, _length=20):
        fasta_in = SeqIO.parse(self.in_path, 'fasta')
        fasta_out = []
        for assembly in fasta_in:
            assembly.seq = assembly.seq + assembly.seq[0:_length-1]
            fasta_out.append(assembly)
        SeqIO.write(fasta_out, self.out_path, 'fasta')

    def database_blast(self):
        database_cmd = NcbimakeblastdbCommandline(
            cmd=os.path.join(self.exec, 'makeblastdb'),
            dbtype='nucl',
            input_file=self.out_path)
        database_cmd()


if __name__ == '__main__':
    import argparse
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
    database_parser.add_argument('-o', '--output_fasta', required=True,
                                 help='<file_path>  The genome fasta output file (a BLAST database would '
                                      'be constructed)')
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

    if args.subcmd == "database":
        db = Database(args.input_fasta, args.output_fasta, args.blast)
        db.database_generate(args.length)
        db.database_blast()

    elif args.subcmd == "iden":
        seq_dict = {_.id: _ for _ in SeqIO.parse(args.database, 'fasta')}
        meta_info = pd.read_table(args.meta, names=['group', 'assembly'])
        meta_info = meta_info.groupby('group')['assembly'].apply(list).reset_index(name='assembly')
        used_assembly_list = [row.assembly[0] for row in meta_info.itertuples()]
        used_seq_dict = {}
        for asm in used_assembly_list:
            try:
                used_seq_dict[asm] = seq_dict[asm]
            except KeyError:
                raise AttributeError('Please Check whether', asm, 'is in the database')
        # release memory
        del seq_dict
        # set temporary directory
        args.tmp = tempfile.mktemp(dir=args.tmp)
        os.makedirs(args.tmp)
 
        for _idx, _row in meta_info.iterrows():
            pb = GenerateProbe(used_seq_dict[_row['assembly'][0]])
            pb.probe_generate(args.length)
            pb.probe_save(os.path.join(args.tmp, _row['group'] + '.fasta'))
            ap = BLASTParse(
                os.path.join(args.tmp, _row['group'] + '.fasta'),
                _row['assembly'],
                args.database,
                args.tmp,
                args.length,
                args.threads,
                args.blast
            )
            ap.blast_cmd()
            result_tb = ap.blast_parse()
            result_tb['assembly'] = _row['assembly'][0]
            result_tb['group'] = _row['group']
            result_tb[['group', 'assembly', 'seq', 'position', 'GC']].\
                to_csv(os.path.join(args.output, _row['group'] + '.txt'),
                       sep='\t',
                       index=False)
            for f in args.tmp:
                os.remove(f)
        