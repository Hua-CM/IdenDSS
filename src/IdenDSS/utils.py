# -*- coding: utf-8 -*-
# @Time : 2021/11/13 10:31
# @Author : Zhongyi Hua
# @FileName: utils.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import logging
import shutil
import os
from sys import platform
from pathlib import Path
import tempfile
from typing import Optional, Union, List
from logging import Logger


class SettingInfo:
    """
    Record parameters for softwares. e.g. threads, execution paths
    """
    def __init__(self,
                 tmp: Optional[Path] = None,
                 bin_dir: Optional[Path] = Path(''),
                 threads: int = 4,
                 logger: Logger = logging.getLogger('IdenDSS')) -> None:
        """
        : param tmp     : The temporary directory path.
        : param blast   : The BLAST execution directory path, contaning blastn, makeblastdb.
        " param threads : The threads number for parallel software.
        """
        self.tmp = Path(tempfile.mktemp(dir=tmp) + str(os.getpid()))
        self.bin_dir = bin_dir
        self.threads = threads
        self.logger = logger
        self._stdout = None
        self._blastn = None
        self._makeblastdb = None
        self._kmc = None
        self._kmcdump = None
        self._kmctools = None

    @staticmethod
    def _func1(unix, win):
        """
        A wrapper to judge platform type
        """
        if platform in ["linux", "linux2", "darwin"]:
            return unix
        if platform == "win32":
            return win

    @property
    def stdout(self) -> Path:
        """
        Generate the equivalent of /dev/stdout on different device.
        """
        if not self._stdout:
            self._stdout = self._func1(Path('/dev/stdout'), Path('CON'))
        return self._stdout

    @property
    def blastn(self) -> Path:
        """
        return blastn exec
        """
        if not self._blastn:
            self._blastn = self._func1(self.bin_dir / 'blastn', self.bin_dir / 'blastn.exe')
        return self._blastn

    @property
    def makeblastdb(self) -> Path:
        """
        return makeblastdb exec
        """
        if not self._makeblastdb:
            self._makeblastdb = self._func1(self.bin_dir / 'makeblastdb',
                                            self.bin_dir / 'makeblastdb.exe')
        return self._makeblastdb

    @property
    def kmc(self) -> Path:
        """
        return kmc exec
        """
        if not self._kmc:
            self._kmc = self._func1(self.bin_dir / 'kmc',
                                    self.bin_dir / 'kmc.exe')
        return self._kmc

    @property
    def kmcdump(self) -> Path:
        """
        return kmc_dump exec
        """
        if not self._kmcdump:
            self._kmcdump = self._func1(self.bin_dir / 'kmc_dump',
                                        self.bin_dir / 'kmc_dump.exe')
        return self._kmcdump
    
    @property
    def kmctools(self) -> Path:
        """
        return kmc_tools exec
        """
        if not self._kmctools:
            self._kmctools = self._func1(self.bin_dir / 'kmc_tools',
                                         self.bin_dir / 'kmc_tools.exe')
        return self._kmctools

    def initiate(self) -> None:
        """
        Make temporary directory
        """
        if not self.tmp.exists():
            self.tmp.mkdir()

    def autoclean(self) -> None:
        """
        Remove temporary directory
        """
        shutil.rmtree(self.tmp)

    @staticmethod
    def clean_file(file: Path) -> None:
        """
        Remove some invisible character for windows (e.g. \\u202A)
        """
        lines = []
        for line in file.read_bytes().decode('utf-8').strip().split('\n'):
            lines.append(line.strip('\u202a\r'))
        file.write_text('\n'.join(lines) + '\n')

class DataInfo:
    """
    Record parameters for data. e.g. meta information, database path
    """
    def __init__(self,
                 in_put: Optional[Path] = None,
                 output: Optional[Path] = None,
                 db: Optional[Union[Path, List[Path]]] = None,
                 db2 : Optional[Union[Path, List[Path]]] = None,
                 length: int = 40,
                 circular: bool = False) -> None:
        """
        main function for identifying DSS using BLAST.
        : param in_put   : Path : fasta file path for *index*; meta information file path for *iden* and *plugin*
        : param output   : Path : Output directory for DSS identification result
        : param db       : Path : fasta database file path
        : param db2      : Path : if a second fasta database is necessary (i.e. validate function)
        : param length   : int  : DSS length. Int
        : param circular : bool : Whether the sequences in database is circular ? e.g. choroplast genome
        """
        self.db  = db
        self.db2 = db2
        self.input = in_put
        self.length = length
        self.circular = circular
        self.output = output

    def inititate(self) -> None:
        """
        Make output directory
        """
        if not self.output.exists():
            self.output.mkdir()


# for RFLP
enzymes_list = [
    'ApaI',
    'BamHI',
    'BglII',
    'EcoRI',
    'HindIII',
    'KpnI',
    'NcoI',
    'NdeI',
    'NheI',
    'NotI',
    'SacI',
    'SalI',
    'SphI',
    'XbaI',
    'XhoI']
