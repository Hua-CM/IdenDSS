# -*- coding: utf-8 -*-
# @Time : 2022/9/23 17:49
# @Author : Zhongyi Hua
# @FileName: validation.py
# @Usage: Using fastq file validate putative DSS
# @Note: dependent on KMC
# @E-mail: njbxhzy@hotmail.com

from pathlib import Path
from typing import Optional, List
import subprocess
import tempfile

"""
In my experience, the KMC indexing a 500 Mbp FASTQ file need ~5 GB memory and ~0.5 min, and
the raw output file is ~200 MB with an ~1.8 GB kmc_dump file. The KMC indexing a 2.5 Gbp
FASTQ file need ~11 Gb and ~2min (The default max amount of RAM for KMC is 12 GB), and the
raw output file is ~3.2 GB.

**The reverse complement sequence have been counted in KMC results.**
"""


class ValiDSS:
    """
    Validate one DSS file using one background speceis FASTQ file and
    one intraspecies FASTQ file.
    """
    def __init__(self, setinfo, datainfo) -> None:
        """
        : param pdss    : set  : putative DSS sequences set, excluding reverse complement sequences.
        : param bg data : Path : The background species HTS dataset
        : param sp_data : Path : The intraspecies HTS dataset
        """
        self.setinfo = setinfo
        self.datainfo = datainfo
        self.run_tmp = Path(tempfile.mktemp(dir=setinfo.tmp) )

    def _gen_dss(self, prefix_dss: str='dss') -> str:
        _fa_path = self.setinfo.tmp / 'dss.fa'
        _lines   = self.datainfo.input.read_text().strip().split('\n')[1:]
        out_lst  = [ f'>seq{_idx}\n' + _line.split('\t')[2]
                   for _idx, _line in enumerate(_lines)]
        _fa_path.write_text('\n'.join(out_lst) + '\n')
        # generate KMC database
        self._run_kmc(_fa_path, prefix_dss, ['-ci1', '-fa'])
        return prefix_dss


    def _gen_bg(self, prefix_bg: str = 'kmc_bg') -> str:
        # for background dataset
        out_file = self.setinfo.tmp / 'background_fq.lst'
        out_lines = []
        for _path in self.datainfo.db:
            out_lines.append(str(_path.resolve()))
        out_file.write_text('\n'.join(out_lines)+'\n')
        self._run_kmc(f'@{str(out_file.resolve())}', prefix_bg)
        return prefix_bg

    def _gen_sp(self, prefix_sp: str = 'kmc_sp') -> str:
        out_file = self.setinfo.tmp / 'species_fq.lst'
        out_lines = []
        for _path in self.datainfo.db2:
            out_lines.append(str(_path.resolve()))
        out_file.write_text('\n'.join(out_lines)+'\n')
        self._run_kmc(f'@{str(out_file.resolve())}', prefix_sp)
        return prefix_sp

    def _run_kmc(self, db_path: Path, _res_prefix: Path, extra_lst: Optional[List] = None) -> None:
        # subprocess run must be a list of **str**!
        _res_path = self.setinfo.tmp / _res_prefix
        arg_lst = [f'-t{self.setinfo.threads}', f'-k{self.datainfo.length}']
        if extra_lst:
            arg_lst += extra_lst
        subprocess.run(
            [str(self.setinfo.kmc)] +
            arg_lst +
            [str(db_path),
             str(_res_path),
             str(self.run_tmp)],
             check=True
        )

    def _run_intersect_kmc(self, _prefix1: str, _prefix2: str) -> str:
        #kmc_tools simple test.kmc kmc_tmp_res intersect test_bg.res
        subprocess.run(
            [str(self.setinfo.kmctools),
             'simple',
             str(self.setinfo.tmp / _prefix1),
             str(self.setinfo.tmp / _prefix2),
             'intersect',
             str(self.setinfo.tmp / (_prefix1 + '_' + _prefix2))],
             check=True)
        return _prefix1 + '_' + _prefix2

    def _run_get_seq_kmc(self, _prefix:str) -> set:
        out_set = set()
        subprocess.run(
            [str(self.setinfo.kmcdump),
             str(self.setinfo.tmp / _prefix),
             str(self.setinfo.tmp / (_prefix + '_tmp'))],
            check=True
        )
        for line in (self.setinfo.tmp / (_prefix + '_tmp')).read_bytes().strip().split(b'\n'):
            out_set.add(line.split(b'\t')[0].decode("utf-8") )
        return out_set

    def validate(self) -> set:
        """
        The public method interface for DSS validation
        """
        self.run_tmp.mkdir()
        # background dataset
        self.setinfo.logger.info('Validate DSS in background dataset.')
        prefix_bg = self._gen_bg()
        # intraspecies dataset
        self.setinfo.logger.info('Validate DSS in intraspecies dataset.')
        prefix_sp = self._gen_sp()
        # DSS dataset
        prefix_dss = self._gen_dss()
        # intersect
        self.setinfo.logger.info('Get intersection results.')
        prefix_bg_dss = self._run_intersect_kmc(prefix_dss, prefix_bg)
        prefix_sp_dss = self._run_intersect_kmc(prefix_dss, prefix_sp)
        # get result
        pdss_bg = self._run_get_seq_kmc(prefix_bg_dss)
        pdss_sp = self._run_get_seq_kmc(prefix_sp_dss)
        self.run_tmp.rmdir()
        self.setinfo.logger.info('Confirm results. Done')
        return pdss_sp - pdss_bg


def v_main(datainfo, setinfo):
    """
    The main interface for DSS validation
    : param datainfo :
        : input  : Path       : The input DSS result path.
        : output : Path       : The output validated DSS result path.
        : db     : List[Path] : The background HTS FASTQ file.
        : db2    : List[Path] : The intraspecies species HTS FASTQ file.
    : param setinfo  :
        : tmp     : The temporary directory path.
        : kmc     : The kmc bin path.
        : kmcdump : The kmc_dump bin path.
        : threads : Parallel threads for KMC3.
    : NOTE : The DSS need generate from the input file !!!
    """
    setinfo.initiate()
    try:
        lines = datainfo.input.read_text().strip().split('\n')
        datainfo.length = len(lines[1].split('\t')[2])
        dss_valiation = ValiDSS(setinfo, datainfo)
        dss_set = dss_valiation.validate()
        out_lines = [lines[0]]
        for line in lines[1:]:
            if line.split('\t')[2] in dss_set:
                out_lines.append(line)
        datainfo.output.write_text('\n'.join(out_lines) + '\n')
    except Exception as e:
        setinfo.logger.error(e)
    finally:
        # remove the tmp directory
        setinfo.autoclean()
