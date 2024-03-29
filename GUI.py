# -*- coding: utf-8 -*-
# @Time : 2020/10/04 21:00
# @Author : Zhongyi Hua
# @FileName: GUI.py
# @Usage: The GUI interface for IdenDSS
# @Note: 
# @E-mail: njbxhzy@hotmail.com

from pathlib import Path
import logging
import PySimpleGUI as sg

from IdenDSS.identification import iden_main
from IdenDSS.plugins import IndexDb, search_rflp, combine_dss, summary_dss, flanking, convert
from IdenDSS.validation import v_main
from IdenDSS.utils import enzymes_list, SettingInfo, DataInfo


def index_tab():
    """
    The tab for indexing database
    """
    row1 = sg.Frame('Data setting', [
        [sg.T('Input path:',  size=15), sg.I(key='-INDEX INPUT-', size=45), sg.FileBrowse(target='-INDEX INPUT-')],
        [sg.T('Output path:', size=15), sg.I(key='-INDEX OUTPUT-', size=45), sg.FolderBrowse(target='-INDEX OUTPUT-')],
        [sg.T('DSS length:', size=15), sg.I('40', key='-INDEX LENGTH-', size=5)]
    ])
    row2 = sg.Frame('Global setting', [
        [sg.T('Bin directory path (Optional):', size=30), sg.I(key='-INDEX BIN_DIR-'), sg.FolderBrowse(target='-INDEX BIN_DIR-')]
    ])
    main_col = [
        [row1],
        [row2]
    ]
    return sg.Tab('index', layout=main_col, expand_x = True)


def iden_tab():
    """
    The tab for identifying DSS
    """
    row1 = sg.Frame('Data setting', [
        [sg.T('Input path:',  size=15), sg.I(key='-IDEN INPUT-', size=45), sg.FileBrowse(target='-IDEN INPUT-')],
        [sg.T('Output path:', size=15), sg.I(key='-IDEN OUTPUT-', size=45), sg.FolderBrowse(target='-IDEN OUTPUT-')],
        [sg.T('Database path:', size=15), sg.I(key='-IDEN DATABASE-', size=45), sg.FileBrowse(target='-IDEN DATABASE-')],
        [sg.T('DSS length:', size=15), sg.I('40', key='-IDEN LENGTH-', size=5)]
    ])
    row2 = sg.Frame('Global setting', [
        [sg.T('Temporary directory path (Optional):', size=30), sg.I(key='-IDEN TMP_DIR-'), sg.FolderBrowse(target='-IDEN TMP_DIR-')],
        [sg.T('Bin directory path (Optional):', size=30), sg.I(key='-IDEN BIN_DIR-'), sg.FolderBrowse(target='-IDEN BIN_DIR-')],
        [sg.T('Threads  (Optional):', size=30), sg.I('4', key='-IDEN THREADS-')]
    ])
    main_col = [
        [row1],
        [row2]
    ]
    return sg.Tab('identify', layout=main_col, expand_x = True)


def plugin_tab():
    """
    The tab for plugin
    """
    row1 = sg.Frame('Option', [
        [sg.Checkbox('Combine', k='-COMBINE-'),
         sg.Checkbox('RFLP', k='-RFLP-'),
         sg.Checkbox('Statistic', k='-STATISTIC-'),
         sg.Checkbox('Convert', k='-CONVERT-')],
        [sg.Checkbox('Flanking', k='-FLANK-'),  sg.I('200', key='-FLANK VALUE-')],
        [sg.Checkbox('Sample', k='-SAMPLE-'),  sg.I('5', key='-SAMPLE VALUE-')]
    ])
    row2 = sg.Frame('Data setting', [
        [sg.T('Input path:',  size=15), sg.I(key='-PLUG INPUT-', size=45), sg.FileBrowse(target='-PLUG INPUT-')],
        [sg.T('Output path:', size=15), sg.I(key='-PLUG OUTPUT-', size=45), sg.FolderBrowse(target='-PLUG OUTPUT-')],
        [sg.T('Database path:', size=15), sg.I(key='-PLUG DATABASE-', size=45), sg.FileBrowse(target='-PLUG DATABASE-')]
    ])
    row3 = sg.Frame('Global setting', [
        [sg.T('Temporary directory path (Optional):', size=30), sg.I(key='-PLUG TMP_DIR-'), sg.FolderBrowse(target='-PLUG TMP_DIR-')]
    ])
    main_col = [
        [row1],
        [row2],
        [row3]
    ]
    return sg.Tab('plugin', layout=main_col, expand_x = True)


def validate_tab():
    """
    The tab for validation
    """
    row1 = sg.Frame('Data setting', [
        [sg.T('Input path:',  size=18),
         sg.I(key='-VALI INPUT-', size=42),
         sg.FileBrowse(target='-VALI INPUT-')],
        [sg.T('Species fastq path:', size=18),
         sg.ML(key='-VALI SP-', size=(40,2)),
         sg.I(visible=False, key='-VALI SP MULTIIN-', enable_events=True),
         sg.FilesBrowse(target='-VALI SP MULTIIN-', file_types=(('fastq', '.fastq .fq'),))],
        [sg.T('Background fastq path:', size=18),
         sg.ML(key='-VALI BG-', size=(40,2)),
         sg.I(visible=False, key='-VALI BG MULTIIN-', enable_events=True),
         sg.FilesBrowse(target='-VALI BG MULTIIN-', file_types=(('fastq', '.fastq .fq'),))],
        [sg.T('Output path:', size=18),
         sg.I(key='-VALI OUTPUT-', size=42),
         sg.FolderBrowse(target='-VALI OUTPUT-')]
    ])
    row2 = sg.Frame('Global setting', [
        [sg.T('Temporary directory path (Optional):', size=30), sg.I(key='-VALI TMP_DIR-'), sg.FolderBrowse(target='-VALI TMP_DIR-')],
        [sg.T('Bin directory path (Optional):', size=30), sg.I(key='-VALI BIN_DIR-'), sg.FolderBrowse(target='-VALI BIN_DIR-')]
    ])
    main_col = [
        [row1],
        [row2]
    ]
    return sg.Tab('validate', layout=main_col, expand_x = True)

def make_window():
    """
    Make the main window
    """
    sg.theme('SystemDefaultForReal')
    tab1 = index_tab()
    tab2 = iden_tab()
    tab3 = plugin_tab()
    tab4 = validate_tab()

    layout = [
        [sg.TabGroup([[tab1, tab2, tab3, tab4]], key='-TASK-')],
        [sg.CB('This is a circular DNA sequence Database', k='-CIRCULAR-')],
        [sg.ML('',
               size=(60,5),
               k='-OUT-',
               write_only=True,
               expand_x=True,
               reroute_stdout=True,
               reroute_stderr=True,
               echo_stdout_stderr=True,
               reroute_cprint=True,
               auto_refresh=True,
               autoscroll=True)],
        [sg.Button('Run', key='run'), sg.Exit()]
        ]
    window = sg.Window('IdenDSS', layout, default_element_size=(30,1))
    return window


def main():
    """
    The main interface for event loop
    """
    window = make_window()
    logger = logging.getLogger('IdenDSS')
    logger.setLevel(level=logging.INFO)
    # Setting a special Handler for Logger to work with PySimpleGUI is necessary
    view_handler = logging.StreamHandler()
    logger.addHandler(view_handler)
    while True:
        event, values = window.read()
        if event in (sg.WINDOW_CLOSED, 'Exit'):
            break
        if event.endswith('MULTIIN-'):
            # for multiline box input. use "-* MULTIIN-" for invisible input box.
            multi_key = ' '.join(event.split(' ')[:-1]) + '-'
            window[multi_key].Update('\n'.join(values[event].split(';')) + '\n', append=True)
        if event == 'run':
            if values['-TASK-'] == 'index':
                set_info = SettingInfo(
                    bin_dir=Path(values['-INDEX BIN_DIR-']),
                    logger=logger)
                data_info = DataInfo(
                    in_put=Path(values['-INDEX INPUT-']),
                    output=Path(values['-INDEX OUTPUT-']),
                    circular=values['-CIRCULAR-'],
                    length=int(values['-INDEX LENGTH-'])
                )
                index_db = IndexDb(data_info, set_info)
                index_db.index()
            elif values['-TASK-'] == 'identify':
                set_info = SettingInfo(
                    tmp=Path(values['-IDEN TMP_DIR-']),
                    bin_dir=Path(values['-IDEN BIN_DIR-']),
                    threads=Path(values['-IDEN THREADS-'],
                    logger=logger)
                )
                data_info = DataInfo(
                    in_put=Path(values['-IDEN INPUT-']),
                    output=Path(values['-IDEN OUTPUT-']),
                    db=Path(values['-IDEN DATABASE-']),
                    length=int(values['-IDEN LENGTH-']),
                    circular=values['-CIRCULAR-']
                )
                iden_main(data_info, set_info)
            elif values['-TASK-'] == 'plugin':
                set_info = SettingInfo(
                    tmp=Path(values['-PLUG TMP_DIR-']),
                    logger=logger)
                set_info.clean_file(Path(values['-PLUG INPUT-']))
                data_info = DataInfo(
                    in_put=Path(values['-PLUG INPUT-']),
                    output=Path(values['-PLUG OUTPUT-']),
                    db=Path(values['-PLUG DATABASE-']),
                    circular=values['-CIRCULAR-'])
                if values['-RFLP-']:
                    search_rflp(data_info, set_info, enzymes_list)
                if values['-COMBINE-']:
                    combine_dss(data_info, set_info)
                if values['-STATISTIC-']:
                    summary_dss(data_info, set_info)
                if values['-FLANK-']:
                    data_info.length = int(values['-FLANK VALUE'])
                    flanking(data_info, set_info)
                if  values['-CONVERT-']:
                    convert(data_info, set_info)
                if  values['-SAMPLE-']:
                    data_info.length = int(values['-SAMPLE VALUE'])
                    convert(data_info, set_info)
            elif values['-TASK-'] == 'validate':
                set_info = SettingInfo(
                    tmp=Path(values['-VALI TMP_DIR-']),
                    bin_dir=Path(values['-VALI BIN_DIR-']),
                    logger=logger)
                data_info = DataInfo(in_put=Path(values['-VALI INPUT-']),
                             output=Path(values['-VALI OUTPUT-']),
                             db=[Path(_) for _ in values['-VALI BG-'].split('\n')],
                             db2=[Path(_) for _ in values['-VALI SP-'].split('\n')])
                v_main(data_info, set_info)
    window.close()

if __name__ == '__main__':
    main()
