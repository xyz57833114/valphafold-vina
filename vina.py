import os
import pathlib
import signal
import time
import sys
import getopt
from typing import Tuple

from absl import app
from absl import logging
import docker
from docker import types



def pdb_to_pdbqt(receptor, out_dir):
    """
    将pdb文件(f'{out_dir}/{receptor}.pdb')转化为pdbqt格式并保存为f'{out_dir}/{receptor}.pdbqt'
    :param receptor: 输入的是蛋白序列文件的名称，无后缀名
    :param out_dir:输出的路径
    :return:
    """
    file_name = os.path.join(out_dir, receptor)  # f'/tmp/alphafold/{receptor}/ranked_0.pdb'
    f = open(file_name, 'r')
    #读取每一行数据
    lines = f.readlines()
    new_lines = lines[:-1]
    f_w = open(os.path.join(out_dir, f'{receptor}.pdbqt'), "w")
    for line in new_lines:
        f_w.write(line)
    f.close()
    f_w.close()

def openbabel(ligand_file, format, out_file):
    """
    :ligand_file:为输入化合物文件名，无后缀名
    :format：为输入化合物文件格式
    """
    #past_file =
    cmd = f'obabel -i {format} {ligand_file}.{format} -opdbqt -O {out_file}'
    return os.popen(cmd, 'r')

def autodock_vina_run(receptor_file, ligand_file, out_file, log_file):
    """
    用于给蛋白质pdbqt格式和化合物pdbqt格式做分子对接 autodock vina
    :param receptor_file:输入的是蛋白序列文件的名称，有后缀名, 绝对路径
    :param ligand_file:输入的是化合物文件的名称，又后缀名, 绝对路径
    :param out_file:输出的对接结果文件，为pdbqt格式
    :param log_file：输出的对接结果打分值，为txt文件
    :return:
    """
    # log_file = '/tmp/autodock_vina/log.txt'
    # out_file = '/tmp/autodock_vina/out.pdbqt'
    # ligand_file = '/tmp/autodock_vina/1.pdbqt'
    # receptor_file = '/tmp/autodock_vina/000001.pdbqt'
    config_file = '/tmp/autodock_vina/config.txt'
    cmd = f'/home/xyzhang/autodock/autodock_vina_1_1_2_linux_x86/bin/vina --config {config_file}' \
        f' --receptor {receptor_file} --ligand {ligand_file} --out {out_file} --log {log_file}'
    # print(cmd)
    return os.popen(cmd, 'r')


def openbabel_vina(receptor, ligand_file, format, out_dir, job_id=0):
    """
    用于给不同格式的化合物进行格式转换
    :param receptor: 输入的是蛋白序列文件的名称，无后缀名
    :param ligand_file: 输入的是化合物文件的名称，无后缀名
    :param format：为输入化合物文件格式
    :param out_dir:输出的路径
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    _, ligand_name = os.path.split(ligand_file)
    _, receptor_name = os.path.split(receptor)
    openbabel(f'{ligand_file}',
              f'{format}',
              os.path.join(out_dir, f'{ligand_name}.pdbqt'))
    pdb_to_pdbqt(receptor, out_dir)
    # receptor_name = receptor_name.replace('.pdbqt', '')
    # ligand_name = ligand_name.replace('.pdbqt', '')
    autodock_vina_run(f'{receptor}.pdbqt',
                      f'{ligand_file}.pdbqt',
                      os.path.join(out_dir, f'{receptor_name}_{ligand_name}.pdbqt'),
                      os.path.join(out_dir, f'{receptor_name}_{ligand_name}.txt'))
    file_name = os.path.join(out_dir, f'{receptor_name}_{ligand_name}.txt')
    while not os.path.exists(file_name):  # 判断文件是否存在
        time.sleep(0.5)
    print('vina finish!')

    generate_html(outdir=out_dir,
                  html_name=f'{receptor_name}_{ligand_name}.html',
                  job_id=job_id,
                  protein_tertiary_structure_PDBQT_format_file=f'{receptor_name}.pdbqt',
                  compound_PDBQT_format_file=f'{ligand_name}.pdbqt',
                  autodock_vina_molecular_docking_result=f'{receptor_name}_{ligand_name}.pdbqt',
                  autodock_vina_molecular_docking_scoring_value=f'{receptor_name}_{ligand_name}.txt'
                  )


def generate_html(outdir,
                  html_name,
                  job_id=0,
                  protein_tertiary_structure_PDBQT_format_file='Y265H/ranked_0.pdb',
                  compound_PDBQT_format_file='1.pdbqt',
                  autodock_vina_molecular_docking_result='Y265H_1.pdbqt',
                  autodock_vina_molecular_docking_scoring_value='Y265H_1.txt'
                  ):
    html_name = os.path.join(outdir, html_name)
    f = open(html_name, 'w')

    message = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>REPORT</title>
        </head>
        <body>
            <h2 align="left">JOB id: %d</h2>
            <h2></h2>
            <h2 align="left">Input file</h2>
            <p align="left">You can find your input fasta file here.</p>
            <h3></h3>
            <h2 align="left">Result files</h2>

            <h3 align="left">Protein tertiary structure PDBQT format file</h3>
            <p align="left">Click <a href="%s">here</a> to download the Protein tertiary structure PDBQT format file.</p>
            <h3 align="left">Compound PDBQT format file</h3>
            <p align="left">Click <a href="%s">here</a> to download the PPDBQT format file.</p>
            <h3 align="left">Autodock vina molecular docking result</h3>
            <p align="left">Click <a href="%s">here</a> to download the Autodock vina molecular docking result file.</p>
            <h3 align="left">Autodock Vina molecular docking scoring value</h3>
            <p align="left">Click <a href="%s">here</a> to download the Autodock Vina molecular docking scoring value file.</p>

        </body>
        </html>
        """ % (job_id,
               protein_tertiary_structure_PDBQT_format_file,
               compound_PDBQT_format_file,
               autodock_vina_molecular_docking_result,
               autodock_vina_molecular_docking_scoring_value)

    f.write(message)
    f.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('receptor', type=str, help='the receptor file, absolute path without file extension, e.g. /tmp/alphafold/Y265H')
    parser.add_argument('ligand', type=str, help='the ligand file, absolute path without file extension, e.g. /tmp/alphafold/1')
    parser.add_argument('formate', type=str, help='the format of ligand file, e.g. mol2')
    parser.add_argument('outdir', type=str, help='the dictionary of output files, e.g. /tmp/alphafold/Y265H_1')
    args = parser.parse_args()
    openbabel_vina(args.receptor, args.ligand, args.formate, args.outdir)

