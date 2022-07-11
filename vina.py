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

    :param receptor: 输入的是蛋白序列文件的名称，无后缀名
    :param out_dir:输出的路径
    :return:
    """
    file_name = os.path.join(out_dir, receptor, 'ranked_0.pdb')  # f'/tmp/alphafold/{receptor}/ranked_0.pdb'
    while not os.path.exists(file_name):  # 判断文件是否存在
        time.sleep(0.5)
    print('alphaflod finish!')
    f = open(file_name, 'r')
    #读取每一行数据
    lines = f.readlines()
    new_lines = lines[:-1]
    f_w = open(os.path.join(out_dir, f'{receptor}.pdbqt'), "w")
    for line in new_lines:
        f_w.write(line)
    f.close()
    f_w.close()


def autodock_vina_run(receptor_file, ligand_file, out_file, log_file):
    """

    :param receptor_file:输入的是蛋白序列文件的名称，有后缀名, 绝对路径
    :param ligand_file:输入的是化合物文件的名称，有后缀名, 绝对路径
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

def pdb_to_pdbqt_autodock_vina_run(receptor_pdb, ligand_file, out_file, log_file):

    pdb_to_pdbqt(receptor_pdb, receptor_pdbqt)
    autodock_vina_run(receptor_pdbqt, ligand_file, out_file, log_file)
    return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('receptor_file', type=str, help='the receptor file, absolute path without file extension, e.g. /tmp/autodock_vina/000001.pdbqt')
    parser.add_argument('ligand_file', type=str, help='the ligand file, absolute path without file extension, e.g. /tmp/autodock_vina/1.pdbqt')
    parser.add_argument('out_file', type=str, help='the format of out file, e.g. /tmp/autodock_vina/out.pdbqt')
    parser.add_argument('log_file', type=str, help='the dictionary of log file, e.g. /tmp/autodock_vina/log.txt')
    args = parser.parse_args()
    autodock_vina_run(args.receptor_file, args.ligand_file, args.out_file, args.log_file)

    
    