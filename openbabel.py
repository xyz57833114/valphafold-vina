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


def openbabel(ligand_file, format, out_file):
    """
    :ligand_file:为输入化合物文件名，无后缀名
    :format：为输入化合物文件格式
    """
    #past_file = 
    cmd = f'obabel -i {format} {ligand_file}.{format} -o pdbqt -O {out_file}'
    return os.popen(cmd, 'r')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('ligand_file', type=str, help='the receptor file, absolute path without file extension, e.g. /tmp/alphafold/1')
    parser.add_argument('formate', type=str, help='the format of ligand file, e.g. mol2')
    parser.add_argument('out_file', type=str, help='the dictionary of output files, e.g. /tmp/alphafold/1.pdbqt')
    args = parser.parse_args()
    openbabel(args.ligand_file, args.formate, args.out_file)

    
    