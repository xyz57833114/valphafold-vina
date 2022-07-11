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



_ROOT_MOUNT_DIRECTORY = '/mnt/'

def _create_mount(mount_name: str, path: str) -> Tuple[types.Mount, str]:
  """Create a mount point for each file and directory used by the model."""
  path = pathlib.Path(path).absolute()
  target_path = pathlib.Path(_ROOT_MOUNT_DIRECTORY, mount_name)

  if path.is_dir():
    source_path = path
    mounted_path = target_path
  else:
    source_path = path.parent
    mounted_path = pathlib.Path(target_path, path.name)
  if not source_path.exists():
    raise ValueError(f'Failed to find source directory "{source_path}" to '
                     'mount in Docker container.')
  logging.info('Mounting %s -> %s', source_path, target_path)
  mount = types.Mount(target=str(target_path), source=str(source_path),
                      type='bind', read_only=True)
  return mount, str(mounted_path)

def docker_service(fasta_paths,
                   alphafold_path='/home/xyzhang/alphafold-main/',
                   use_gpu=True,
                   run_relax=True,
                   enable_gpu_relax=True,
                   gpu_devices='all',
                   output_dir='/tmp/alphafold',
                   data_dir='/mnt/nfs02/57t02/af2/download/',
                   docker_image_name='alphafold',
                   max_template_date='2020-05-14',
                   db_preset='full_dbs',
                   model_preset='monomer',
                   num_multimer_predictions_per_model=5,
                   benchmark=False,
                   use_precomputed_msas=False,
                   docker_user=f'{os.geteuid()}:{os.getegid()}'):
    assert model_preset in ['monomer', 'monomer_casp14', 'monomer_ptm', 'multimer'], print(
        'model preset is not in available values')
    assert db_preset in ['full_dbs', 'reduced_dbs'], print(
        'db_preset is not in available values')


    # You can individually override the following paths if you have placed the
    # data in locations other than the FLAGS.data_dir.

    # Path to the Uniref90 database for use by JackHMMER.
    uniref90_database_path = os.path.join(
        data_dir, 'uniref90', 'uniref90.fasta')

    # Path to the Uniprot database for use by JackHMMER.
    uniprot_database_path = os.path.join(
        data_dir, 'uniprot', 'uniprot.fasta')

    # Path to the MGnify database for use by JackHMMER.
    mgnify_database_path = os.path.join(
        data_dir, 'mgnify', 'mgy_clusters_2018_12.fa')

    # Path to the BFD database for use by HHblits.
    bfd_database_path = os.path.join(
        data_dir, 'bfd',
        'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt')

    # Path to the Small BFD database for use by JackHMMER.
    small_bfd_database_path = os.path.join(
        data_dir, 'small_bfd', 'bfd-first_non_consensus_sequences.fasta')

    # Path to the Uniclust30 database for use by HHblits.
    uniclust30_database_path = os.path.join(
        data_dir, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08')

    # Path to the PDB70 database for use by HHsearch.
    pdb70_database_path = os.path.join(data_dir, 'pdb70', 'pdb70')

    # Path to the PDB seqres database for use by hmmsearch.
    pdb_seqres_database_path = os.path.join(
        data_dir, 'pdb_seqres', 'pdb_seqres.txt')

    # Path to a directory with template mmCIF structures, each named <pdb_id>.cif.
    template_mmcif_dir = os.path.join(data_dir, 'pdb_mmcif', 'mmcif_files')

    # Path to a file mapping obsolete PDB IDs to their replacements.
    obsolete_pdbs_path = os.path.join(data_dir, 'pdb_mmcif', 'obsolete.dat')

    data_dir_path = pathlib.Path(data_dir)
    if alphafold_path == data_dir_path or alphafold_path in data_dir_path.parents:
        raise app.UsageError(
            f'The download directory {data_dir} should not be a subdirectory '
            f'in the AlphaFold repository directory. If it is, the Docker build is '
            f'slow since the large databases are copied during the image creation.')

    mounts = []
    command_args = []

    # Mount each fasta path as a unique target directory.
    target_fasta_paths = []
    for i, fasta_path in enumerate(fasta_paths):
        mount, target_path = _create_mount(f'fasta_path_{i}', fasta_path)
        mounts.append(mount)
        target_fasta_paths.append(target_path)
    command_args.append(f'--fasta_paths={",".join(target_fasta_paths)}')

    database_paths = [
        ('uniref90_database_path', uniref90_database_path),
        ('mgnify_database_path', mgnify_database_path),
        ('data_dir', data_dir),
        ('template_mmcif_dir', template_mmcif_dir),
        ('obsolete_pdbs_path', obsolete_pdbs_path),
    ]

    if model_preset == 'multimer':
        database_paths.append(('uniprot_database_path', uniprot_database_path))
        database_paths.append(('pdb_seqres_database_path',
                               pdb_seqres_database_path))
    else:
        database_paths.append(('pdb70_database_path', pdb70_database_path))

    if db_preset == 'reduced_dbs':
        database_paths.append(('small_bfd_database_path', small_bfd_database_path))
    else:
        database_paths.extend([
            ('uniclust30_database_path', uniclust30_database_path),
            ('bfd_database_path', bfd_database_path),
        ])
    for name, path in database_paths:
        if path:
            mount, target_path = _create_mount(name, path)
            mounts.append(mount)
            command_args.append(f'--{name}={target_path}')

    output_target_path = os.path.join(_ROOT_MOUNT_DIRECTORY, 'output')
    mounts.append(types.Mount(output_target_path, output_dir, type='bind'))

    use_gpu_relax = enable_gpu_relax and use_gpu

    command_args.extend([
        f'--output_dir={output_target_path}',
        f'--max_template_date={max_template_date}',
        f'--db_preset={db_preset}',
        f'--model_preset={model_preset}',
        f'--benchmark={benchmark}',
        f'--use_precomputed_msas={use_precomputed_msas}',
        f'--num_multimer_predictions_per_model={num_multimer_predictions_per_model}',
        f'--run_relax={run_relax}',
        f'--use_gpu_relax={use_gpu_relax}',
        '--logtostderr',
    ])

    client = docker.from_env()
    device_requests = [
        docker.types.DeviceRequest(driver='nvidia', capabilities=[['gpu']])
    ] if use_gpu else None

    container = client.containers.run(
        image=docker_image_name,
        command=command_args,
        device_requests=device_requests,
        remove=True,
        detach=True,
        mounts=mounts,
        user=docker_user,
        environment={
            'NVIDIA_VISIBLE_DEVICES': gpu_devices,
            # The following flags allow us to make predictions on proteins that
            # would typically be too long to fit into GPU memory.
            'TF_FORCE_UNIFIED_MEMORY': '1',
            'XLA_PYTHON_CLIENT_MEM_FRACTION': '4.0',
        })

    # Add signal handler to ensure CTRL+C also stops the running container.
    signal.signal(signal.SIGINT,
                  lambda unused_sig, unused_frame: container.kill())

    for line in container.logs(stream=True):
        logging.info(line.strip().decode('utf-8'))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_paths', type=str, help='the fasta_paths file, absolute path without file extension, e.g. /tmp/alphafold/Y265H.fasta')
    parser.add_argument('output_dir', type=str, help='the output_dir file, absolute path without file extension, e.g. /tmp/alphafold')
    args = parser.parse_args()
    docker_service(fasta_paths=[args.fasta_paths], output_dir=args.output_dir)
