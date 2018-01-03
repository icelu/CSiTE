#!/usr/bin/env python3

#########################################################################
# Author: Bingxin Lu
# Created Time: 2017-12-18
# File Name: fa2wes.py
# Description: Simulate WES reads from whole genome sequences
#########################################################################

import sys
import os
import argparse
import numpy
import logging
import pyfaidx
import subprocess
from csite.phylovar import check_seed, random_int
import shutil
# handle the error below
# python | head == IOError: [Errno 32] Broken pipe
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)


MAX_INT = 1e8   # CapSim cannot accept a seed that is too large

def check_input(args):
    # there must be two haplotype fasta in the normal dir
    assert os.path.isdir(
        args.normal), '{} is not exist or not a folder.'.format(args.normal)
    for parental in 0, 1:
        assert os.path.isfile('{}/normal.parental_{}.fa'.format(args.normal, parental)),\
            'Can not find normal.parental_{}.fa under the normal directory: {}'.format(
                parental, args.normal)

# tumor directory and chain directory must exist.
# also file chain_dir/tip_node_sample.count.
    assert os.path.isdir(
        args.chain), '{} is not exist or not a folder.'.format(args.chain)
    assert os.path.isfile(args.chain + '/tip_node_sample.count'),\
        'Can not find tip_node_sample.count under the chain directory: {}'.format(
            args.chain)
    assert os.path.isdir(
        args.tumor), '{} is not exist or not a folder.'.format(args.chain)


def tip_node_leaves_counting(f=None):
    '''
    Return a dictionay with structure:
    {tip_node1:leaves_count1,tip_node2:leaves_count2,...}
    '''
    tip_node_leaves = {}
    with open(f, 'r') as input:
        for line in input:
            if not line.startswith('#'):
                tip_node, leaves = line.split()
                tip_node_leaves[tip_node] = int(leaves)
    return tip_node_leaves


# TODO: Extract the size of exome
def genomesize(fasta=None):
    '''
    Extract genome size from .fa file.
    '''
    fa = pyfaidx.Faidx(fasta)
    gsize = 0
    for chroms in fa.index.keys():
        gsize += fa.index[chroms].rlen
    return gsize


def compute_normal_gsize(normal_dir):
    normal_gsize = 0
    for parental in 0, 1:
        normal_gsize += genomesize(
            fasta='{}/normal.parental_{}.fa'.format(normal_dir, parental))
    return normal_gsize


def compute_tumor_dna(tumor_dir, tip_node_leaves):
    tumor_dna = 0
    tip_node_gsize = {}
    for tip_node, leaves in tip_node_leaves.items():
        # The value of tip_node_gsize[tip_node] is a list of three elements:
        # 0)genomesize of parental 0
        # 1)genomesize of parental 1
        # 2)the sum of parental 0 and 1
        tip_node_gsize[tip_node] = []

        for parental in 0, 1:
            assert os.path.isfile('{}/{}.parental_{}.fa'.format(tumor_dir, tip_node, parental)),\
                'Can not find {}.parental_{}.fa under the tumor directory: {}'.format(
                    tip_node, parental, tumor_dir)
            tip_node_gsize[tip_node].append(genomesize(
                fasta='{}/{}.parental_{}.fa'.format(tumor_dir, tip_node, parental)))

        tip_node_gsize[tip_node].append(
            tip_node_gsize[tip_node][0] + tip_node_gsize[tip_node][1])
        tumor_dna += tip_node_gsize[tip_node][2] * tip_node_leaves[tip_node]

    return tip_node_gsize, tumor_dna


def prepare_sample_normal(sample_file, args, normal_gsize):
    '''
    Create a configuration file for running snakemake
    '''
    with open(sample_file, 'w') as fout:
        fout.write('directory: normal\n')
        fout.write("samples:\n")

        # two normal cell haplotypes
        for parental in 0, 1:
            fout.write("  normal.parental_{}:\n".format(parental))

            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            fullname = os.path.abspath(ref)
            fout.write('    genome: {}\n'.format(fullname))

            proportion = genomesize(fasta=ref) / normal_gsize
            fout.write('    proportion: {}\n'.format(str(proportion)))

            readnum = int((proportion * args.normal_depth *
                       args.target_size) / args.read_length)
            readnum = int(readnum / args.capture_efficiency)
            fout.write('    readnum: {}\n'.format(str(readnum)))


def prepare_sample_tumor(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize):
    '''
    Create a configuration file for running snakemake
    '''
    with open(sample_file, 'w') as fout:
        fout.write('directory: tumor\n')
        fout.write("samples:\n")

        # two normal cell haplotypes
        for parental in 0, 1:
            fout.write("  normal.parental_{}:\n".format(parental))

            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            fullname = os.path.abspath(ref)
            fout.write('    genome: {}\n'.format(fullname))

            cell_proportion = normal_cells / total_cells
            fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
            proportion = cell_proportion * genomesize(fasta=ref) / normal_gsize
            fout.write('    proportion: {}\n'.format(str(proportion)))

            readnum = int((proportion * args.depth *
                       args.target_size) / args.read_length)
            readnum = int(readnum / args.capture_efficiency)
            fout.write('    readnum: {}\n'.format(str(readnum)))

        # tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            for parental in 0, 1:
                fout.write("  {}.parental_{}:\n".format(tip_node, parental))

                ref = '{}/{}.parental_{}.fa'.format(
                    args.tumor, tip_node, parental)
                fullname = os.path.abspath(ref)
                fout.write('    genome: {}\n'.format(fullname))

                cell_proportion = tip_node_leaves[tip_node] / total_cells
                fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                proportion = cell_proportion * tip_node_gsize[tip_node][parental] / tip_node_gsize[tip_node][2]
                fout.write('    proportion: {}\n'.format(str(proportion)))

                readnum = int((proportion * args.depth *
                           args.target_size) / args.read_length)
                readnum = int(readnum / args.capture_efficiency)
                fout.write('    readnum: {}\n'.format(str(readnum)))

def run_snakemake(outdir, args, jobs, sample_file, snake_file, cluster_file):
    stddir = os.path.join(outdir, 'stdout')
    if not os.path.exists(stddir):
        os.makedirs(stddir)

    snake_file_copy = os.path.join(outdir, 'config/Snakefile')
    cluster_file_copy = os.path.join(outdir, 'config/cluster.yaml')
    # Copy Snakefile and cluster.yaml to the output folder
    shutil.copyfile(snake_file, snake_file_copy)
    shutil.copyfile(cluster_file, cluster_file_copy)

    orig_params = args.capsim.split()
    config = ''
    if '--config' in orig_params:
        cfg_index = orig_params.index('--config')
        config += orig_params[cfg_index+1]
        del orig_params[cfg_index]
        del orig_params[cfg_index]
    config += ' rlen=' + str(args.read_length) + ' seed=' + str(numpy.random.randint(MAX_INT))
    cluster = '\"qsub -V -l mem_free={cluster.mem},h_rt={cluster.time} -pe OpenMP {cluster.n} -o ' + os.path.abspath(stddir) + ' -e '  + os.path.abspath(stddir) +'\"'

    final_cmd_params =  orig_params + ['-s', os.path.abspath(snake_file_copy), '-d', os.path.abspath(outdir), '--cluster-config', cluster_file_copy,  '--cluster', cluster, '--configfile', os.path.abspath(sample_file), '--jobs', str(jobs),  '--config', config]
    logging.info(' Command: %s', ' '.join(final_cmd_params))

    os.system(' '.join(final_cmd_params))



def main(progname=None):
    parse = argparse.ArgumentParser(
        description='a wrapper of simulating targeted capture sequencing from reference genome files',
        prog=progname if progname else sys.argv[0])
    parse.add_argument('-n', '--normal', required=True,
                       help='the directory of the normal fasta')
    parse.add_argument('-t', '--tumor', required=True,
                       help='the directory of the tumor fasta')
    parse.add_argument('-c', '--chain', required=True,
                       help='the directory of the tumor chain')
    default = 'wes_reads'
    parse.add_argument('-o', '--output', type=str, default=default,
                       help='output directory [{}]'.format(default))
    group = parse.add_mutually_exclusive_group()
    default = 0
    group.add_argument('-d', '--depth', type=float, default=default,
                       help='the mean depth of tumor for simulating short reads [{}]'.format(default))
    default = 0
    group.add_argument('-D', '--normal_depth', type=float, default=default,
                       help='the mean depth of normal for simulating short reads [{}]'.format(default))
    default = 51189318
    parse.add_argument('--target_size', type=int, default=default,
                       help='the size of target regions for simulating short reads [{}]'.format(default))
    default = 0.5
    parse.add_argument('--capture_efficiency', type=float, default=default,
                       help='the capture efficiency of the capture kit [{}]'.format(default))
    default = 150
    parse.add_argument('--read_length', type=int, default=default,
                       help='Illumina: read length [{}]'.format(default))
    default = 0.5
    parse.add_argument('-p', '--purity', type=float, default=default,
                       help='the proportion of tumor cells in simulated sample [{}]'.format(default))
    default = None
    parse.add_argument('-s', '--random_seed', type=check_seed,
                       help='the seed for random number generator [{}]'.format(default))
    default = 'fa2wes.log'
    parse.add_argument('-g', '--log', type=str, default=default,
                       help='the log file to save the settings of each command [{}]'.format(default))
    default = 'snakemake --rerun-incomplete -k --latency-wait 120 --config fmedian=500'
    parse.add_argument('--capsim', type=str, default=default,
                       help='the parameters for calling CapSim program [{}]'.format(default))
    args = parse.parse_args()

    # logging and random seed setting
    logging.basicConfig(filename=args.log,
                        filemode='w', format='[%(asctime)s] %(levelname)s: %(message)s',
                        datefmt='%m-%d %H:%M:%S', level='INFO')
    logging.info(' Command: %s', ' '.join(sys.argv))
    if args.random_seed == None:
        seed = random_int()
    else:
        seed = args.random_seed
    logging.info(' Random seed: %s', seed)
    numpy.random.seed(seed)

    check_input(args)

    snake_file = os.path.join(os.path.dirname(sys.argv[0]), 'config/Snakefile')
    assert os.path.isfile(snake_file), 'Cannot find Snakefile under the program directory'
    cluster_file = os.path.join(os.path.dirname(sys.argv[0]), 'config/cluster.yaml')
    assert os.path.isfile(cluster_file), 'Cannot find cluster.yaml under the program directory'

    normal_gsize = compute_normal_gsize(args.normal)

    # Separate the simulation of tumor and normal samples
    if args.depth > 0:
        outdir = os.path.join(args.output, 'tumor')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        configdir = os.path.join(outdir, 'config')
        if not os.path.exists(configdir):
            os.makedirs(configdir)

        tip_node_leaves = tip_node_leaves_counting(
            f='{}/tip_node_sample.count'.format(args.chain))
        tumor_cells = sum(tip_node_leaves.values())
        total_cells = tumor_cells / args.purity
        logging.info(' Number of total cells: %d', total_cells)
        normal_cells = total_cells - tumor_cells
        logging.info(' Number of normal cells: %d', normal_cells)
        normal_dna = normal_gsize * normal_cells
        tip_node_gsize, tumor_dna= compute_tumor_dna(args.tumor, tip_node_leaves)
        total_dna = (normal_dna + tumor_dna)

        sample_file = os.path.join(outdir, 'config/sample.yaml')
        prepare_sample_tumor(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize)

        jobs = 4 * (len(tip_node_leaves) * 2 + 2)

        run_snakemake(outdir, args, jobs, sample_file, snake_file, cluster_file)

    if args.normal_depth > 0:
        outdir = os.path.join(args.output, 'normal')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        configdir = os.path.join(outdir, 'config')
        if not os.path.exists(configdir):
            os.makedirs(configdir)

        sample_file = os.path.join(outdir, 'config/sample.yaml')
        prepare_sample_normal(sample_file, args, normal_gsize)

        jobs = 8

        run_snakemake(outdir, args, jobs, sample_file, snake_file, cluster_file)
