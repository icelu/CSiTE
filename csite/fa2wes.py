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
# from csite.phylovar import check_seed, random_int
import shutil
# handle the error below
# python | head == IOError: [Errno 32] Broken pipe
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)


# MAX_INT = 1e8   # CapSim cannot accept a seed that is too large

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

    assert args.wes in ["capsim","wessim","capgem"], "The specified simulator {} is not supported".format(args.wes)


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
    gsize_file = os.path.join(normal_dir, 'normal_gsize.txt')
    if os.path.exists(gsize_file):
        with open(gsize_file,'r') as fin:
            normal_gsize = int(fin.read().strip())
    else:
        for parental in 0, 1:
            normal_gsize += genomesize(
                fasta='{}/normal.parental_{}.fa'.format(normal_dir, parental))
        with open(gsize_file,'w') as fout:
            line = '{}\n'.format(normal_gsize)
            fout.write(line)

    return normal_gsize


def compute_tumor_dna(tumor_dir, tip_node_leaves):
    tumor_dna = 0
    tip_node_gsize = {}
    gsize_file = os.path.join(tumor_dir, 'tumor_gsize.txt')
    if os.path.exists(gsize_file):
        with open(gsize_file,'r') as fin:
            tumor_dna = int(fin.readline().strip())
            for line in fin:
                fields = line.strip().split('\t')
                tip_node_gsize[fields[0]] = [int(fields[1]), int(fields[2]), int(fields[3])]
    else:
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

        with open(gsize_file,'w') as fout:
            line = '{}\n'.format(tumor_dna)
            fout.write(line)
            for key, val in tip_node_gsize.items():
                line = '{}\t{}\t{}\t{}\n'.format(key, val[0], val[1], val[2])
                fout.write(line)

    return tip_node_gsize, tumor_dna


def clean_output(level, outdir):
    '''
    Remove intermediate output of WES simulators according to the specified levels.
    Level 0: keep all the files.
    Level 1: remove "stdout", ".snakemake", "log", "XXX_reads".
    Level 2: keep "config", "genome_index", "mapping", "frags"(capgem), "merged".
    Level 3: keep only "merged".
    '''
    if level == 0:
        return
    elif level == 1:
        # Remove tracking files while running snakemake
        dirs_del = ["stdout", ".snakemake"]
        for entry in os.scandir(outdir):
            if entry.is_dir():
                if entry.name in dirs_del or "reads" in entry.name:
                    shutil.rmtree(entry.path)
    elif level == 2:
        # Used to rerun based on previous mapping results
        dirs_keep = ["config", "genome_index", "mapping", "frags", "merged"]
        for entry in os.scandir(outdir):
            if entry.is_dir():
                if entry.name not in dirs_keep:
                    shutil.rmtree(entry.path)
    elif level == 3:
        # Only keep the final reads
        for entry in os.scandir(outdir):
            if entry.is_dir():
                if entry.name != "merged":
                    shutil.rmtree(entry.path)


def prepare_sample_normal(sample_file, args, normal_gsize):
    '''
    Create a configuration file for running snakemake
    '''
    with open(sample_file, 'w') as fout:
        fout.write('probe: {}\n'.format(args.probe))
        fout.write('directory: normal\n')

        fout.write("genomes:\n")
        # two normal cell haplotypes
        for parental in 0, 1:
            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            fullname = os.path.abspath(ref)
            fout.write("  normal.parental_{}: {}\n".format(parental, fullname))

        fout.write("samples:\n")
        total_num_splits = 0
        # two normal cell haplotypes
        for parental in 0, 1:
            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            proportion = genomesize(fasta=ref) / normal_gsize
            readnum = int((proportion * args.normal_depth *
                       args.target_size) / args.read_length)
            readnum = int(readnum / args.capture_efficiency)

            if readnum > args.max_readnum:
                num_splits = int(numpy.ceil(readnum / args.max_readnum))
                total_num_splits += num_splits
                for split in range(1, num_splits+1):
                    fout.write("  normal.parental_{}_{}:\n".format(parental, str(split)))
                    fout.write('    gid: normal.parental_{}\n'.format(parental))
                    fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                    fout.write('    split: {}\n'.format(str(split)))
                    split_readnum = int(numpy.ceil(readnum/num_splits))
                    fout.write('    readnum: {}\n'.format(str(split_readnum)))
            else:
                total_num_splits += 1
                fout.write("  normal.parental_{}:\n".format(parental))
                fout.write('    gid: normal.parental_{}\n'.format(parental))
                fout.write('    proportion: {}\n'.format(str(proportion)))
                fout.write('    readnum: {}\n'.format(str(readnum)))

    return total_num_splits


def prepare_sample_tumor(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize):
    '''
    Create a configuration file for running snakemake
    '''
    with open(sample_file, 'w') as fout:
        fout.write('probe: {}\n'.format(args.probe))
        fout.write('directory: tumor\n')
        fout.write("genomes:\n")
        # two normal cell haplotypes
        for parental in 0, 1:
            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            fullname = os.path.abspath(ref)
            fout.write("  normal.parental_{}: {}\n".format(parental, fullname))
        # tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            for parental in 0, 1:
                ref = '{}/{}.parental_{}.fa'.format(
                    args.tumor, tip_node, parental)
                fullname = os.path.abspath(ref)
                fout.write('  {}.parental_{}: {}\n'.format(tip_node, parental, fullname))


        fout.write("samples:\n")
        total_num_splits = 0
        # two normal cell haplotypes
        for parental in 0, 1:
            ref = '{}/normal.parental_{}.fa'.format(args.normal, parental)
            fullname = os.path.abspath(ref)
            cell_proportion = normal_cells / total_cells
            proportion = cell_proportion * genomesize(fasta=ref) / normal_gsize
            readnum = int((proportion * args.depth *
                       args.target_size) / args.read_length)
            readnum = int(readnum / args.capture_efficiency)

            if readnum > args.max_readnum:
                num_splits = int(numpy.ceil(readnum / args.max_readnum))
                total_num_splits += num_splits
                for split in range(1, num_splits+1):
                    fout.write("  normal.parental_{}_{}:\n".format(parental, str(split)))
                    fout.write('    gid: normal.parental_{}\n'.format(parental))
                    fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                    fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                    fout.write('    split: {}\n'.format(str(split)))
                    split_readnum = int(numpy.ceil(readnum/num_splits))
                    fout.write('    readnum: {}\n'.format(str(split_readnum)))
            else:
                total_num_splits += 1
                fout.write("  normal.parental_{}:\n".format(parental))
                fout.write('    gid: normal.parental_{}\n'.format(parental))
                fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                fout.write('    proportion: {}\n'.format(str(proportion)))
                fout.write('    readnum: {}\n'.format(str(readnum)))

        # tumor cells haplotypes
        for tip_node in sorted(tip_node_leaves.keys()):
            for parental in 0, 1:
                ref = '{}/{}.parental_{}.fa'.format(
                    args.tumor, tip_node, parental)
                fullname = os.path.abspath(ref)
                cell_proportion = tip_node_leaves[tip_node] / total_cells
                proportion = cell_proportion * tip_node_gsize[tip_node][parental] / tip_node_gsize[tip_node][2]
                readnum = int((proportion * args.depth *
                           args.target_size) / args.read_length)
                readnum = int(readnum / args.capture_efficiency)

                if readnum > args.max_readnum:
                    num_splits = int(numpy.ceil(readnum / args.max_readnum))
                    total_num_splits += num_splits
                    for split in range(1, num_splits+1):
                        fout.write("  {}.parental_{}_{}:\n".format(tip_node, parental, str(split)))
                        fout.write('    gid: {}.parental_{}\n'.format(tip_node, parental))
                        fout.write('    proportion: {}\n'.format(str(proportion/num_splits)))
                        fout.write('    split: {}\n'.format(str(split)))
                        split_readnum = int(numpy.ceil(readnum/num_splits))
                        fout.write('    readnum: {}\n'.format(str(split_readnum)))
                else:
                    total_num_splits += 1
                    fout.write("  {}.parental_{}:\n".format(tip_node, parental))
                    fout.write('    gid: {}.parental_{}\n'.format(tip_node, parental))
                    fout.write('    cell_proportion: {}\n'.format(str(cell_proportion)))
                    fout.write('    proportion: {}\n'.format(str(proportion)))
                    fout.write('    readnum: {}\n'.format(str(readnum)))
    return total_num_splits



def run_snakemake(outdir, args, jobs, sample_file, snake_file):
    stddir = os.path.join(outdir, 'stdout')
    if not os.path.exists(stddir):
        os.makedirs(stddir)

    if args.wes == "capsim":
        snake_file_copy = os.path.join(outdir, 'config/Snakefile_capsim')
    elif args.wes == "wessim":
        snake_file_copy = os.path.join(outdir, 'config/Snakefile_wessim')
    else:
        snake_file_copy = os.path.join(outdir, 'config/Snakefile')
    # Copy Snakefile  to the output folder
    shutil.copyfile(snake_file, snake_file_copy)

    orig_params = args.snakemake.split()
    config = ''
    if '--config' in orig_params:
        cfg_index = orig_params.index('--config')
        config += orig_params[cfg_index+1]
        del orig_params[cfg_index]
        del orig_params[cfg_index]
    config += ' rlen=' + str(args.read_length)
     # + ' seed=' + str(numpy.random.randint(MAX_INT))
    if args.use_cluster:
        cluster_file = os.path.join(os.path.dirname(sys.argv[0]), 'config/cluster.yaml')
        assert os.path.isfile(cluster_file), 'Cannot find cluster.yaml under the program directory'
        cluster_file_copy = os.path.join(outdir, 'config/cluster.yaml')
        shutil.copyfile(cluster_file, cluster_file_copy)
        cluster = '\"qsub -V -l mem_free={cluster.mem},h_rt={cluster.time} -pe OpenMP {cluster.n} -o ' + os.path.abspath(stddir) + ' -e '  + os.path.abspath(stddir) +'\"'
        final_cmd_params =  orig_params + ['-s', os.path.abspath(snake_file_copy), '-d', os.path.abspath(outdir), '--cluster-config', cluster_file_copy,  '--cluster', cluster, '--configfile', os.path.abspath(sample_file), '--jobs', str(jobs),  '--config', config]
    else:
        final_cmd_params =  orig_params + ['-s', os.path.abspath(snake_file_copy), '-d', os.path.abspath(outdir), '--configfile', os.path.abspath(sample_file), '--jobs', str(jobs),  '--config', config]
    logging.info(' Command: %s', ' '.join(final_cmd_params))

    os.system(' '.join(final_cmd_params))
    clean_output(args.out_level, outdir)



def main(progname=None):
    parser = argparse.ArgumentParser(
        description='a wrapper of simulating targeted capture sequencing from reference genome files',
        prog=progname if progname else sys.argv[0])

    group1 = parser.add_argument_group('Input options')
    group1.add_argument('-n', '--normal', metavar='DIR', required=True,
                       help='The directory of the fasta files of normal genomes')
    group1.add_argument('-t', '--tumor', metavar='DIR', required=True,
                       help='The directory of the fasta files of tumor genomes')
    group1.add_argument('-c', '--chain', metavar='DIR', required=True,
                       help='The directory of the tumor chain files')
    group1.add_argument( '--probe', metavar='FILE', required=True,
                       help='The file containing the sequences of target region')

    group2 = parser.add_argument_group('Parameters for sequencing')
    # group = group2.add_mutually_exclusive_group()
    default = 100
    group2.add_argument('-d', '--depth', metavar='FLOAT', type=float, default=default,
                       help='The mean depth of tumor for simulating short reads [{}]'.format(default))
    default = 100
    group2.add_argument('-D', '--normal_depth', metavar='FLOAT', type=float, default=default,
                       help='The mean depth of normal for simulating short reads [{}]'.format(default))
    default = 0.5
    group2.add_argument('-p', '--purity', metavar='FLOAT', type=float, default=default,
                       help='The proportion of tumor cells in simulated sample [{}]'.format(default))
    default = 150
    group2.add_argument('--read_length', metavar='INT', type=int, default=default,
                       help='Illumina: read length [{}]'.format(default))

    # TODO: Comute target size from probe sequences
    default = 51189318
    group2.add_argument('--target_size', metavar='INT', type=int, default=default,
                       help='The size of target regions for simulating short reads [{}]'.format(default))
    default = 0.5
    group2.add_argument('--capture_efficiency', metavar='FLOAT', type=float, default=default,
                       help='The capture efficiency of the capture kit [{}]'.format(default))
    default = 1e6
    group2.add_argument('--max_readnum', metavar='INT', type=int, default=default,
                       help='The number of maximum short reads [{}] for a single run of simulation'.format(default))
    # default = None
    # group2.add_argument('-s', '--random_seed', type=check_seed,
    #                    help='The seed for random number generator [{}]'.format(default))
    default = 'wessim'
    group2.add_argument('--wes', nargs='?', default=default, choices=['wessim', 'capsim', 'capgem'],
                       help='The whole-exome sequencing simulator used for simulating short reads [{}]'.format(default))
    default = 'snakemake --rerun-incomplete -k --latency-wait 120 --config fmedian=500'
    group2.add_argument('--snakemake', metavar='STR', type=str, default=default,
                       help='The command used for calling a whole-exome sequencing simulator [{}]'.format(default))
    default = False
    group2.add_argument('--use_cluster', action='store_true', default=default,
                   help='Run simulation on a cluster [{}]'.format(default))

    group3 = parser.add_argument_group('Output options')
    default = 'wes_reads'
    group3.add_argument('-o', '--output', metavar='DIR', type=str, default=default,
                       help='The output directory [{}]'.format(default))
    default = 'fa2wes.log'
    group3.add_argument('-g', '--log', metavar='FILE', type=str, default=default,
                       help='The log file to save the settings of each command [{}]'.format(default))
    default = 1
    group3.add_argument('--out_level', metavar='INT', type=int, default=default,
                       help='The level used to indicate how many intermediate output files are kept [{}]. \
                       Level 0: keep all the files.\
                           Level 1: remove "stdout", ".snakemake", "log", "XXX_reads". \
                           Level 2: keep "config", "genome_index", "mapping", "frags"(capgem), "merged".\
                           Level 3: keep only "merged".'.format(default))

    args = parser.parse_args()

    # logging and random seed setting
    logging.basicConfig(filename=args.log,
                        filemode='w', format='[%(asctime)s] %(levelname)s: %(message)s',
                        datefmt='%m-%d %H:%M:%S', level='INFO')
    logging.info(' Command: %s', ' '.join(sys.argv))

    # if args.random_seed == None:
    #     seed = random_int()
    # else:
    #     seed = args.random_seed
    # logging.info(' Random seed: %s', seed)
    # numpy.random.seed(seed)

    check_input(args)

    # Add folder wes to system path for simulation
    wes_dir = os.path.join(os.path.dirname(sys.argv[0]), 'wes')
    os.environ["PATH"] += os.pathsep + wes_dir

    if args.wes == "capsim":
        snake_file = os.path.join(os.path.dirname(sys.argv[0]), 'config/Snakefile_capsim')
    elif args.wes == "wessim":
        snake_file = os.path.join(os.path.dirname(sys.argv[0]), 'config/Snakefile_wessim')
        wessim_dir = os.path.join(wes_dir, 'wessim')
        os.environ["PATH"] += os.pathsep + wessim_dir
        file_path = os.path.join(wessim_dir, 'src/Wessim2.py')
        subprocess.call(['chmod', '0755', file_path])
    else:
        snake_file = os.path.join(os.path.dirname(sys.argv[0]), 'config/Snakefile')
        capgem_dir = os.path.join(wes_dir, 'capgem')
        os.environ["PATH"] += os.pathsep + os.path.join(capgem_dir, 'bin')
        os.environ["PATH"] += os.pathsep + os.path.join(capgem_dir, 'src')
        file_path = os.path.join(capgem_dir, 'src/frag2read.py')
        subprocess.call(['chmod', '0755', file_path])

    assert os.path.isfile(snake_file), 'Cannot find Snakefile under the program directory'

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
        total_num_splits = prepare_sample_tumor(sample_file, args, total_cells, normal_cells, normal_gsize, tip_node_leaves, tip_node_gsize)

        # jobs = 4 * (len(tip_node_leaves) * 2 + 2)
        jobs = total_num_splits
        run_snakemake(outdir, args, jobs, sample_file, snake_file)

    if args.normal_depth > 0:
        outdir = os.path.join(args.output, 'normal')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        configdir = os.path.join(outdir, 'config')
        if not os.path.exists(configdir):
            os.makedirs(configdir)

        sample_file = os.path.join(outdir, 'config/sample.yaml')
        total_num_splits = prepare_sample_normal(sample_file, args, normal_gsize)

        jobs = total_num_splits

        run_snakemake(outdir, args, jobs, sample_file, snake_file)
