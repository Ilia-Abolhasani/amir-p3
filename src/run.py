#!/usr/bin/env python
import os
import sys
import argparse
from enum import Enum


class SecondaryStructure(Enum):
    mfold = 'mfold'
    viennarna = 'viennarna'
    contrafold = 'contrafold'
    mxfold2 = 'mxfold2'

    def __str__(self):
        return self.value


class ProteinCodingEliminationMethod(Enum):
    dimond = 'dimond'
    blastx = 'blastx'

    def __str__(self):
        return self.value


program_name = "AmiR-P3"

description = """This is a miRNA prediction pipline.
Input should be a FASTA files"""

parser = argparse.ArgumentParser(prog=program_name,
                                 description=description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=False)

# arguments
parser.add_argument('--input',
                    type=str,
                    required=True,
                    help='Path to input genome')


parser.add_argument('--ncpu',
                    default=4,
                    type=int,
                    required=False,
                    help='Number of cpu cores')


parser.add_argument('--nc',
                    default=3,
                    type=int,
                    required=False,
                    help='Number of nonconformity')


parser.add_argument('--fv',
                    default=200,
                    type=int,
                    required=False,
                    help='Number of flanking value')


parser.add_argument('--ss',
                    default=2,
                    type=int,
                    required=False,
                    help='Start seed position.')


parser.add_argument('--se',
                    default=13,
                    type=int,
                    required=False,
                    help='Eed seed position.')


parser.add_argument('--ht',
                    default=0.8,
                    type=float,
                    required=False,
                    help='Hit threshold jaccard similarity')


parser.add_argument('--pt',
                    default=0.8,
                    type=float,
                    required=False,
                    help='Precursor threshold jaccard similarity')


parser.add_argument('--bt',
                    default=0.8,
                    type=float,
                    required=False,
                    help='BOI threshold jaccard similarity')


parser.add_argument('--ssm',
                    default=SecondaryStructure.viennarna,
                    type=SecondaryStructure,
                    choices=list(SecondaryStructure),
                    required=False,
                    help='Secondary structure method')


parser.add_argument('--pce',
                    default=True,
                    type=bool,
                    required=False,
                    help='Protein coding elimination')


parser.add_argument('--pcem',
                    default=ProteinCodingEliminationMethod.dimond,
                    type=ProteinCodingEliminationMethod,
                    choices=list(ProteinCodingEliminationMethod),
                    required=False,
                    help='Protein coding elimination method')


parser.add_argument('--ft',
                    default=22,
                    type=int,
                    required=False,
                    help='Folding temperature')

parser.add_argument('--nr',
                    default=None,
                    type=int,
                    required=False,
                    help='RefSeq non-redundant proteins database path')


if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
    parser.print_help()
    sys.exit()

# config
experiment = "A.thaliana"
input_genome_name = "GCF_000001735.4_TAIR10.1_genomic.fna"
experiment_dir = "./Experiment"
mirbase_dir = "./miRBase"
secondary_structure_method = "mfold"


input_genome_path = f"{experiment_dir}/{experiment}/{input_genome_name}"
temp_path = f"{experiment_dir}/{experiment}/Temp"
result_path = f"{experiment_dir}/{experiment}/Result"
current_path = os.getcwd()
