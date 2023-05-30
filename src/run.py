#!/usr/bin/env python
import os
import sys
import argparse
import pipeline
from enum import Enum


class SecondaryStructure(Enum):
    mfold = 'mfold'
    viennarna = 'viennarna'
    contrafold = 'contrafold'
    mxfold2 = 'mxfold2'

    def __str__(self):
        return self.value


class ProteinCodingEliminationMethod(Enum):
    diamond = 'diamond'
    blastx = 'blastx'

    def __str__(self):
        return self.value


program_name = "AmiR-P3"

description = "AmiR-P3 is an advanced ab initio plant miRNA prediction pipeline in Python, offering customizable prediction criteria and leveraging powerful computational techniques."

parser = argparse.ArgumentParser(prog=program_name,
                                 description=description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=False)

# arguments
parser.add_argument('--input',
                    type=str,
                    required=True,
                    help='Path to the input genome file in FASTA format.')

parser.add_argument('--experiment',
                    type=str,
                    required=True,
                    help='Name of the experiment or analysis.')

parser.add_argument('--ncpu',
                    default=4,
                    type=int,
                    required=False,
                    help='Number of CPU cores to be used for the computation.')


parser.add_argument('--nc',
                    default=3,
                    type=int,
                    required=False,
                    help='Number of nonconformity.')


parser.add_argument('--fv',
                    default=200,
                    type=int,
                    required=False,
                    help='Number of flanking value.')


parser.add_argument('--ssm',
                    default=SecondaryStructure.viennarna,
                    type=SecondaryStructure,
                    choices=list(SecondaryStructure),
                    required=False,
                    help='Secondary structure prediction method.')


parser.add_argument('--ft',
                    default=22,
                    type=int,
                    required=False,
                    help='Folding temperature in degrees. Default is 22.')


parser.add_argument('--ss',
                    default=2,
                    type=int,
                    required=False,
                    help='Start seed position.')


parser.add_argument('--se',
                    default=13,
                    type=int,
                    required=False,
                    help='End seed position.')


parser.add_argument('--ht',
                    default=0.8,
                    type=float,
                    required=False,
                    help='Hit threshold for Jaccard similarity. Default is 0.8.')


parser.add_argument('--pt',
                    default=0.8,
                    type=float,
                    required=False,
                    help='Precursor threshold for Jaccard similarity. Default is 0.8.')


parser.add_argument('--bt',
                    default=0.8,
                    type=float,
                    required=False,
                    help='BOI (Branch of Interest) threshold for Jaccard similarity. Default is 0.8.')


parser.add_argument('--pce',
                    default=False,
                    type=bool,
                    required=False,
                    help='Enable or disable protein coding elimination. Default is False.')


parser.add_argument('--pcem',
                    default=ProteinCodingEliminationMethod.diamond,
                    type=ProteinCodingEliminationMethod,
                    choices=list(ProteinCodingEliminationMethod),
                    required=False,
                    help='Protein coding elimination method. Available options: diamond, blastx. Default is diamond.')


parser.add_argument('--nr',
                    default=None,
                    type=str,
                    required=False,
                    help='Path to the RefSeq non-redundant proteins database.')


parser.add_argument('--diamonddb',
                    default=None,
                    type=str,
                    required=False,
                    help='Path to the Diamond database.')


if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
    parser.print_help()
    sys.exit()

args = parser.parse_args()

pipeline.start(
    input_genome_path=args.input,
    experiment=args.experiment,
    experiment_dir="./experiment",
    mirbase_dir="./data/miRBase",
    nonconformity=args.nc,
    flanking_value=args.fv,
    folding_temperature=args.ft,
    seed_start=args.ss,
    seed_end=args.se,
    hit_threshold=args.ht,
    precursor_threshold=args.pt,
    boi_threshold=args.bt,
    num_cpus=args.ncpu,
    secondary_structure_method=args.ssm,
    apply_protein_coding_elimination=args.pce,
    protein_coding_elimination_method=args.pcem,
    diamond_nr_path=args.nr,
    diamond_db_path=args.diamonddb,
)
