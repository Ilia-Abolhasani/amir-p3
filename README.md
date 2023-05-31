# AmiR-P3: Ab Initio Plant miRNA Prediction Pipeline

AmiR-P3 is a powerful ab initio plant miRNA prediction pipeline written in Python 3. It leverages various computational techniques and integrates several tools to predict putative miRNAs from genomic or transcriptomic sequences in plants.

## Background

MicroRNAs (miRNAs) are small noncoding RNAs that play important post-transcriptional regulatory roles in animals and plants. However, the inherent complexity of miRNA biogenesis in plants hampers the application of standard miRNA prediction tools, which are often optimized for animal sequences. Therefore, computational approaches to predict putative miRNAs from genomic sequences in plants, regardless of their expression levels or tissue specificity, are of great interest.

## Features

- AmiR-P3 is a novel ab initio plant miRNA prediction pipeline.
- Leverages the strengths of various utilities for its key computational steps.
- Allows users to adjust the prediction criteria based on the state-of-the-art biological knowledge of plant miRNA properties.
- Starts with finding the potential homologs of known plant miRNAs in the input sequence(s) and ensures they do not overlap with protein-coding regions.
- Computes the minimum free energy structure of the presumed RNA sequence.
- Uses a pre-trained deep learning classification model to predict potential miRNAs based on the computed structures.
- Applies a set of criteria to select the most likely miRNAs from the set of predicted miRNAs.
- Yields acceptable predictions in a variety of plant species.

## Installation

To use AmiR-P3, it is recommended to use the provided Docker image. The Docker image contains all the necessary software, making the installation process straightforward. Follow these steps to install and run AmiR-P3 using Docker:

1. Install Docker on your system by following the official Docker installation guide for your operating system.

2. Pull the AmiR-P3 Docker image from Docker Hub by running the following command:

   ```shell
      docker pull micrornaproject/amir-p3
   ```
3. Once the image is downloaded, you can run the AmiR-P3 pipeline using the following command:
   ```shell
      docker run -v /path/to/input:/data -v /path/to/output:/output docker.com/amir-p3 python3 amiR-P3.py -i /data/input_sequence.fasta -o /output/predicted_miRNAs.fasta
   ```
   Make sure to replace /path/to/input and /path/to/output with the actual paths to your input sequence file and desired output directory, respectively. The input_sequence.fasta should contain the genomic sequence(s) from which you want to predict miRNAs.
   
   Alternatively, you can connect to the AmiR-P3 Docker container in interactive mode, which allows you to manually execute commands and interact with the pipeline. Follow the instructions below to use the interactive mode:
   
   ```shell
   sudo docker run -it micrornaproject/amir-p3:latest /bin/bash   
   ```
   ```shell
   # Run the AmiR-P3 pipeline inside the Docker container
   python3 amiR-P3.py --input ./data/example/example_genome.fasta --experiment test
   ```
4. AmiR-P3 will generate the predicted miRNAs in the specified output directory.

Note: The Docker image contains all the necessary software, except for mxfold2, which needs to be manually installed if selected as the second structure prediction software. To install mxfold2 manually, follow these steps:

1. Download mxfold2:
   ```shell
   wget https://github.com/keio-bioinformatics/mxfold2/releases/download/v0.1.1/mxfold2-0.1.1.tar.gz
2. Install mxfold2:
   ```shell
   pip3 install mxfold2-0.1.1.tar.gz

## NR Database Download

To use the Diamond or Blastx tools for finding non-coding sequences, you will need to download the NR (non-redundant) proteins database. Follow these steps to download the NR database:

1. Visit the following link: [NR Database Download](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz)

2. Download the `nr.gz` file to your local machine.

3. Extract the contents of the `nr.gz` file. You can use any suitable tool or command for extracting gzip archives.

After extracting the database, you will need to provide the path to the NR database file when running the AmiR-P3 pipeline using the `--nr` command-line argument.

Please note that the NR database is a large file, and the download process may take some time depending on your internet connection speed.

## Command-Line Arguments

The AmiR-P3 pipeline supports the following command-line arguments:

- `--input`: Path to the input genome (required).
- `--experiment`: Experiment name (required).
- `--ncpu`: Number of CPU cores to use (default: 4).
- `--nc`: Number of nonconformity (default: 3).
- `--fv`: Number of flanking value (default: 200).
- `--ssm`: Secondary structure method (default: viennarna). Options: mfold, viennarna, contrafold, mxfold2.
- `--ft`: Folding temperature (default: 22).
- `--ss`: Start seed position (default: 2).
- `--se`: Seed seed position (default: 13).
- `--ht`: Hit threshold Jaccard similarity (default: 0.8).
- `--pt`: Precursor threshold Jaccard similarity (default: 0.8).
- `--bt`: BOI threshold Jaccard similarity (default: 0.8).
- `--pce`: Protein coding elimination (default: True).
- `--pcem`: Protein coding elimination method (default: diamond). Options: diamond, blastx.
- `--nr`: Path to the RefSeq non-redundant proteins database.
- `--diamonddb`: Path to the diamond database.

Please note that some arguments have default values and are not required. However, the `--input` and `--experiment` arguments are mandatory for running the pipeline.

## Citation

If you use AmiR-P3 in your research or find it useful, please consider citing: []

## License
AmiR-P3 is released under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC-BY-NC-SA 4.0).

This means that you are free to use, share, and adapt AmiR-P3 for non-commercial purposes, as long as you provide attribution, use the same license for any modifications, and do not use it for commercial purposes.

To view a copy of the license, please visit [License](LICENSE.md).

If you are interested in using AmiR-P3 for commercial purposes, please contact us for further licensing options.

## Contact Information

For any questions, issues, or feedback, please feel free to contact the developer:

- [Ilia-Abolhasani](https://github.com/Ilia-Abolhasani/)

