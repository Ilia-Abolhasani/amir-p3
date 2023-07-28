#!/bin/bash
directory="./data/rfam"
download() {
    accession="$1"    
    wget "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/${accession}.fa.gz" -P "${directory}"    
    gunzip "${directory}/${accession}.fa.gz"
}

# Create the target directory if it doesn't exist
mkdir -p "${directory}"


### accession, ID, Type, description
# RF00001, 5S_rRNA, rRNA, 5S ribosomal RNA
download "RF00001" 
# RF00002, 5_8S_rRNA, rRNA, 5.8S ribosomal RNA
download "RF00002" 
# RF00005, tRNA, tRNA, tRNA
download "RF00005" 
# RF01960, SSU_rRNA_eukarya, rRNA, Eukaryotic small subunit ribosomal RNA
download "RF01960" 
# RF02543, LSU_rRNA_eukarya, rRNA, Eukaryotic large subunit ribosomal RNA
download "RF02543"
