# This will create a CSV file with the paths to the SP GPSC fasta files

#### for reads 

#!/bin/bash

# Directory where GPSCs are stored
GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/TB_reps/"

# Output CSV file
output_file="TB_reps.tsv"

# Write the header line
echo -e "sample_id\tseq_path" > $output_file

# Loop over the fasta files in GPSCs dir
for fasta_file in $GPSC_dir*.fa; do
    # Extract the sample_id from the file name
    sample_id=$(basename $fasta_file .fa)
    # Write a line to the CSV file with the path to the fasta file
    echo -e "$sample_id\t$fasta_file" >> $output_file
done


#### for contigs

#!/bin/bash

# Directory where the seqs folder is stored
seqs_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/amr_serotyping/contigs"

# Output CSV file
output_file="SP_czidseqs.tsv"

# Write the header line
echo -e "sample_id\tseq_path" > $output_file

# Loop over the .fasta files in the seqs dir
find $seqs_dir -type f -name "*.fasta" | while read fasta_file; do
    # Extract the sample_id from the file name using regex
    sample_id=$(basename "$fasta_file" | sed -E 's/(_[0-9]+_contigs_nh)?\.fasta//')
    # Define the new file name
    new_fasta_file="${seqs_dir}/${sample_id}.fasta"
    # Rename the file
    mv "$fasta_file" "$new_fasta_file"
    # Write a line to the CSV file with the path to the new fasta file
    echo -e "$sample_id\t$new_fasta_file" >> $output_file
done

#### for reads

#!/bin/bash

# Directory where GPSCs are stored
GPSC_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/pangenome_alignment/GitHub/pipelines/genome_annotation/seqs/TB_reps/"

# Output text file
output_file="TB_reps.txt"

# Loop over the fasta files in GPSCs dir and write their paths to the output file
for fasta_file in $GPSC_dir*.fa; do
    echo "$fasta_file" >> $output_file
done



#### for contigs 

#!/bin/bash

# Directory where the seqs folder is stored
seqs_dir="/vast/palmer/scratch/turner/flg9/snakemake_workflows/1_pangenome/seqs/czid_seqs/*"

# Output text file
output_file="SP_czidseqs.txt"

# Loop over the contigs.fa files in the seqs dir and write their paths to the output file
find $seqs_dir -type f -name "contigs.fa" | while read fasta_file; do
    echo "$fasta_file" >> $output_file
done