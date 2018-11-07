#!/usr/bin/env bash

evalue="0.01"
num_iterations="3"

db_file=$1
input_fasta_dir=$2
output_pssm_dir=$3
output_blast_dir=$4

count=0
total_count=1
for file in "$input_fasta_dir"*; do
    filename=$(basename -- "$file")
    echo "$total_count process $file"
    psiblast -query "$file" -db "$db_file" -evalue "$evalue" -num_iterations "$num_iterations" -out_ascii_pssm "$output_pssm_dir${filename%.*}".pssm -out "$output_blast_dir${filename%.*}".alns.blast &
    ((total_count++))
    ((count++))
    if [ $count -gt 5 ]
    then
        echo batch
        count=0
        wait
    fi
done