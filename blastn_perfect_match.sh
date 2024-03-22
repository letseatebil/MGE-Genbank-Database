#!/bin/bash

# Define paths to files and databases
fasta_file="/home/xander/Desktop/Honours/Database/Database_V4/Original_concatenated_FASTA_database.fna"
blast_db="/home/xander/Desktop/Honours/Database/Database_V6/Blast_Database_V6"
inter_blast="/home/xander/Desktop/Honours/Database/Database_V6/interim_blast_results.txt"
output_blast='/home/xander/Desktop/Honours/Database/Database_V6/final_blast_results.txt'
inter_blast_temp='/home/xander/Desktop/Honours/Database/Database_V6/temp_blast.txt'
non_matching_accessions='/home/xander/Desktop/Honours/Database/Database_V6/non_matching_accession.txt'

# First part of the bash script is to perform BLAST on the FASTA database
echo "Initiating BLAST"
blastn_output=$(blastn -query "$fasta_file" -db "$blast_db" -evalue 1e-150 -outfmt '6 std slen' -perc_identity 93)
echo "BLASTn completed"

awk_output=$(echo "$blastn_output" | awk '{print $1, $2, $3, $4, $9, $10, $12, $4/$13;}')
echo "First awk command completed"

# Sort the blast hits by its bit score and coverage, both in a numerical reverse order
# Then make the output in a tab-delimited format where it is output into a text file
sort_output=$(echo "$awk_output" | sort -k7,7nr -k8,8nr | sed 's/  */\t/g' > "$inter_blast")

echo "Sorting completed"
echo "Proceeding to force BLAST hits with blast output file"

#Second part of the bash script is to force match accession numbers between query_id and subject_id
while read -r field1 field2 rest_of_line;
do

    # Extract the string after the hyphen in field1
    # Remove version number in field2
    accession_number=$(echo "$field1" | awk -F '-' '{print $NF}' | awk '{sub(/\..*$/, ""); print}')
    modified_field2=$(echo "$field2" | awk '{sub(/\..*$/, ""); print}')

    # Print progress
    echo "Accession number is: $accession_number"
    echo "Subject ID is: $modified_field2"

    if [ "$accession_number" = "$modified_field2" ]; then

        #Print matching hits and output into an interim temporary text file
        echo "$field1\t$field2\t$rest_of_line" >> "$inter_blast_temp"

        #Remove any duplicate hits of different start and end position and output to final blast output file
        cat $inter_blast_temp | awk '!seen[$1]++ {print $1, $2, $5, $6}' > $output_blast
        echo "Match found for $accession_number. Result added to $output_blast."
        echo "------------------------"
    else
        echo "Match not found for $accession_number, non-matching $accession_number output to text file for analysis"
        echo "$accession_number" >> "$non_matching_accessions"
    fi

done < "$inter_blast"

rm $inter_blast_temp

echo "Force filtering completed"
