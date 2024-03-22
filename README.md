# MGE Genbank Database

This workflow can be used to generate a bacterial Mobile Genetic Element (MGE) database in the Genbank format. 

## Dependencies 

Requirements:

- Python 3.11.5 (biopython)
- [BLAST+] (https://www.ncbi.nlm.nih.gov/books/NBK279690/) (tested on version 2.10.1+)

## Running the workflow

### Basic usage

1. Extract accession number from your fasta database and run it through entrez_genbank_downloader.py
2. Concatenate the downloaded genbank files into a singular genbank file

### Skip this step if you do not require TnCentral Genbank annotations

___________________________________________________________________________________
3. Concatenate TnCentral annotation file Final_TnCentral_Genbank_v2.gbff with your singular genbank file, making sure that the TnCentral genbank annotation file is the first argument
___________________________________________________________________________________

3. Convert your concatenate genbank file back to FASTA format using the Genbank_Fasta_Converter.py script
4. Remove any duplicates in the FASTA file
5. Run the newly converted FASTA file into ```makeblastdb``` to generate your local BLAST database
6. BLAST your original fasta database against the local database using the blastn_perfect_match.sh BASH script
7. Parse the BLAST output file into the Genbank_Database_Generation.py script to generate your final curated Genbank database

