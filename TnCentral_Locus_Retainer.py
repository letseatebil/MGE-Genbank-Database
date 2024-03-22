import re

#This python script aims to remove any information that are not accession number in the VERSION field of the genbank file 

def process_genbank_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Define a regular expression pattern to match the version field
            pattern = re.compile(r'^\s*ACCESSION\s+(\S+)')
            match = pattern.match(line)

            # Check if the line contains the version field
            if match:
                locus_field = match.group(1)
                accession_number = locus_field.split('-')[-1]
                # Replace the original version field with the accession number
                modified_line = pattern.sub(r'ACCESSION   {}'.format(accession_number), line)
                outfile.write(modified_line)
            else:
                # Keep other lines unchanged
                outfile.write(line)

# Example usage
input_genbank_file = "/Users/xanderlee/Desktop/Honours/Database/Archive_Database/TnCentral_Database/TnCentral_Genbank_Files/Final_TnCentral_Genbank_v1.gbff"
output_genbank_file = "/Users/xanderlee/Desktop/Honours/Database/Archive_Database/TnCentral_Database/TnCentral_Genbank_Files/Final_TnCentral_Genbank_v2.gbff"
process_genbank_file(input_genbank_file, output_genbank_file)
