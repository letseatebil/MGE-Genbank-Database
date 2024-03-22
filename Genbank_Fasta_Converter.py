from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Define the input GenBank file with annotations
multigenbank_file = "/home/xander/Desktop/Honours/Database/Database_V6/Database_V6.gbff"

# Open the input GenBank file and create an output FASTA file for writing
with open(multigenbank_file) as input_handle, open("/home/xander/Desktop/Honours/Database/Database_V6/Converted_Database_V6.fna", "w") as output_handle, open("/home/xander/Desktop/Honours/Database/Database_V6/Skipped_Records.txt", "w") as skipped_handle:

    # Parse sequences from the input GenBank file
    sequences = SeqIO.parse(input_handle, "genbank")
    # Counter for successfully converted records
    count = 0
    # Counter for skipped records
    skipped = 0  
    for record in sequences:
        try:
            # Try to access the sequence content
            sequence_length = len(record.seq)
            if sequence_length > 0:
                # Write the record to the output FASTA file
                SeqIO.write(record, output_handle, "fasta")
                count += 1
        except Exception as e:
            # Write details of skipped record to the skipped records file
            skipped_handle.write(f"Skipped record due to error: {e}\n")
            skipped_handle.write(f"Record ID: {record.id}\n")
            skipped_handle.write(f"Record Description: {record.description}\n\n")
            skipped += 1

# Print the total number of skipped and converted records
print(f"Skipped {skipped} records. Details written to Skipped_Records.txt")
print("Converted %i records" % count)
