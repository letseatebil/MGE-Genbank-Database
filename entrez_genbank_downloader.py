from Bio import Entrez

# Read the file containing accession numbers
with open("/Users/xanderlee/Desktop/Honours/Database/Database_V3/unique_accessions_part_2.txt", "r") as file:
    accession_numbers = [line.strip() for line in file.readlines()]

# Entrez requires an email address for accessing their services
Entrez.email = "a1809437@student.adelaide.edu.au"

# Download GenBank annotations for each accession number
for acc in accession_numbers:
    count = 0
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gbwithparts", retmode="text")
        data = handle.read()
        with open(f"{acc}.gb", "w") as output_file:
            output_file.write(data)
        handle.close()
        print(f"Downloaded {acc}.gb")
    except Exception as e:
        print(f"Failed to download {acc}: {e}")
        count+=1

print(f"Skipped {count} entries")
