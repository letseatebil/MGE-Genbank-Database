from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

TnCentral_Genbank = '/Users/xanderlee/Desktop/Honours/Database/Database_V6/TnCentral.gbff'
TnCentral_Modified_Genbank = '/Users/xanderlee/Desktop/Honours/Database/Database_V6/TnCentral_v2.gbff'

#This python script extracts the string from the 'keywords' key and replaces the string in the locus key
extracted_records = []
# Explicitly specify the encoding as 'latin1'
with open(TnCentral_Genbank, 'r', encoding='latin1') as file:
    for record in SeqIO.parse(file, "genbank"):
        # Extract the keywords from the annotations dictionary
        keywords = record.annotations.get('keywords', [])

        # Use string found in keyword as the ID
        locus_tag = "_".join(keywords)

        # Remove element name in the locus field
        locus_tag_parts = locus_tag.split('-')
        modified_locus_tag = locus_tag_parts[-1]

        extracted_record = SeqRecord(
            seq=record.seq,
            id=f'{modified_locus_tag}',
            annotations=record.annotations,
            features=record.features,
        )
        extracted_records.append(extracted_record)

with open(TnCentral_Modified_Genbank, "w") as output_file:
    SeqIO.write(extracted_records, output_file, "genbank")
