{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Skipped 18 records. Details written to Skipped_Records.txt\n",
      "Converted 3978 records\n"
     ]
    }
   ],
   "source": [
    "from Bio import SeqIO\n",
    "from Bio.SeqRecord import SeqRecord\n",
    "\n",
    "# Define the input GenBank file with annotations\n",
    "multigenbank_file = \"/home/xander/Desktop/Honours/Database/Database_V6/Database_V6.gbff\"\n",
    "\n",
    "# Open the input GenBank file and create an output FASTA file for writing\n",
    "with open(multigenbank_file) as input_handle, open(\"/home/xander/Desktop/Honours/Database/Database_V6/Converted_Database_V6.fna\", \"w\") as output_handle, open(\"/home/xander/Desktop/Honours/Database/Database_V6/Skipped_Records.txt\", \"w\") as skipped_handle:\n",
    "\n",
    "    # Parse sequences from the input GenBank file\n",
    "    sequences = SeqIO.parse(input_handle, \"genbank\")\n",
    "    # Counter for successfully converted records\n",
    "    count = 0\n",
    "    # Counter for skipped records\n",
    "    skipped = 0  \n",
    "    for record in sequences:\n",
    "        try:\n",
    "            # Try to access the sequence content\n",
    "            sequence_length = len(record.seq)\n",
    "            if sequence_length > 0:\n",
    "                # Write the record to the output FASTA file\n",
    "                SeqIO.write(record, output_handle, \"fasta\")\n",
    "                count += 1\n",
    "        except Exception as e:\n",
    "            # Write details of skipped record to the skipped records file\n",
    "            skipped_handle.write(f\"Skipped record due to error: {e}\\n\")\n",
    "            skipped_handle.write(f\"Record ID: {record.id}\\n\")\n",
    "            skipped_handle.write(f\"Record Description: {record.description}\\n\\n\")\n",
    "            skipped += 1\n",
    "\n",
    "# Print the total number of skipped and converted records\n",
    "print(f\"Skipped {skipped} records. Details written to Skipped_Records.txt\")\n",
    "print(\"Converted %i records\" % count)\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "bioinf",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
