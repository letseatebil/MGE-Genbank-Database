{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "from Bio import SeqIO\n",
    "from Bio.SeqRecord import SeqRecord\n",
    "\n",
    "TnCentral_Genbank = '/Users/xanderlee/Desktop/Honours/Database/Database_V6/TnCentral.gbff'\n",
    "TnCentral_Modified_Genbank = '/Users/xanderlee/Desktop/Honours/Database/Database_V6/TnCentral_v2.gbff'\n",
    "\n",
    "#This python script extracts the string from the 'keywords' key and replaces the string in the locus key\n",
    "extracted_records = []\n",
    "# Explicitly specify the encoding as 'latin1'\n",
    "with open(TnCentral_Genbank, 'r', encoding='latin1') as file:\n",
    "    for record in SeqIO.parse(file, \"genbank\"):\n",
    "        # Extract the keywords from the annotations dictionary\n",
    "        keywords = record.annotations.get('keywords', [])\n",
    "\n",
    "        # Use string found in keyword as the ID\n",
    "        locus_tag = \"_\".join(keywords)\n",
    "\n",
    "        # Remove element name in the locus field\n",
    "        locus_tag_parts = locus_tag.split('-')\n",
    "        modified_locus_tag = locus_tag_parts[-1]\n",
    "\n",
    "        extracted_record = SeqRecord(\n",
    "            seq=record.seq,\n",
    "            id=f'{modified_locus_tag}',\n",
    "            annotations=record.annotations,\n",
    "            features=record.features,\n",
    "        )\n",
    "        extracted_records.append(extracted_record)\n",
    "\n",
    "with open(TnCentral_Modified_Genbank, \"w\") as output_file:\n",
    "    SeqIO.write(extracted_records, output_file, \"genbank\")\n"
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
   "version": "3.11.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
