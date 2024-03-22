{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [],
   "source": [
    "import re\n",
    "\n",
    "#This python script aims to remove any information that are not accession number in the VERSION field of the genbank file \n",
    "\n",
    "def process_genbank_file(input_file, output_file):\n",
    "    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:\n",
    "        for line in infile:\n",
    "            # Define a regular expression pattern to match the version field\n",
    "            pattern = re.compile(r'^\\s*ACCESSION\\s+(\\S+)')\n",
    "            match = pattern.match(line)\n",
    "\n",
    "            # Check if the line contains the version field\n",
    "            if match:\n",
    "                locus_field = match.group(1)\n",
    "                accession_number = locus_field.split('-')[-1]\n",
    "                # Replace the original version field with the accession number\n",
    "                modified_line = pattern.sub(r'ACCESSION   {}'.format(accession_number), line)\n",
    "                outfile.write(modified_line)\n",
    "            else:\n",
    "                # Keep other lines unchanged\n",
    "                outfile.write(line)\n",
    "\n",
    "# Example usage\n",
    "input_genbank_file = \"/Users/xanderlee/Desktop/Honours/Database/Archive_Database/TnCentral_Database/TnCentral_Genbank_Files/Final_TnCentral_Genbank_v1.gbff\"\n",
    "output_genbank_file = \"/Users/xanderlee/Desktop/Honours/Database/Archive_Database/TnCentral_Database/TnCentral_Genbank_Files/Final_TnCentral_Genbank_v2.gbff\"\n",
    "process_genbank_file(input_genbank_file, output_genbank_file)\n"
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
