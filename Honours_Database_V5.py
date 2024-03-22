from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import warnings
import time
from datetime import date
import re

# Ignore warnings that are appearing in the command line
warnings.filterwarnings("ignore")

# Get current script directory
script_dir = os.path.dirname(os.path.realpath(__file__))
blast_output_file = "/home/xander/Desktop/Honours/Database/Database_V6/93_percidentity_final_blast_results.txt"

# Assign genbank database to variable
genbank_database = "/home/xander/Desktop/Honours/Database/Database_V6/Database_V6.gbff"
genbank_db_path = os.path.join(script_dir, genbank_database)

# Get today's date so it will output the date the genbank database was generated
today = date.today()
date4 = today.strftime('%d-%b-%Y').upper()

script_output = "/home/xander/Desktop/Honours/Database/Database_V6/Final_Genbank_Database_V6.gbff"

def extract_sequences():
    # Create a dictionary to store subject positions and query_id for each accession
    subject_positions = {}

    with open(blast_output_file, "r") as blast_file:

        # Read Blast output file and store subject start, stop positions, and query_id for each accession
        for line in blast_file:

            #Strips any leading or trailing whitespaces from the current line
            #Then splits it into a list of fields based on spaces
            fields = line.strip().split(" ")
            
            #Unpacks the first four elements of the fields list into separate variables
            query_id, subject_id, subject_start, subject_end = fields[:4]

            #Checks if the subject start position is greater than the subject end position
            #If true, swap the values to ensure proper order
            if int(subject_start) > int(subject_end):
                subject_start, subject_end = subject_end, subject_start

            #Checks if the subject_id is not already a key in the subject_positions dictionary
            if subject_id not in subject_positions:
                
                #If the subject_id is not in the dictionary, add a new entry with an empty list for 'query_ids' and 'positions'
                subject_positions[subject_id] = {'query_ids': [], 'positions': []}

            #Appends the query_id to the list associated with the 'query_ids' key in the dictionary
            subject_positions[subject_id]['query_ids'].append(query_id)

            #Appends a tuple (int(subject_start), int(subject_end)) to the list associated with the 'positions' key in the dictionary
            subject_positions[subject_id]['positions'].append((int(subject_start), int(subject_end)))

    extracted_records = {}

    for record in SeqIO.parse(genbank_db_path, "genbank"):
        accession = record.id

        # Assigning fixed values to keys in the SeqRecord annotation dictionary
        record.annotations['topology'] = "linear"
        record.annotations['references'] = "."
        record.annotations['data_file_division'] = "BCT"
        record.annotations['molecule_type'] = "DNA"

        print(f"Accession Number is: {record.id}")
        # Check if the accession is in subject_positions dictionary
        if accession in subject_positions:

            print(f"{accession} matches subject position")

            for i, (start, end) in enumerate(subject_positions[accession]['positions']):

                query_id = subject_positions[accession]['query_ids'][i]
                print(f"Query ID is {query_id}")

                # Extract sequence based on start and end positions
                sequence = record.seq[min(start, end) - 1: max(start, end)]

                # Use a tuple of (accession, query_id) as the key to ensure uniqueness
                key = (accession, query_id)
                print(f"Accession & Query ID are: {key}")
                if key not in extracted_records:
                    extracted_records[key] = SeqRecord(
                        seq=sequence,
                        id=f"{accession} REGION: {min(start, end) - 1}..{max(start, end)}",
                        name=query_id,
                        description=f"Annotations filtered from MEFinder, TnCentral & MGEDB\t",
                        annotations={
                            'topology': record.annotations['topology'],
                            'molecule_type': record.annotations['molecule_type'],
                            'date': date4,
                            'data_file_division': record.annotations['data_file_division'],
                            'references': record.annotations['references'],
                            'organism': record.annotations['organism'],
                        },
                    )

                    print("Generating new Genbank entries for", record.name)
                    print("Genbank record generation complete")
                    print("Generating features for", record.name)

                    # Add custom annotation based on specific query_id
                    priority_element = ["Tn3-V00613", "Tn21-AF071413", "Tn1721-X61367", "Tn1331-NC_003486", "IS26-X00011"]
                    if query_id in priority_element:
                        extracted_records[key].description += "*Priority Element*"

                # Filter features within the specified start and end positions
                for feature in record.features:
                    
                    #Declaring variables
                    #----------------------------------------------------------------------------------------------------------

                    #Adjust to 0-based indexing for python
                    feature_start = feature.location.start + 1
                    feature_end = feature.location.end

                    #Extract element name from locus field
                    element_name = re.split("-", query_id)

                    #Check for presence of /note or /label keys in qualifiers dictionary 
                    note_value = feature.qualifiers.get('note', [''])[0] if 'note' in feature.qualifiers else ''
                    label_value = feature.qualifiers.get('label', [''])[0] if 'label' in feature.qualifiers else ''
                    
                    #Filtering feature type for internal elements
                    #desired_feature_types = ['mobile_element', 'repeat_region', 'CDS', 'misc_feature']

                    #----------------------------------------------------------------------------------------------------------
                    #Checking for internal elements and annotating it as such
                    #Identify criteria for internal elements (case insensitive)
                    '''is_internal_element = element_name[0].casefold() not in note_value or element_name[0].casefold() not in label_value
                    
                    if is_internal_element and feature.type in desired_feature_types:
                            
                        internal_element_tag = 'Internal element = yes'

                        # Check if "Internal element = yes" is already present in note_value 
                        if internal_element_tag not in note_value:

                            # Add a custom qualifier to indicate an internal element
                            #print(element_name[0].casefold())
                            feature.qualifiers['note'] = [f'{note_value}; {internal_element_tag}']
                        else:

                            # If it's not an internal element, remove any existing "Internal element = yes" from the note
                            if 'Internal element = yes' in note_value:
                                feature.qualifiers['note'] = [note_value.replace('Internal element = yes', '').strip()]
                    '''
                    #Filtering annotations based on start and end position
                    #----------------------------------------------------------------------------------------------------------
                    if (
                        feature_start >= min(start, end)
                        and feature_end <= max(start, end)
                    ):
                        #Check for keywords specific to inverted repeats on both bounderies of element and append 'transposition feature' in /note qualifier
                        if (
                            'note' in feature.qualifiers
                            or 'label' in feature.qualifiers
                        ):

                            inverted_repeat_keywords = ['inverted repeat', 'IRR', 'IRL', 'IR', 'inverted terminal repeat', 'terminal IR of', 'IRt']

                            #----------------------------------------------------------------------------------------------------------
                            #Check if any of the keywords are found in /note qualifier
                            if any(keyword in note_value or keyword in label_value for keyword in inverted_repeat_keywords):

                                #Check if the keywords are found in boundary inverted repeats
                                #If yes,append annotation to these boundary inverted repeats
                                #Otherwise, replace 'Transposition feature = no' with blankspace
                                if (
                                    feature_start == start
                                    or feature_end == end
                                ):
                                    if 'Transposition feature' not in note_value:
                                        feature.qualifiers['note'] = [f'Transposition feature = yes; {note_value}']
                                else:
                                    feature.qualifiers['note'] = [note_value.replace('Transposition feature = no', '')]

                        #----------------------------------------------------------------------------------------------------------
                        # Check for transposase and resolvase in the various annotations
                        keywords = ['tnpA', 'TnpA', 'tnp', 'Tnp', 'tnpR', 'TnpR', 'Transposase', 'transposase', 'Resolvase', 'resolvase']

                        is_transposase_or_resolvase = any(
                            keywords in feature.qualifiers.get(qualifier, '')
                            for feature in [feature]
                            for qualifier in ['gene', 'function', 'note', 'product']
                            if qualifier in feature.qualifiers
                        )

                        if is_transposase_or_resolvase:
                            print(f'Element name is: {element_name[0]}')
                            print('Main Transposase or Resolvase found!')

                            # Check whether /note qualifier is present in the annotation, if it is absent, add new /note key and annotation to qualifiers dictionary
                            # Otherwise, just append new annotation to current /note key
                            if 'note' not in feature.qualifiers:
                                feature.qualifiers |= {'note': 'Transposition feature = yes'}
                            else:
                                if 'Transposition feature' not in feature.qualifiers['note'][0]:
                                    feature.qualifiers['note'] = [f'Transposition feature = yes; {note_value}']
                        
                        #----------------------------------------------------------------------------------------------------------
                        # Specify start and end location of the feature
                        feature_location = FeatureLocation(
                            feature.location.start,
                            feature.location.end,
                            strand=feature.location.strand,
                        )

                        # Represent a Sequence Feature holding info about a part of a sequence
                        extracted_feature = SeqFeature(
                            location=feature_location,
                            type=feature.type,
                            qualifiers=feature.qualifiers,
                        )
                        extracted_records[key].features.append(extracted_feature)

                #print(feature.qualifiers)
                print("Feature generation complete")
                print("Appending records")
                print("---------------------------------")
        else:
            print(f"Record for {record.id} is invalid or not found")

    # Write extracted sequences with features to output file in GenBank format
    with open(script_output, "w") as output_file:
        SeqIO.write([record for record in extracted_records.values()], output_file, "genbank")

def main():
    # Log down time taken to blast and generate genbank file from the query file
    start_time = time.time()

    extract_sequences()

    end_time = time.time()
    elapsed_time = end_time - start_time

    # Print the time taken to run the script
    print(f"Time taken: {elapsed_time} seconds")

    os.system(f"espeak 'Script Executed'")

if __name__ == "__main__":
    main()
