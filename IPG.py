#!/usr/bin/env python
# coding: utf-8

# In[41]:


# import pandas as pd
# import requests

# # Read Excel file
# df = pd.read_excel('./data/bca_individual_species 2.xlsx')

# # Extract unique protein accessions
# accession_values = df['Accession'].unique().tolist()

# len(accession_values)


# In[16]:


# from Bio import Entrez
# import pandas as pd

# # Set up Entrez
# Entrez.email = "hodane@utu.fi"
# Entrez.api_key = "..."

# def get_refseq_data(protein_accession):
#     try:
#         # Search for RefSeq entries using the protein accession code
#         handle = Entrez.esearch(db="protein", term=protein_accession + '[Accession]', retmode='xml')
#         records = Entrez.read(handle)

#         # Check if any records were found
#         if 'IdList' in records and len(records['IdList']) > 0:
#             refseq_id = records['IdList'][0]

#             # Retrieve RefSeq entry
#             handle = Entrez.efetch(db="protein", id=refseq_id, rettype="ipg", retmode="text")
#             refseq_data = handle.read()
#             return refseq_data
#         else:
#             print(f"No RefSeq entry found for accession {protein_accession}")
#             return None
#     except Exception as e:
#         print(f"Error fetching RefSeq data: {e}")
#         return None

# def parse_refseq_data(refseq_data):
#     if refseq_data:
#         # Decode the bytes data to a string
#         refseq_data_str = refseq_data.decode("utf-8")

#         # Extract IPG entries (if available)
#         refseq_entries = refseq_data_str.split("\n\n")
#         return refseq_entries
#     else:
#         return []

# # Example list of protein accession values
# accession_values = ['WP_046254879.1',
#  'WP_046255806.1',
#  'WP_046254117.1']

# # Initialize an empty list to store all entries
# all_entries = []

# for protein_accession in accession_values:
#     refseq_data = get_refseq_data(protein_accession)
#     if refseq_data:
#         refseq_entries = parse_refseq_data(refseq_data)
#         all_entries.extend(refseq_entries)

# # Save the entries to a CSV file without header
# with open("refseq_data2.csv", "w") as f:
#     # Write the header only once
#     f.write("Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\tProtein Name\tOrganism\tStrain\tAssembly\n")

#     # Write the entries
#     for entry in all_entries:
#         f.write(entry + "\n")

# print("CSV file saved successfully.")


# In[15]:


# from Bio import Entrez
# import pandas as pd

# # Set up Entrez
# Entrez.email = "hodane@utu.fi"
# Entrez.api_key = "..."

# def get_refseq_data(protein_accession):
#     try:
#         # Search for RefSeq entries using the protein accession code
#         handle = Entrez.esearch(db="protein", term=protein_accession + '[Accession]', retmode='xml')
#         records = Entrez.read(handle)

#         # Check if any records were found
#         if 'IdList' in records and len(records['IdList']) > 0:
#             refseq_id = records['IdList'][0]

#             # Retrieve RefSeq entry
#             handle = Entrez.efetch(db="protein", id=refseq_id, rettype="ipg", retmode="text")
#             refseq_data = handle.read()
#             return refseq_data
#         else:
#             print(f"No RefSeq entry found for accession {protein_accession}")
#             return None
#     except Exception as e:
#         print(f"Error fetching RefSeq data: {e}")
#         return None

# def parse_refseq_data(refseq_data):
#     if refseq_data:
#         # Decode the bytes data to a string
#         refseq_data_str = refseq_data.decode("utf-8")

#         # Extract IPG entries (if available)
#         refseq_entries = refseq_data_str.split("\n\n")
#         return refseq_entries
#     else:
#         return []

# # Example list of protein accession values
# accession_values = ['WP_046254879.1', 'WP_046255806.1', 'WP_046254117.1']

# # Initialize an empty list to store all entries
# all_entries = []

# for protein_accession in accession_values:
#     refseq_data = get_refseq_data(protein_accession)
#     if refseq_data:
#         refseq_entries = parse_refseq_data(refseq_data)
#         all_entries.extend(refseq_entries)

# # Create a DataFrame from the collected entries
# df2 = pd.DataFrame(all_entries)

# # Save the DataFrame to a CSV file without header
# df2.to_csv("refseq_data.csv", header=False, index=False)

# print("CSV file saved successfully.")


# ### final results

# In[17]:


import pandas as pd
import requests
df = pd.read_excel('./data/bca_individual_species 2.xlsx')
# Extract unique protein accessions
accession_values = df['Accession'].unique().tolist()

from Bio import Entrez
import pandas as pd

Entrez.email = "hodane@utu.fi"
Entrez.api_key = "..."

def get_refseq_data(protein_accession):
    try:
        # Search for RefSeq entries using the protein accession code
        handle = Entrez.esearch(db="protein", term=protein_accession + '[Accession]', retmode='xml')
        records = Entrez.read(handle)

        # Check if any records were found
        if 'IdList' in records and len(records['IdList']) > 0:
            refseq_id = records['IdList'][0]

            # Retrieve RefSeq entry
            handle = Entrez.efetch(db="protein", id=refseq_id, rettype="ipg", retmode="text")
            refseq_data = handle.read()
            return refseq_data
        else:
            print(f"No RefSeq entry found for accession {protein_accession}")
            return None
    except Exception as e:
        print(f"Error fetching RefSeq data: {e}")
        return None

def parse_refseq_data(refseq_data):
    if refseq_data:
        # Decode the bytes data to a string
        refseq_data_str = refseq_data.decode("utf-8")

        # Extract IPG entries (if available)
        refseq_entries = refseq_data_str.split("\n\n")
        return refseq_entries
    else:
        return []

# Example list of protein accession values
# accession_values = ['WP_046254879.1',#  'WP_046255806.1',#  'WP_046254117.1']

# Initialize an empty list to store all entries
all_entries = []

for protein_accession in accession_values:
    refseq_data = get_refseq_data(protein_accession)
    if refseq_data:
        refseq_entries = parse_refseq_data(refseq_data)
        all_entries.extend(refseq_entries)

# Save the entries to a CSV file without header
with open("refseq_data_final.csv", "w") as f:
    # Write the header only once
    f.write("Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\tProtein Name\tOrganism\tStrain\tAssembly\n")

    # Write the entries
    for entry in all_entries:
        f.write(entry + "\n")

print("CSV file saved successfully.")


# In[42]:


# import pandas as pd

# # Read the CSV file skipping the first and last rows
# df2 = pd.read_csv("refseq_data_final.csv", skiprows=2, header=None, delimiter="\t")

# # Rename the columns
# df2.columns = ["Id", "Source", "Nucleotide Accession", "Start", "Stop", "Strand", 
#               "Protein", "Protein Name", "Organism", "Strain", "Assembly"]

# # Display the DataFrame
# df2


# In[40]:


import pandas as pd

# Read the CSV file skipping the first and last rows
df_test = pd.read_csv("refseq_data_final.csv", skiprows=2, header=None, delimiter="\t")
# Rename the columns
df_test .columns = ["Id", "Source", "Nucleotide Accession", "Start", "Stop", "Strand", 
              "Protein", "Protein Name", "Organism", "Strain", "Assembly"]
# Exclude rows where "Id" column contains the string "Id"
df_test  = df_test [df_test ["Id"] != "Id"]
# Reset the index after filtering
df_test  = df_test .reset_index(drop=True)
# Display the DataFrame
df_test


# ## next task (29-4-2024):

# In[2]:


import pandas as pd
file_path = './data/filtered_ipgdata.xlsx'
df1 = pd.read_excel(file_path)
df1.head()


# In[3]:


file_path2 = './data/MB genomes data.xlsx'
df2 = pd.read_excel(file_path2)
df2.head()


# In[4]:


merged_df = pd.merge(df1, df2[['Assembly Accession']], how='inner', left_on='Assembly', right_on='Assembly Accession')
merged_df


# In[5]:


import pandas as pd
df = pd.read_excel("./data/bca_individual_species 2.xlsx")
# a dictionary to store species names and their corresponding accessions
species_accessions = {}
# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    species_name = row["species_individual"]
    accession = row["Accession"]
    ortholog = row["all_orth"]
    
    # Initialize an empty list for the species if not already present
    if species_name not in species_accessions:
        species_accessions[species_name] = {
            "CanA": [],
            "CanB": [],
            "SulP-CA": [],
            "BCA-5": [],
            "B-CARP": [],
            "Novel": []
        }
    
    # Append the accession to the corresponding ortholog list
    species_accessions[species_name][ortholog].append(accession)

# Create a DataFrame from the dictionary
ortholog_df  = pd.DataFrame(species_accessions)

# Transpose the DataFrame to have species names as rows and orthologs as columns
ortholog_df  = ortholog_df.T
ortholog_df


# In[13]:


import pandas as pd
# 'merged_df' and 'ortholog_df' loaded
# a dictionary to map protein accessions to ortholog names
protein_to_ortholog = {}
for species_name, row in ortholog_df.iterrows():
    for ortholog, accessions in row.items():
        for accession in accessions:
            protein_to_ortholog[accession] = ortholog

# Add an "Ortholog" column to merged_df
merged_df["Ortholog"] = merged_df["Protein"].map(protein_to_ortholog)
merged_df


# bases on your comment 2/5/2024:

# In[6]:


output_file_path = "./data/560rows×13.xlsx" 
merged_df.to_excel(output_file_path, index=False)


# In[12]:


# Create a dictionary to store species names and their corresponding accessions
species_accessions = {}

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    species_name = row["Organism"]
    accession = row["Protein"]  # Use the "Protein" column
    ortholog = row["Ortholog"]
    
    # Initialize an empty list for the species if not already present
    if species_name not in species_accessions:
        species_accessions[species_name] = {
            "CanA": [],
            "CanB": [],
            "SulP-CA": [],
            "BCA-5": [],
            "B-CARP": [],
            "Novel": []
        }
    
    # Check for NaN values in the ortholog column
    if not pd.isna(ortholog):
        # Append the accession to the corresponding ortholog list
        species_accessions[species_name][ortholog].append(accession)

# Print the updated species accessions
for species, accessions in species_accessions.items():
    print(f"Species: {species}")
    for ortholog, accession_list in accessions.items():
        print(f"{ortholog}: {accession_list}")
    print("\n")

