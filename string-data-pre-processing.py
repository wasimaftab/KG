#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:45:32 2024

@author: wasim
"""
import pandas as pd
import pdb
from neo4j import GraphDatabase

## Load STRING PPIs
ppi = pd.read_csv("7227.protein.physical.links.v12.0.txt", sep = "\t")

ppi['protein1'] = ppi['protein1'].str.replace('7227.', '', regex=False)
ppi['protein2'] = ppi['protein2'].str.replace('7227.', '', regex=False)

## Load FBgn to FBpp mapping form Flybase       
flybase_map = pd.read_csv("fbgn_fbtr_fbpp_fb_2024_02.tsv", sep="\t")

## Now load the Uniprot to FBgn id converion table for chromatin proteins
uniprot_fbgn = pd.read_csv("Chomatin_Uniprot_FBgn_idmapping_2024_05_08.tsv", sep = "\t")
 
## Subset fbgn_uniprot to only contain FBgn uniprot mapping for chromatin proteins
filtered_map = flybase_map[flybase_map['FlyBase_FBgn'].isin(uniprot_fbgn.To.to_list())].copy()

protein = list(set(filtered_map.FlyBase_FBpp.to_list()))
matches = ppi[ppi['protein1'].isin(protein) & ppi['protein2'].isin(protein)]

# Create a mapping from FBgn to UniProt IDs using a dictionary
fbgn_to_uniprot = dict(zip(uniprot_fbgn['To'], uniprot_fbgn['From']))

# Create a mapping from FBpp to UniProt by using the FBgn to UniProt map on the filtered map
filtered_map.loc[:, 'UniProt'] = filtered_map['FlyBase_FBgn'].map(fbgn_to_uniprot)

# Create dictionary for quick lookup
fbpp_to_uniprot = pd.Series(filtered_map['UniProt'].values, index=filtered_map['FlyBase_FBpp']).to_dict()

# pdb.set_trace()
# Replace FBpp IDs in 'matches' with UniProt IDs
matches.loc[:, 'protein1'] = matches['protein1'].map(fbpp_to_uniprot).fillna(matches['protein1'])
matches.loc[:, 'protein2'] = matches['protein2'].map(fbpp_to_uniprot).fillna(matches['protein2'])

# Print the updated 'matches' DataFrame
print(matches.head())

# pdb.set_trace()

# ## To neo4j
# from neo4j import GraphDatabase

# # Connection to the Neo4j database
# uri = "bolt://localhost:7687"  # Default URI
# user = "neo4j"  # Default user
# password = ""  # Replace with your actual password

# driver = GraphDatabase.driver(uri, auth=(user, password))

# def add_protein_interaction(tx, protein1, protein2):
#     tx.run("MERGE (p1:Protein {id: $protein1}) "
#             "MERGE (p2:Protein {id: $protein2}) "
#             "MERGE (p1)-[:INTERACTS_WITH]->(p2)", protein1=protein1, protein2=protein2)

# # Assuming 'df' is your DataFrame
# with driver.session() as session:
#     for index, row in matches.iterrows():
#         session.write_transaction(add_protein_interaction, row['protein1'], row['protein2'])

# driver.close()
# pdb.set_trace()

## Annotate further
quickGO = pd.read_excel('QuickGO-annotations-Chromatin-Dmel.xlsx')
df_descriptions = pd.DataFrame({
    'GENE PRODUCT ID': quickGO["GENE PRODUCT ID"].to_list(),
    'GO NAME': quickGO["GO NAME"].to_list()
})

# Group by 'GENE PRODUCT ID' and join the 'GO NAME' entries
df_descriptions = df_descriptions.groupby('GENE PRODUCT ID')['GO NAME'].apply('; '.join).reset_index()

# Remove duplicate entries within the 'GO NAME' field
df_descriptions['GO NAME'] = df_descriptions['GO NAME'].apply(lambda x: '; '.join(sorted(set(x.split('; ')), key=x.split('; ').index)))

print(df_descriptions)
# pdb.set_trace()

## Update the Neo4j Database
# Connection to the Neo4j database
uri = "bolt://localhost:7687"  # Default URI
user = "neo4j"  # Default user
password = ""  # Replace with your actual password

driver = GraphDatabase.driver(uri, auth=(user, password))

def add_description(tx, protein_id, description):
    tx.run("MERGE (p:Protein {id: $protein_id}) "
           "SET p.description = $description", protein_id=protein_id, description=description)

with driver.session() as session:
    for index, row in df_descriptions.iterrows():
        session.write_transaction(add_description, row['GENE PRODUCT ID'], row['GO NAME'])

driver.close()

