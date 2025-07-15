import os
import numpy as np
import pandas as pd
import requests
import tarfile
import networkx as nx

# Set up display preferences for viewing pandas dataframes.
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

# Set up dataset directories
paths = ['data/mondo',
         'data/mondo2omop',
         'data/omop']

for path in paths:
    if not os.path.exists(path):
        os.makedirs(path)


# *** Indicate the release date of the version of Mondo you wish to download from ***
# https://kg-hub.berkeleybop.io/kg-obo/mondo/ in the format YYYY-MM-DD:
mondo_version = '2025-04-01'

# Get and extract the Mondo KGX files
mondo_url = 'https://kg-hub.berkeleybop.io/kg-obo/mondo/' + mondo_version + '/mondo_kgx_tsv.tar.gz'
filename = 'data/mondo/mondo_kgx_tsv.tar.gz'
with open(filename, "wb") as f:
    r = requests.get(mondo_url)
    f.write(r.content)
file = tarfile.open(filename)
file.extractall(path='data/mondo/')
file.close()

# Load the Mondo edges and nodes files
mondo_edges_df = pd.read_csv("data/mondo/mondo_kgx_tsv_edges.tsv", sep='\t', header=0, low_memory=False)
mondo_nodes_df = pd.read_csv("data/mondo/mondo_kgx_tsv_nodes.tsv", sep='\t', header=0, low_memory=False)


# *** Build the Mondo graph from the Mondo edges and nodes files. ***
edges_df = mondo_edges_df.copy()
edges_df = edges_df[['subject', 'object', 'predicate']]

nodes_df = mondo_nodes_df.copy()
nodes_df = nodes_df[['id', 'category', 'name']]
mondo_graph = pd.merge(edges_df, nodes_df, left_on='subject', right_on='id', how='left')
mondo_graph = pd.merge(mondo_graph, nodes_df, left_on='object', right_on='id', how='left', suffixes=('', '_obj'))

# Filter the Mondo graph to needed categories and predicates, remove obsolete nodes, and retain just the subject and
# object columns.
mondo_graph = mondo_graph[(mondo_graph.category == 'biolink:Disease') & (mondo_graph.category_obj == 'biolink:Disease') & (mondo_graph.predicate == 'biolink:subclass_of')]
mondo_graph = mondo_graph[~mondo_graph['name'].str.contains('obsolete')]
mondo_graph = mondo_graph[~mondo_graph['name_obj'].str.contains('obsolete')]
mondo_graph = mondo_graph[['subject', 'object']]

# Build the mondo graph with NetworkX so that it can be used to identify term descendants.
G = nx.DiGraph()
for index, row in mondo_graph.iterrows():
    G.add_edge(row['object'], row['subject'])


# Filter the nodes dataframe to just the disease nodes, explode the same_as column to multiple rows.
nodes = mondo_nodes_df.copy()
nodes = nodes[nodes['category'] == 'biolink:Disease']

# Remove any obsolete terms
nodes.loc[:, ('name')] = nodes.loc[:, ('name')].astype(str)
nodes = nodes[~nodes['name'].str.contains(r'obsolete(?!$)')]

# Remove any terms not a descendant of human disease, or that are descendants of disease susceptibility,
# disease characteristic, or injury. Filter the nodes dataframe to only the descendants of human disease (MONDO:0700096)
human_disease_descendants = list(nx.descendants(G, 'MONDO:0700096'))
nodes = nodes[nodes['id'].isin(human_disease_descendants)]

# Filter out the descendants of disease susceptibility (MONDO:0042489)
disease_susceptibility_descendants = list(nx.descendants(G, 'MONDO:0042489'))
nodes = nodes[~nodes['id'].isin(disease_susceptibility_descendants)]

# Filter out the descendants of disease characteristic (MONDO:0021125)
disease_characteristic_descendants = list(nx.descendants(G, 'MONDO:0021125'))
nodes = nodes[~nodes['id'].isin(disease_characteristic_descendants)]

# Filter out the descendants of injury (MONDO:0021178)
injury_descendants = list(nx.descendants(G, 'MONDO:0021178'))
nodes = nodes[~nodes['id'].isin(injury_descendants)]


# *** Parse the 'same as' column ***
same_as_df = nodes.copy()
same_as_df['same_as_array'] = same_as_df['same_as'].str.split('|')
same_as_df = same_as_df.explode('same_as_array')

# Retain just the columns of interest.
same_as_df = same_as_df[['id', 'name', 'description', 'category', 'same_as_array']]
# Filter the same_as_array to just the identifiers of interest (ICD10, SNOMED, MeSH)
substrings = ['http://identifiers.org/snomedct/', 'http://identifiers.org/mesh/',
              'http://purl.bioontology.org/ontology/ICD10CM/']
same_as_df = same_as_df[same_as_df['same_as_array'].str.contains('|'.join(substrings), na=False)]

# Add a 'concept_code' column with the prefixes from the vocabulary codes removed.
same_as_df['concept_code'] = (same_as_df['same_as_array'].replace('http://purl.bioontology.org/ontology/ICD10CM/', '', regex=True)
                                                         .replace('http://identifiers.org/snomedct/', '', regex=True)
                                                         .replace('http://identifiers.org/mesh/', '', regex=True))

# Set up a list of conditions and values to be used to add a 'vocabulary' column, with the vocabulary ID corresponding
# to the indicated vocabulary.
conditions = [
    (same_as_df['same_as_array'].str.contains(r'http://identifiers.org/snomedct/')),
    (same_as_df['same_as_array'].str.contains(r'http://identifiers.org/mesh/')),
    (same_as_df['same_as_array'].str.contains(r'http://purl.bioontology.org/ontology/ICD10CM/'))
]
values = ['SNOMED', 'MeSH', 'ICD10CM']
same_as_df['vocabulary'] = np.select(conditions, values, default='Other')
same_as_df = same_as_df.rename(columns={'same_as_array': 'same_as'})


# *** Parse the 'subsets' column ***
subsets_df = nodes.copy()
subsets_exploded_df = nodes.copy()
subsets_exploded_df['subsets_array'] = subsets_exploded_df['subsets'].str.split('|')
subsets_exploded_df = subsets_exploded_df.explode('subsets_array')

# Add a binary column for each 'rare' designation.
rare_subset = subsets_exploded_df.copy()
rare_subset = rare_subset[rare_subset['subsets_array'] == 'rare']
rare_subset = rare_subset[['id']]
rare_subset = rare_subset.rename(columns={'id': 'rare_id'})
subsets_df = pd.merge(subsets_df, rare_subset, left_on='id', right_on='rare_id', how='left')
subsets_df['rare'] = np.where(subsets_df['rare_id'].notnull(), 1, 0)
subsets_df.drop('rare_id', axis=1, inplace=True)

gard_rare = subsets_exploded_df.copy()
gard_rare = gard_rare[gard_rare['subsets_array'] == 'gard_rare']
gard_rare = gard_rare[['id']]
gard_rare = gard_rare.rename(columns={'id': 'rare_id'})
subsets_df = pd.merge(subsets_df, gard_rare, left_on='id', right_on='rare_id', how='left')
subsets_df['gard_rare'] = np.where(subsets_df['rare_id'].notnull(), 1, 0)
subsets_df.drop('rare_id', axis=1, inplace=True)

nord_rare = subsets_exploded_df.copy()
nord_rare = nord_rare[nord_rare['subsets_array'] == 'nord_rare']
nord_rare = nord_rare[['id']]
nord_rare = nord_rare.rename(columns={'id': 'rare_id'})
subsets_df = pd.merge(subsets_df, nord_rare, left_on='id', right_on='rare_id', how='left')
subsets_df['nord_rare'] = np.where(subsets_df['rare_id'].notnull(), 1, 0)
subsets_df.drop('rare_id', axis=1, inplace=True)

orphanet_rare = subsets_exploded_df.copy()
orphanet_rare = orphanet_rare[orphanet_rare['subsets_array'] == 'orphanet_rare']
orphanet_rare = orphanet_rare[['id']]
orphanet_rare = orphanet_rare.rename(columns={'id': 'rare_id'})
subsets_df = pd.merge(subsets_df, orphanet_rare, left_on='id', right_on='rare_id', how='left')
subsets_df['orphanet_rare'] = np.where(subsets_df['rare_id'].notnull(), 1, 0)
subsets_df.drop('rare_id', axis=1, inplace=True)

inferred_rare = subsets_exploded_df.copy()
inferred_rare = inferred_rare[inferred_rare['subsets_array'] == 'inferred_rare']
inferred_rare = inferred_rare[['id']]
inferred_rare = inferred_rare.rename(columns={'id': 'rare_id'})
subsets_df = pd.merge(subsets_df, inferred_rare, left_on='id', right_on='rare_id', how='left')
subsets_df['inferred_rare'] = np.where(subsets_df['rare_id'].notnull(), 1, 0)
subsets_df.drop('rare_id', axis=1, inplace=True)

mondo_rare = subsets_exploded_df.copy()
mondo_rare = mondo_rare[mondo_rare['subsets_array'] == 'mondo_rare']
mondo_rare = mondo_rare[['id']]
mondo_rare = mondo_rare.rename(columns={'id': 'rare_id'})
subsets_df = pd.merge(subsets_df, mondo_rare, left_on='id', right_on='rare_id', how='left')
subsets_df['mondo_rare'] = np.where(subsets_df['rare_id'].notnull(), 1, 0)
subsets_df.drop('rare_id', axis=1, inplace=True)

# Select the final subsets columns
subsets_df = subsets_df[['id', 'rare', 'gard_rare', 'nord_rare', 'orphanet_rare', 'inferred_rare', 'mondo_rare']]

# *** Map Mondo to OMOP ***
# Load the OMOP concept table.
concept_df = pd.read_csv("data/omop/CONCEPT.csv", sep='\t', header=0, low_memory=False)

# Load the OMOP concept_relationship table
concept_relationship_df = pd.read_csv("data/omop/CONCEPT_RELATIONSHIP.csv", sep='\t', header=0, low_memory=False)

# Merge the same_as dataframe with the OMOP concept and concept_relationship tables.
mondo_to_omop_df = same_as_df.copy()
mondo_to_omop_df = mondo_to_omop_df.rename(columns={'vocabulary': 'vocabulary_id'})
mondo_to_omop_df = pd.merge(mondo_to_omop_df, concept_df,
                            on=['concept_code', 'vocabulary_id'],
                            how='inner')
maps_to_df = concept_relationship_df.copy()
maps_to_df = maps_to_df[maps_to_df['relationship_id'] == 'Maps to']
maps_to_df = maps_to_df[['concept_id_1', 'concept_id_2', 'relationship_id']]

mondo_to_omop_df = pd.merge(mondo_to_omop_df, maps_to_df, left_on='concept_id', right_on='concept_id_1', how='inner')

# Prep a copy of the concept table for mapping to standard concepts.
standard_concept_df = concept_df.copy()
# Pre-filter to OMOP 'standard' concepts (standard_concept == 'S').
standard_concept_df = standard_concept_df[standard_concept_df['standard_concept'] == 'S']
standard_concept_df = standard_concept_df.rename(columns={'concept_id': 'standard_concept_id',
                                                          'concept_name': 'standard_concept_name',
                                                          'vocabulary_id': 'standard_vocabulary_id',
                                                          'domain_id': 'standard_domain_id',
                                                          'concept_code': 'standard_concept_code'})
standard_concept_df = standard_concept_df[['standard_concept_id', 'standard_concept_name', 'standard_vocabulary_id', 'standard_domain_id', 'standard_concept_code']]
standard_concept_df = standard_concept_df[standard_concept_df['standard_domain_id'] == 'Condition']

# Merge in mappings to standard concepts.
mondo_to_omop_df = pd.merge(mondo_to_omop_df, standard_concept_df, left_on='concept_id_2', right_on='standard_concept_id', how='inner')
mondo_to_omop_df.drop('concept_id_1', axis=1, inplace=True)
mondo_to_omop_df.drop('concept_id_2', axis=1, inplace=True)

# Join in the rare disease subset designation columns
mondo_to_omop_df = pd.merge(mondo_to_omop_df, subsets_df, left_on='id', right_on='id', how='left')

# Save the final Mondo-to-OMOP dataframe to a tab-delimited file.
mondo_to_omop_df.to_csv("data/mondo2omop/MONDO2OMOP.tsv", sep='\t', index=False)
