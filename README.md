# mondo2omop
This repository provides an example script for combining the Mondo Disease Ontology with the OMOP Concept and Concept Relationship tables, utilizing the precise mappings from Mondo to other terminologies, to produce a Mondo-to-OMOP mapping table that can be utilized to connect Mondo to an OMOP instance and generate Mondo patient cohorts.

# Prerequisites
You will need to have previously downloaded and prepared, or otherwise made available, the OMOP CONCEPT and CONCEPT_RELATIONSHIP tables, or made copies available from your local OMOP instance, and place copies of those tables in data/omop. If you don't yet have these OMOP tables available, they can be downloaded from https://athena.ohdsi.org/vocabulary/list. Be sure that the ICD10CM, SNOMED, and MeSH vocabularies are included in the download in order to utilize those mappings from Mondo.

Within the mondo_to_omop.py script, you will need to set the 'mondo_version' variable to indicate the release date of your desired Mondo version. Mondo KGX release versions can be found at https://kg-hub.berkeleybop.io/kg-obo/mondo/.

# Citation
TBD

# References
The Mondo-to-OMOP script makes use of NetworkX for graph-related operations:
Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008