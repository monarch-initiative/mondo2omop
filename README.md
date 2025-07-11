# mondo2omop
This repository provides an example script for combining the Mondo Disease Ontology with the OMOP Concept and Concept Relationship tables, utilizing the precise mappings from Mondo to other terminologies, to produce a Mondo-to-OMOP mapping table that can be utilized to connect Mondo to an OMOP instance and generate Mondo patient cohorts.

# Prerequisites
You will need to have previously downloaded and prepared, or otherwise made available, the OMOP CONCEPT and CONCEPT_RELATIONSHIP tables, or made copies available from your local OMOP instance, and place copies of those tables in data/omop. If you don't yet have these OMOP tables available, they can be downloaded from https://athena.ohdsi.org/vocabulary/list. Be sure that the ICD10CM, SNOMED, and MeSH vocabularies are included in the download in order to utilize those mappings from Mondo.

