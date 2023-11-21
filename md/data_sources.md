# Data Sources

## Influenza Virus Sequence Data Processing

All influenza virus sequence data utilized in Flu-GLUE were sourced from GenBank.

### Validation and Normalization

For each influenza A virus species, an extensive validation and normalization process was applied to the sequence and isolate data associated with GenBank entries. The following steps were performed:

1. **Grouping Entries**: Entries from the same isolate were grouped together.
2. **Validation of Segment Numbers**: The segment number of each individual sequence was validated, addressing potential missing or incorrect information in GenBank.

### Isolate-Associated Information

- **Recovery Process**: In cases where isolate-associated information was missing, a comprehensive investigation of sequences and their associated documentation (e.g., scientific publications) was conducted to recover this information.

### Standardization of Host Species Information

- **Taxonomic Nomenclature**: Information about the host species from which influenza A virus (IAV) isolates were obtained, available in GenBank, was standardized using correct taxonomic nomenclature for the species. This was done in accordance with taxonomic definitions in the NCBI taxonomy database.

- **Ambiguous Species Identification**: In instances where information about the species of isolation was ambiguous (e.g., the host is given as 'duck' without specifying the exact duck species), correct scientific nomenclature for the species group was used. For example, if the species was listed as 'duck,' 'Anatidae' - the taxonomic family to which ducks belong - was substituted to the highest level of accuracy possible from the available GenBank entry information.

### Tabular File Production

A tabular file containing information about all IAV complete genome sequences was produced. This includes isolate information and the GenBank accession number of each segment.

### Subtyping and Clade Definition

- **IAV Isolate Subtyping**: IAV isolates were divided into subtypes based on the hemagglutinin (H) and neuraminidase (N) proteins on the virus surface. A BLAST-based genotyping procedure was performed to validate the subtype of each IAV isolate in the complete genome set.

- **Clade Definition**: Clades were defined for the other six IAV segments, providing taxonomic information for all segments of complete genome IAV isolates.

This thorough data processing ensures the integrity, accuracy, and completeness of the influenza virus sequence data used in Flu-GLUE for comparative genomic analysis.
