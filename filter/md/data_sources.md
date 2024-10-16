# Data Sources

All influenza virus sequence data utilized in Flu-GLUE were sourced from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore).

## Influenza Virus Sequence Data Processing

To improve the data quality in Flu-GLUE, we performed extensive processing of GenBank sequence data, as described below.

### Validation and Normalization

For each influenza virus species, an extensive validation and normalization process was applied to the sequence and isolate data associated with GenBank entries. The following steps were performed:

1. **Grouping Entries**: Entries from the same isolate were grouped together.
2. **Validation of Segment Numbers**: The segment number of each individual sequence was validated, addressing potential missing or incorrect information in GenBank.

### Isolate-Associated Information

- **Recovery Process**: In cases where isolate-associated information was missing, a comprehensive investigation of sequences and their associated documentation (e.g., scientific publications) was conducted to recover this information.

### Standardization of Host Species Information

- **Taxonomic Nomenclature**: Information about the host species from which influenza A virus (IAV) isolates were obtained, as recorded in GenBank entries, was standardized using correct taxonomic nomenclature for the species. This was done in accordance with taxonomic definitions in the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) database.

- **Ambiguous Species Identification**: In instances where information about the species of isolation was ambiguous (e.g., the host is given as 'duck' without specifying the exact duck species), correct scientific nomenclature for the species group was used. For example, where the species was listed as 'duck,' the taxonomic family to which ducks belong ('Anatidae') was substituted.

### Tabular File Production

A tabular file containing information about all IAV complete genome sequences was produced for each influenza virus species - see [this directory](https://github.com/giffordlabcvr/Flu-GLUE/tree/main/tabular/genus). The files include standardized isolate information and the GenBank accession number of each segment.

### Subtyping and Clade Definition for IAV

- **IAV Isolate Subtyping**: IAV isolates were divided into subtypes based on the hemagglutinin (H) and neuraminidase (N) genes, following established conventions. A BLAST-based genotyping procedure was performed to validate the subtype of each IAV isolate in the complete genome set.

- **Clade Definition**: For each of the other IAV segments (i.e. those that do not encode hemagglutinin or neuraminidase) clades were defined on the basis of an extensive phylogenetic investigation. A BLAST-based genotyping was performed for each complete genome isolate to provide taxonomic information for these six segments.

This thorough data processing ensures the integrity, accuracy, and completeness of the influenza virus sequence data contained in Flu-GLUE.

