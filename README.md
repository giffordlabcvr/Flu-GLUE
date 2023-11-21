# Flu-GLUE: Comparative Genomic Analysis of Influenza Viruses

## Overview

Flu-GLUE is a sequence-oriented resource for comparative genomic analysis of influenza viruses, developed using the GLUE software framework. It provides a freely accessible, public-facing platform that offers a phylogenetically-structured collection of all publicly available influenza sequence data. Sequence data are richly annotated with gene features and isolate-associated information.

GLUE (**G**enes **L**inked by **U**nderlying **E**volution) is a data-centric bioinformatics environment for virus sequence data, with a focus on variation, evolution and sequence interpretation.

You can learn more about it on the [GLUE web site](http://glue-tools.cvr.gla.ac.uk).

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Usage](#usage)
- [Data Sources](#data-sources)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Key Features

- **GLUE Framework Integration**: Built on the GLUE software framework, Flu-GLUE offers an extensible platform for efficient, standardized, and reproducible computational genomic analysis of influenza viruses.

- **Phylogenetic Structure**: The data in Flu-GLUE is organized in a phylogenetically-structured manner, allowing users to explore evolutionary relationships easily.

- **Rich Annotations**: Each sequence is annotated with gene features, enabling rigorous comparative genomic analysis related to conservation, adaptation, structural context, and genotype-to-phenotype associations.

## Installation

If you have not done so already, install the GLUE software framework by following the [installation instructions](http://glue-tools.cvr.gla.ac.uk/#/installation) on the GLUE web site. 

Assuming you have installed GLUE correctly, you can now build the Flu-GLUE project. Note the Flu-GLUE project has a layered structure. This approach simplifies project management because it allows data items that are likely to be used across a wide range of analysis contexts to be maintained separately from those only required for more specialized purposes. The ‘base’ layer of Flu-GLUE contains only a minimal set of essential data items required for comparative analysis of influenzaviruses.

To build the base (or 'core') project, start the GLUE command line interpreter, and at the GLUE command prompt, run the 'buildCoreProject.glue' file as follows:

`GLUE> run file buildCoreProject.glue`

Once the core project has been constructed, influenza virus species-specific extension layers containing larger numbers of sequences can be constructed. 

First, however, you should download these sequences from GenBank via GLUE. For three influenza virus species - influenza B virus (IBV),  influenza C virus (ICV), and  influenza D virus (IDV), the number of GenBank entries is small enough that all sequences can be downloaded. For example:

`GLUE> project flu run file glue/build/genus/icv/icvDownloadCurated.glue`

For influenza A virus (IAV) the number of available GenBank sequence entire is much larger than for other influenza virus species. Accordingly, we have implemented a ‘plug-and-play’ approach to IAV analysis wherein researchers separately download sequence entries for each IAV subtype, using the files in [this directory](https://github.com/giffordlabcvr/Flu-GLUE/tree/main/glue/build/genus/iav/download/).

`GLUE> project flu run file glue/build/genus/iav/download/downloadNcbiSequencesIavH5N1.glue`

Once segment sequences have been downloaded, they can be exported and stored on your hard drive. Use the following GLUE commands to export sources for each segment:

`GLUE> project flu export source icv-ncbi-curated-segment-1`

`GLUE> project flu export source icv-ncbi-curated-segment-2`

...etc. Remember that ICV & IDV have seven segments, whereas IAV and IBV have eight!

These commands will export directories containing the curated, segment-specific sequences for each influenza virus species or IAV subtype. Note that for IAV, curated sequences are restricted to complete genome isolates (i.e. isolates for which each individual segment has been sequenced).

Store the exported source directories in a sensible place on your hard drive (e.g. the 'sources' folder in the Flu-GLUE project folder), then update the input 'LoadSources.glue' file for the relevant influenza virus species so that when the extension layer is built, these sequences will be loaded. 

For example, for ICV, update [this file](https://github.com/giffordlabcvr/Flu-GLUE/blob/main/glue/build/genus/icv/icvLoadSources.glue). The 'import source' commands for each segment should point to the relevant directory. For example, to import the 'segment 1' sequences for ICV:

`import source path/to/source/folder/icv-ncbi-curated-segment-1`

Once the 'LoadSources.glue' file for the relevant influenza virus species has been updated, the extension layer for that virus can be built by running the relevant extension layer build file. For example, for ICV:

`run file glue/build/genus/icv/buildIcvExtension.glue`

## Usage

For detailed instructions on how to use Flu-GLUE for comparative genomic analysis of influenza viruses, refer to the GLUE software environment's [reference documentation](http://glue-tools.cvr.gla.ac.uk/).

## Data Sources

Flu-GLUE relies on the following data sources:

- [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore)
- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

Additional flu sequence data can be imported from FASTA files. Each FASTA file should have a unique sequence ID and be stored in a separate file. FASTA data can be imported by using GLUE's '[fastaImporter](http://glue-tools.cvr.gla.ac.uk/#/moduleReference/moduleType/fastaImporter)' module type.

For more information on the data sources and how they are integrated, see [Data Sources](./md/data_sources.md).

## Contributing

We welcome contributions from the community! If you're interested in contributing to Flu-GLUE, please review our [Contribution Guidelines](./md/CONTRIBUTING.md).

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](./md/code_of_conduct.md) 

## License

The project is licensed under the [GNU Affero General Public License v. 3.0](https://www.gnu.org/licenses/agpl-3.0.en.html)

## Contact

For questions, issues, or feedback, please contact us at [robert.gifford@glasgow.ac.uk](mailto:robert.gifford@glasgow.ac.uk) or open an issue on the [GitHub repository](https://github.com/giffordlabcvr/Flu-GLUE/issues).


