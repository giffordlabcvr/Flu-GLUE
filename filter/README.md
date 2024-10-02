# Flu NCBI Filter

<img src="md/flu-filter-logo.png" align="right" alt="" width="260" />

Welcome to the **Flu NCBI Filter** repository on GitHub!

Large quantities of influenza virus genome sequence data are publicly available via [**GenBank**](https://www.ncbi.nlm.nih.gov/nuccore). However, GenBank serves as a comprehensive repository for genetic information across all domains of life and was not developed specifically for viral data. Consequently, Influenza virus researchers face a number of data integration and comparison challenges when using GenBank, including:

- **Separate Entries for Each Segment**:  GenBank stores sequences for each influenza genome segment separately, which can make it cumbersome for researchers to access and analyse all relevant genomic information for a particular strain.
- **Lack of Standardised Isolate Information**: GenBank lacks a standardized format for recording isolate-associated information, such as the geographical location of sample collection, host species, date of isolation, and clinical data. As a result, researchers may encounter inconsistencies or incomplete information across different entries. This hampers efforts to conduct comprehensive analyses, track viral spread, and understand the epidemiology of influenza strains.
- **Quality Control and Data Verification**: With the vast amount of data submitted to GenBank, maintaining quality control and verifying the accuracy of submitted sequences and associated metadata can be an overwhelming task. Inaccurate or poorly annotated entries may lead to erroneous conclusions in research studies and impede the progress of influenza virus research.

This repository contains code and scripts used to derive influenza virus sequence databases from GenBank. The code introduces a higher level of order to influenza virus seqeuence data, capturing the links between sequences and isolates, and allows validation and standardisation of sequence-associated data.

## Contents of this repository

This repository comprises two main components: 

1. **A console-based PERL program** that uses cleaned, validated influenza virus sequence data, exported from one of the GLUE 'filtering' projects described below, to generate a database of influenza isolates. This program will also:
 - Deal with redundant sequences for individual isolate segments (selects best representative)
 - Check the consistency of isolate-associated metadata across all sequences.
 - Identify isolates that are incomplete (one or more segment sequences are missing) and export these to a separate tabular file.
 - Generate GLUE module definitions and console code that can be used in to selectively import sets of isolates into **[Flu-GLUE](https://github.com/giffordlabcvr/Flu-GLUE)**. 
2. **A set of [GLUE](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2459-9) projects**, one for each influenza virus species, that can be used to:
 - Download influenza virus sequence data from GenBank (in GenBank XML format).
 - Extract isolate and sequence-associated metadata from GenBank XML files.
 - Validate and standardise sequence data fields.
 - Perform rapid genotyping for all genome segments.
 - Export validated, standardised metadata for GenBank influenza virus sequence entries.

## PERL program

Two PERL programs are provided - one for processing IAV and IBV, and one for processing ICV and IDV. Run the programs with the -h option to see help information.

`./iav-or-ibv-processing.pl -h`



## GLUE projects

The GLUE projects included in this repository can be used to download and process GenBank sequence entries for IAV, IBV, ICV, and IDV. The cleaned, processed data exported from these GLUE projects is included with this respoitory. Therefore, the only reason you may want to build these GLUE projects is if you wish to update the cleaned data to include GenBank entries published more recently and not included in the exported data.

### Building the GLUE projects

To build the GLUE projects, you must first install [GLUE](glue-tools.cvr.gla.ac.uk).

Once GLUE has been installed, local instances of the filter projects can be created by starting the GLUE console and running the relevant build file as follows:

```
Mode path: /
GLUE> run file buildIavFilterProject.glue
```

This will build an empty version of the project (i.e. containing no sequences). 

Sequences can be imported by first editing and then running the iav source download file. Uncomment one or more of the following lines:

```
#module iavNcbiImporter1970_1990 import
#module iavNcbiImporter1990_2005 import
#module iavNcbiImporter2005_2010 import
#module iavNcbiImporter2010_2015 import
#module iavNcbiImporter2015_2020 import
```

Then run the file from the GLUE console, as follows:

```
Mode path: /
GLUE> project iav_filter 
OK
Mode path: /project/iav_filter
GLUE> run glue/filter/iav/iavDownloadSources.glue
```

...to download batches of sequences with publication dates falling within the year range specified by the module name as shown above. For ICV and IDV, all sequences are downloaded in a single batch.

Because GenBank contains large numbers of entries for IAV and IBV (hundreds of thousands in the case of IAV), sequences are batched by year range for these viruses. In general, it is advisable to limit the number of sequences that are downloaded or processed at once, especially if working on a PC with limited disk space and/or RAM. 

Once sources have been downloaded, they should be exported as shown below, and stored in an appropriate filesystem location.

```
Mode path: /
GLUE> run file glue/idv/idvExportProcessedData.glue
```

Now edit the file `glue/filter/alphainfluenzavirus/iavFilterLoadSources` so that the `import source` statements point to this location. For example:

```
import source /Users/JoeBloggs/sources/flu/ncbi-iav-nuccore-1970-1990
```

Now re-run the project build file from the GLUE console (as shown above) to create a GLUE project that imports and processses sequence sources.

### Processing steps performed in the GLUE project build

The GLUE project build process entails the following processing steps, described in greater detail below:

1. Reference sequence declarations.
2. Import of sequence sources (i.e. folders containing GenBank XML files, downloaded as described above)
3. Extraction of sequence and isolate metadata from GenBank XML.
4. Segment recognition.
5. IAV and IBV only: genotyping of individual segment sequences (IAV), or of segment 4 sequences (IBV).

#### Reference sequence declarations

The GLUE framework uses `ReferenceSequence` objects to organise, link and interpret sequence data within a project. Each of these objects is based on specific sequence included in the GLUE project. Each of the four GLUE projects defined in this repository includes a set of representative reference sequences for the corresponding influenza virus species (IAV, IBV, ICV, or IDV). Here, reference sequences are primarily used for segment recognition and genotyping.

#### Extraction of GenBank data

Isolate and sequence metadata is extracted from GenBank XML using a customised version of GLUE's `GenBankPopulator` module, which can be viewed [here](https://github.com/giffordlabcvr/Flu-NCBI-Filter/blob/main/modules/filter/alphainfluenzavirus/iavFilterGenbankXmlPopulator.xml).

#### Segment recognition

The segment number associated with influenza virus sequences is typically recorded in the 'notes' section of a GenBank XML file. Occasionally, influenza virus sequences are associated with incorrect segment numbers. Each of the filtering projects entails a step in which a customised version of GLUE's `SequenceRecogniser` module is used to independently confirm the segment that each sequence derives from.

#### Genotyping

For IAV and IBV, a genotyping step is performed. For IAV, genotyping is performed on all segments, while for IBV it is only performed on segment 4 sequences.

In IAV, genotypes have traditionally been defined by genotyping of the two virion surface proteins, hemagglutinin (H) and neuraminidase (N) (e.g. ‘H1N1’), which are encoded by segments 4 and 6, respectively.

#### Data export

To export the process GenBank sequence entries, execute the following command on the GLUE console:

```
Mode path: /
GLUE> run file glue/idv/idvExportProcessedData.glue
```

## License

The project is licensed under the [GNU Affero General Public License v. 3.0](https://www.gnu.org/licenses/agpl-3.0.en.html)
