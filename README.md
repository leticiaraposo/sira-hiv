# SIRA-HIV

SIRA-HIV is an application for identifying HIV drug resistance from a next generation sequencing data. This tool can detect minority mutations in HIV sequences with a 1% cut-off. [Segminator II](http://www.bioinf.manchester.ac.uk/segminator/) is used for the characterization of viral read data. This software supports the next generation sequencing platforms including: 454 Life Sciences, Illumina, Ion Torrent and Pacific Biosciences. The output of SIRA-HIV is a drug resistance report which provides the frequencies of the amino acids found in the regions of gene Pol and the susceptibility of HIV-1 variants to antiretroviral drugs by genotypic interpretation systems.
                           
### Prerequisites

Before running SIRA-HIV it is necessary to have RStudio and Java installed on your computer. Click [here](https://www.rstudio.com/products/rstudio/download/) to download RStudio and [here](https://www.java.com/download/) to download Java.

### Installing

To install SIRA-HIV, just download the zip file here in GitHub, and then extract the files into a folder.
Inside the folder, open the file ```SIRA-HIV.R``` and click on ```Run App```.
SIRA-HIV will automatically open in a new window. It is recommended to open SIRA-HIV in a browser to have a complete experience of the program. 

![runapp](https://user-images.githubusercontent.com/14352761/49885138-a628d280-fe1d-11e8-9855-ea9989be780c.png)

## Running 

SIRA-HIV can be run from raw data output of next-generation sequencing (NGS) :one: or from pre-processed data by the Segminator II :two:.

![drugresistancepositions](https://user-images.githubusercontent.com/14352761/49885571-c60cc600-fe1e-11e8-93b3-22660236bef1.png)

### Raw Data Output of NGS
Segminator II must be initialized.
Before inserting the files, it is recommended to adjust the alignment parameters. These can be found in the ```Parameters``` menu. It is recommended to select the parameter ```Replace Template With Con. During Mapping``` before inserting the files so that the data is remapped with the generated consensus sequence.

![parameters](https://user-images.githubusercontent.com/14352761/49884565-4b42ab80-fe1c-11e8-9dd7-7b955e700dd2.png)

:one: In the ```Project``` menu, the user provides the [FASTA file](sira-hiv/data/ReferenceSequence.fasta) containing the reference genome and in the ```Data``` menu, the raw FASTQ file is inserted.

:two: Selecting the ```Project``` and the ```Data```, Segminator II performs an automatic alignment.

:three: The results can be exported using the ```Tools -> VEME Table``` menu.

![segminator](https://user-images.githubusercontent.com/14352761/49884407-f1da7c80-fe1b-11e8-8490-2752f6510a73.png)

### Preprocessed Data by Segminator II
After obtaining the VEME Table file, SIRA-HIV is ready to generate its results. 
It is just necessary to insert the VEME Table file and choose the region of the pol gene to be analyzed.
SIRA-HIV performs an automatic analysis of the data.

There is a [VEME Table](sira-hiv/data/VEMETableExample.txt) and [FASTQ](sira-hiv/data/SampleExample.fastq) example files available in sira-hiv/data.

## Authors

* **Letícia Martins Raposo** 
* **Flavio Fonseca Nobre** 

## Acknowledgments

* Biomedical Engineering Program, Alberto Luiz Coimbra Institute for Graduate Studies and Research in Engineering, Federal University of Rio de Janeiro (UFRJ)
* Laboratory of Molecular Virology, Dept. of Genetics, Federal University of Rio de Janeiro (UFRJ)
* Coordination for the Improvement of Higher Education Personnel (CAPES) 
* National Council for Scientific and Technological Development (CNPq)
* Carlos Chagas Filho Foundation for Research Support in the State of Rio de Janeiro (FAPERJ)
