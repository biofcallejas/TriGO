# TriGO (*Trichomonas vaginalis* GO enrichment analysis)

GO enrichment analysis for Trichomonas vaginalis using hypergeometric distribution ([LightGOEA](https://github.com/biofcallejas/LightGOEA) function).
This script was written to calculate GO enrichment for *Trichomonas vaginalis* datasets. However, the script will work with other species files that meet the input formats.

## Installation

It runs in Python 3 and has been tested in version 3.8.12. Using a conda environment is encouraged.
The main script does not require installation but has some Python dependencies (all of them can be easily installed using pip or conda). 

* [Pronto](https://pronto.readthedocs.io/en/stable/) tested on version 2.4.4
* [Numpy](https://numpy.org/) tested on version 1.22.4.
* [Scipy](https://scipy.org/) tested on version 1.12.
* [Pandas](https://pandas.pydata.org/) tested on version 1.3.4.
* [Tabulate](https://pypi.org/project/tabulate/) tested on version 0.8.9.

## GO database

The script needs the [GO database](https://geneontology.org/docs/download-ontology/) (.obo format) to be available in the working directory. 
It can be downloaded [here](https://geneontology.org/docs/download-ontology/) (**go.obo**).

## Usage

* The input file **(-i)** is a two-column TSV file. The first column is the gene_id, and the second column is the list of GO terms associated with it (separated by commas). No headers, no extra/blank spaces.
This file contains all the genes with one or more GO terms. **This is the database, not the test file**. 
* The geneset file **(-g)** is a single-column file with the list of gene IDs to calculate the enrichment. One gene_id per line (no headers or blank spaces).
This is the **test** file.
* The output prefix **(-o)** is a prefix for the output file. 
 
```
Usage enrichmentAnalysisGO.py [-h] -i INPUT -g GENESET -o OUTPUT

	Calculate GO enrichment analysis using a hypergeometric distribution.

Optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUT, --input INPUT
                        Input database, two columns; GO terms and list of geneIDs associated (no headers)
			
  -g GENESET, --geneset GENESET
  
                        List of geneIDs to calculate enrichment, one geneID per line (no headers)
  -o OUTPUT, --output OUTPUT
  
                        Prefix: Output text file. Results are also printed on the screen
```

## Results

A table with the results is printed on the screen.

The output is a seven-column TSV file summarizing the results. 

* **GO** enriched GO term. 
* **GO_level** GO level; Molecular Function (MF), Cellular Component (CC), Biological Process (BP).
* **Description** GO name/description.
* **Enrichment** Enrichment value.
* **fdr_q_value** False Discover Rate (*q*value)
* **p_value** *p*-value
* **genes** Genes enriched with this term

## Additional info

Newer versions of Scipy and Numpy might cause running problems.
The script was tested on MACOS and UNIX OS.

