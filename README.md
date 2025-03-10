# TriGO
GO enrichment analysis for Trichomonas vaginalis

usage: enrichmentAnalysisGO.py [-h] -i INPUT -g GENESET -o OUTPUT

	Calculate GO enrichment analysis using hypergeometric distribution
Version 2.0 for python 3 or higher

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUT, --input INPUT
                        Input database, 2 columns; GO terms and list of geneIDs associated (no headers)
  -g GENESET, --geneset GENESET
                        List of geneIDs to calculate enrichment, one geneID per line (no headers)
  -o OUTPUT, --output OUTPUT
                        Prefix: Output text file, results are also printed on screeen
