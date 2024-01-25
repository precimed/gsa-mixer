https://www.syngoportal.org 
SynGO dataset version/release: 20210225

syngo_genes.xlsx contains all genes in the SynGO database.

syngo_annotations.xlsx lists all individual annotation statements.

syngo_ontologies.xlsx lists all ontology terms and their respective annotated genes.

In case you are matching/integrating this data with other GO resources, note that some SynGO ontology terms are not part of the GO database, their identifiers are easily recognised by the "SYNGO:" prefix as opposed to "GO:". Annotations against these terms is integrated into GO, but translated into an annotation+extensions as some of the SynGO terms do not comply with GO guidelines (eg; a term like "presynaptic processes" is too generic for GO as they'd end up with thousands of suchs terms).

The .json files are intended for bioinformatic application of the SynGO dataset.
The SynGO_geneset_CC/BP.json files contain all SynGO annotated genes for each synaptic ontology term. On the first level of this nested list you can find 'direct' and 'aggregate' annotations, you will want to use the latter for geneset analyses as these contain (for each term) both the annotated genes for an ontology term and all genes annotated against (recursive) child terms. For instance, in the aggregate dataset the term 'presynapse' also contains all genes annotated against its child term 'active zone'. Further, note that we originally mapped the annotated proteins from various species to HGNC identifiers. The ensembl/entrez/MGI/RGD mappings (from HGNC ID) in these JSON files were provided by the genenames.org webservice, so if you focus on Ensembl genes consider mapping the HGNC IDs to the exact Ensembl build that you are using.
Do consider using the SynGO color-coding tool to visualize your own SynGO ontology term summary scores.