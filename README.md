#Intercellular communincation based on scRNAseq data in healthy and UC colon
Codes developed by Lejla Gul and Dezso Modos

##Aim
Looking for protein-protein interactions between cells where the expression of one participant is
changing during UC, therefor selecting the potential differences in cell-cell communication in diseased state

##Required packages
[pandas](https://pandas.pydata.org/),
[itertools](https://docs.python.org/3.1/library/itertools.html)

##Workflow
1. Importing pre-processed single-cell RNAseq data

2. Selecting expressed genes in each cell type (separately in healthy and non-inflamed UC conditions)

3. Filtering to genes/proteins participating in intercellular communication based on OmniPath

4. Pairwise comparison of cells

##Pre-processing scRNAseq data
Raw 10X data published by [Smillie et al](https://pubmed.ncbi.nlm.nih.gov/31348891/), pre-processing steps described in
the article.

##Transformation of pre-processed, average expression data
To select the expressed genes, we have used a published
[zFPKM method](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-778) to get the expressed genes
in different cell types. Then, to avoid the huge differences among values, log2 transformation has been applied.

##Get interaction from OmniPath using OmniPathR
import_intercell_network() function downloads the interactions.
With the get_intercell_classes() function information is available about type of intercellular proteins.
Here we selected membrane based or secreted proteins from these classes,
only extracellular matrix, matrix adhesion, matrix adhesion regulator,ecm regulator, ligand regulator
and receptor regulator categories have been discarded.

##Input of the pipeline
1. Tables describing the average gene expression in epithelial cells, immune cells and fibroblasts
(healthy, inflamed and non-inflamed UC)

2. Intercellular interactions from OmniPath

##Output of the pipeline
The outputs of the cell-cell communication are .txt files which describe the condition specific interactions
(one part of the intercellular communication appears only in one condition - healthy or non-inflamed UC).
