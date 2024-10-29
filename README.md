# Random Walk with Restart on Multiplex Network
This project employs the Random Walk with Restart (RWR) algorithm on a multiplex network, with -log10(P-value) used as the initial score vector. Each network consists of three layers:

1. Protein-Protein Interaction (PPI) Layer: Contains genes/proteins directly interacting with disease-associated genes (1st shell) and those interacting with 1st shell genes/proteins (2nd shell).
2. KEGG Pathways Layer: Represents pathways in which the disease genes participate.
3. Co-expression Layer: Derived from RNA expression data from 109 immune cell samples and myositis/immune-related tissues, including bone marrow, spleen, thymus, tonsils, skin, lymph node, and skeletal muscle.

## Permutation Testing
We conducted 1000 RWR iterations by randomly permuting seed nodes. Candidate p-values were calculated based on the proportion of permuted RWR scores that were equal to or exceeded the observed RWR score.

## Running Sample Code
To perform the RWR and Permutation Test, use the following command structure:
```
Rscript RWR_with_Permutation_Test.R <disease_name> <PPI_input_file> <KEGG_input_file> <coexpression_network_input_file> <output_directory>
```
## Parameters
 - disease_name: Name of the disease for analysis (e.g., "PM")
 - PPI_input_file: Path to the PPI adjacency matrix file
 - KEGG_input_file: Path to the KEGG pathway gene adjacency matrix file
coexpression_network_input_file: Path to the co-expression network adjacency matrix file
 - output_directory: Directory where results will be saved
 - 
## Algorithm Reference
Our paper, Meta-analyses Uncover the Genetic Architecture of Idiopathic Inflammatory Myopathies, forms the foundation of this analysis (accepted). Citation will be provided soon.

The RWR algorithm implementation is based on:

A Valdeolivas, L Tichit, C Navarro, S Perrin, G Odelin, N Levy, P Cau, E Remy, and A Baudot. 2018. “Random walk with restart on multiplex and heterogeneous biological networks.” Bioinformatics 35 (3).

Our code includes modifications based on the [source code](https://github.com/alberto-valdeolivas/RWR-MH/tree/master/RWR-M.zip), with adjustments documented in the Supplementary Information, such as:

 - Layer Normalization: Performed on each layer of the multiplex network.
 - Degree-degree Spearman Correlation: Calculated for identical nodes across layers.
 - Permutation Test: Conducted to assess the significance of RWR scores.
