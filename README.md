# Random Walk with Restart on Multiplex Network
We utilized -log10(P-value) as the initial score vector for the Random Walk with Restart(RWR) Algorithm on Multiplex Network. For the total IIMs, PM, DM, JDM, and myositis with anti-Jo1, we construted a multiplex network for each of those. Each Multiplex Network was compsoed of three layers including:

1. Protein-Protein Interaction of genes/proteins that directly interacted with the disease genes(1st shell) and those that interacted with the disease genes or proteins from the 1st shell(2nd shell)
2. Co-expression Layer generated from RNA transcript expression data of 109 immune cell samples and myositis and immune-system related tissues, including bone marrow, spleen, thymus, tonsils, skin, lymph node, and skeletal muscle
3. KEGG Pathways in which disease genes participate.

The RWR code was based upon the algorithm published in:

A Valdeolivas, L Tichit, C Navarro, S Perrin, G Odelin, N Levy, P Cau, E Remy, and A Baudot. 2018. “Random walk with restart on multiplex and heterogeneous biological networks.” Bioinformatics 35 (3). DOI: https://doi.org/10.1093/bioinformatics/bty637

Modifications were performed on the source code (https://github.com/alberto-valdeolivas/RWR-MH/tree/master/RWR-M.zip) as described in the Supplementary, including:
1. Normalization on each layer of the multiplex network
2. Degree-degree Spearman correlation between the same nodes in different layers

# Permutation Test
We conducted 1000 iterations of RWR by randomly permuting the seed nodes. The p-values of the candidates were determined by calculating the proportion of permuted RWR scores that were equal to or greater than the observed RWR score.
