# allele.varibility
The Rscript allele.variability.R quantifies allele variability of KEGG annotated gene orthologs across multiple metagenome data and estimate the influence of meta community size, average per cell gene copy number and sub cellular location on allele variability using bootstrapped statistics as detailed in Beier et al. 2020.

#### Parameters to set:
```bash
n.min #minimum number of reads per metagenome mapping onto each considered gene orthodox
n.boot #number of permutations for bootstrapping statistics (warning: high number of permutations cause enhanced computing time and correspondingly the number of threads should be adjusted)
n.thread #number of threads for parallelization of the bootstrapping approach
```

#### Input file 1: 
tab-delimited count-table  with the first column indicating gene ortholog IDs, the second column indicating allele IDs and the third to nth column representing the number of reads mapped per allele ID in each metagenome (e.g. counts.tab)

#### Input file 2: 
tab-delimited text file with information about the average per cell copy number and the occurrence of gene orthologs in prokaryotic genomes (e.g. TableS3.tab)

#### Citation
Beier S, Andersson AF, Galand PE, Hochart C, Logue JB, McMahon K, Bertilsson S. 2020. The environment drives microbial trait variability in aquatic habitats. Mol Ecol 29:4605–4617.

