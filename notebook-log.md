# Siberian Wolly Mammoth Phylogenetic Analysis
- Name of Article: Intraspecific phylogenetic analysis of Siberian woolly mammoths using complete mitochondrial genomes
- DOI:  https://doi.org/10.1073/pnas.0802315105
- Link to Article: https://www.pnas.org/doi/10.1073/pnas.0802315105#supplementary-materials

## Species Information
The Siberian Wolly Mammoth is an extinct species of elephant that roamed parts of the Northern Hemisphere during the Pleistocene epoch. Their remains have been remarkably well-preserved in the permafrost of Siberia, allowing for easier study of their DNA. In this study, mitochondrial DNA (mtDNA) data sets collected from hair shafts will be used to complete a phylogenetic analysis. Samples from thirteen species will be used.

## Data Source
I retrieved complete mitochondrial genome sequences for various Siberian woolly mammoth (Mammuthus primigenius) specimens from the NCBI GenBank database. Each genome was downloaded individually in FASTA format based on accession numbers provided in the published article. The sequences were then compiled into a single raw file and saved as wolly_mammoth_data_raw.fasta.

## Alignment Process
To align the complete mitochondrial genomes, I used the multiple sequence alignment tool ClustalW in the terminal with the following command:
```
clustalw -ALIGN -INFILE=wolly_mammoth_data_raw.fasta -OUTFILE=wolly_mammoth_aligned_final.fasta -OUTPUT=FASTA
```

Alignment Score: 2983150
Output File: wolly_mammoth_aligned_final.fasta

ClustalW
- Description: ClustalW is a widely used progressive alignment algorithm that aligns biological sequences by first conducting pairwise alignments, constructing a guide tree based on sequence similarity, and then progressively aligning sequences according to the guide tree.
- Assumptions: Sequences with higher similarity are more closely related evolutionarily. Evolutionary events like substitutions, insertions, and deletions occur independently. Gap penalties approximate biological realities.
- Limitations: Errors introduced in early pairwise alignments can propagate through the final alignment. Predefined gap penalties may not always reflect true evolutionary constraints. For large datasets, ClustalW is generally slower and less accurate compared to newer tools like MUSCLE.
- Source for ClustalW: https://doi.org/10.1093/nar/22.22.4673

## Phylogenetic Tree Construction
### Method 1: Distance Based Tree Using R
I constructed a phylogenetic tree of woolly mammoth specimens using Neighbor-Joining (NJ) methods in R. This approach involved maximum likelihood estimation using the GTR+G model, which accounts for nucleotide substitution and gamma-distributed rate heterogeneity across sites. Due to file size and computational limitations, the dataset was split into two parts for analysis.

Code: 
```r
library(ape)                              
part1_dnabin <- as.DNAbin(part1)
part2_dnabin <- as.DNAbin(part2)
dist1 <- dist.dna(part1_dnabin, model = "raw")
dist2 <- dist.dna(part2_dnabin, model = "raw")
tree1 <- nj(dist1)
tree2 <- nj(dist2)
par(mfrow = c(1, 2))
plot(tree1, main = "Tree from Part 1")
plot(tree2, main = "Tree from Part 2")
combined_dnabin <- c(part1_dnabin, part2_dnabin)
dist_combined <- dist.dna(combined_dnabin, model = "raw")
tree_combined <- nj(dist_combined)
plot(tree_combined, main = "Consensus Tree from All Sequences")
```
Tree saved as "wolly_mammoth_NJ_tree.pdf"

- Assumptions: All input sequences are properly aligned and of mitochondrial origin. NJ method assumes equal evolutionary rates across lineages (no molecular clock). The raw distance model assumes equal mutation weights and does not correct for multiple substitutions.
- Limitations: File size (215 KB) made direct computation slow; splitting into parts was necessary. Less accurate than model-based methods (e.g., Maximum Likelihood or Bayesian). If sequence quality or alignment errors are present, trees may reflect noise instead of true evolutionary relationships.
- Strengths: NJ is computationally inexpensive and appropriate for large alignments. Splitting the dataset allowed analysis to proceed despite memory constraints. Side-by-side trees allow comparison of topology and potential clustering.

## Method 2: RAxML and R
I constructed a phylogenetic tree of woolly mammoth specimens using RAxML, a tool for maximum likelihood-based phylogenetic analysis. The tree was estimated using the GTR+G model, which accounts for nucleotide substitution and gamma-distributed rate heterogeneity across sites. The analysis was performed on the full dataset of aligned sequences.

RAxML in Terminal:
The sequence data was first modified to ensure unique sequence names and then processed using RAxML in the terminal. The following command was used to generate the tree:
```
/usr/local/Cellar/raxml/8.2.12/bin/raxmlHPC-AVX -s /Users/graceswenson/Desktop/563/modified_wolly_mammoth_aligned_final_unique.fasta -n my_tree_run3 -m GTRGAMMA -p 12345
```

Code: R Visualization of Phylogenetic Tree
Once the tree was generated by RAxML, it was visualized in R using the ape package. The following steps were taken:
```
install.packages("ape") 
library(ape) tree <- read.tree("/Users/graceswenson/Desktop/563/RAxML_bestTree.my_tree_run3")
plot(tree)
```

Tree saved as: "wolly_mammoth_raxml_tree.pdf"

- Assumptions: The input sequences are aligned and of mitochondrial origin. The GTR+G model was selected to account for nucleotide substitutions and site-rate heterogeneity. RAxML uses maximum likelihood estimation, which provides accurate tree topologies for large datasets.
- Limitations: The quality of the alignment is crucial. Any misalignment or sequence errors can affect the results. The RAxML analysis assumes no molecular clock, meaning it does not assume constant rates of evolution across all branches. Computationally intensive for large datasets, but the model offers a good balance of accuracy and computational feasibility.
- Strengths: RAxML provides robust maximum likelihood trees, widely used for large-scale phylogenetic analysis. The GTR+G model accounts for rate heterogeneity, making it suitable for diverse evolutionary scenarios. Visualizing the tree in R helps assess the quality and topology of the tree and provides a straightforward interpretation.

## Bayesian Inference
In addition to maximum likelihood analysis, I chose to implement a Bayesian inference approach using MrBayes to construct a phylogenetic tree. Bayesian methods offer an alternative statistical framework for estimating phylogenies and allow for the incorporation of prior probabilities and model parameters in a probabilistic way. This method complements the results from RAxML-NG and provides a deeper understanding of the phylogenetic relationships among woolly mammoth mitochondrial genomes.

Method:
1. To prepare for MrBayes, I first converted my aligned FASTA file (wolly_mammoth_aligned_final.fasta) into the NEXUS format, which is the required format for MrBayes. I used the seqinr package in R to do this. 
2. Once I had the .nex file ready, I used the following MrBayes commands in the terminal or within the MrBayes shell:
mb
execute wolly_mammoth_aligned_final.nex;
lset nst=6 rates=gamma;
prset statefreqpr=dirichlet(1,1,1,1);
mcmc ngen=1000000 samplefreq=100 printfreq=100 diagnfreq=1000 nchains=4 temp=0.2;
sump burnin=2500;
sumt burnin=2500;
This code runs the MCMC algorithm for 1 million generations, sampling every 100 generations. I then discarded the first 2,500 trees as burn-in and summarized the remaining trees to construct a consensus phylogeny with posterior probabilities.

Assumptions: Bayesian inference requires setting priors for evolutionary parameters. I used non-informative (flat) priors such as a uniform Dirichlet distribution for state frequencies, assuming no prior knowledge about base frequency distribution.  I assumed the MCMC chains converged and mixed well based on log-likelihood plots and diagnostic statistics like the potential scale reduction factor (PSRF) and effective sample size (ESS). I assumed the GTR+G model (nst=6, rates=gamma) was appropriate for my mitochondrial sequence data, as this is widely accepted for DNA sequence evolution.

Strengths: MrBayes provides posterior probabilities for each node, which can be more interpretable than bootstrap values because they directly represent the probability of clades given the model and data. MrBayes allows for complex models of evolution, including mixture models and codon-based models, which could be useful in future studies with more complex datasets.

Limitations: MrBayes is computationally intensive, especially with large datasets or when running very long chains. For this analysis, I needed to allow significant time for convergence. Choosing an appropriate burn-in is crucial. If the chains have not converged, then the posterior estimates may be biased. I visually inspected the likelihood trace and checked convergence diagnostics but this step requires careful interpretation. As with other methods, the results can be sensitive to the chosen substitution model. While I used GTR+G, different models may yield different topologies.
