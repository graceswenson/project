# Siberian Wolly Mammoth Phylogenetic Analysis
- Name of Article: Intraspecific phylogenetic analysis of Siberian woolly mammoths using complete mitochondrial genomes
- DOI:  https://doi.org/10.1073/pnas.0802315105
- Link to Article: https://www.pnas.org/doi/10.1073/pnas.0802315105#supplementary-materials

## Species Information
The Siberian Wolly Mammoth is an extinct species of elephant that roamed parts of the Northern Hemisphere during the Pleistocene epoch. Their remains have been remarkably well-preserved in the permafrost of Siberia, allowing for easier study of their DNA. In this study, mitochondrial DNA (mtDNA) data sets collected from hair shafts will be used to complete a phylogenetic analysis. Samples from __________ species will be used.

## Data Source
I retrieved complete mitochondrial genome sequences for various Siberian woolly mammoth (Mammuthus primigenius) specimens from the NCBI GenBank database. Each genome was downloaded individually in FASTA format based on accession numbers provided in the published article. The sequences were then compiled into a single raw file and saved as wolly_mammoth_data_raw.fasta.

## Alignment Process
To align the complete mitochondrial genomes, I used the multiple sequence alignment tool ClustalW in the terminal with the following command:
clustalw -ALIGN -INFILE=wolly_mammoth_data_raw.fasta -OUTFILE=wolly_mammoth_aligned_final.fasta -OUTPUT=FASTA

Alignment Score: 2983150
Output File: wolly_mammoth_aligned_final.fasta

ClustalW
- Description: ClustalW is a widely used progressive alignment algorithm that aligns biological sequences by first conducting pairwise alignments, constructing a guide tree based on sequence similarity, and then progressively aligning sequences according to the guide tree.
- Assumptions: Sequences with higher similarity are more closely related evolutionarily. Evolutionary events like substitutions, insertions, and deletions occur independently. Gap penalties approximate biological realities.
- Limitations: Errors introduced in early pairwise alignments can propagate through the final alignment. Predefined gap penalties may not always reflect true evolutionary constraints. For large datasets, ClustalW is generally slower and less accurate compared to newer tools like MUSCLE.
- Source for ClustalW: https://doi.org/10.1093/nar/22.22.4673

## Phylogenetic Tree Construction
I constructed the phylogenetic tree of wolly mammoth speciments using RAxML-NG and visualized the tree using R. The process involved obtaining a multiple sequence alignment (MSA) of mitochondrial DNA sequences, constructing the phylogenetic tree using maximum likelihood (ML) estimation, and visualizing the tree for further analysis. For tree construction, I chose to use RAxML-NG because it offers a robust implementation of maximum likelihood estimation. I applied the GTR+G model to infer the best-fit tree, accounting for nucleotide substitution and gamma-distributed rate heterogeneity across sites. Once I had inferred the phylogenetic tree, I used R to visualize it. Specifically, I used the ape package, which is excellent for handling and plotting phylogenetic trees. I read the tree file into R, rooted it with an outgroup, and added labels for nodes and edges to enhance the tree’s clarity.

Code: 
./raxml-ng --check --msa wolly_mammoth_aligned_final.fasta --model GTR+G
./raxml-ng --msa wolly_mammoth_aligned_final.fasta --model GTR+G --prefix wolly_mammoth_tree --threads 2 --seed 2
install.packages("ape")
library(ape)
tree <- read.tree("wolly_mammoth_tree.raxml.bestTree")
plot(tree, main = "Woolly Mammoth Phylogenetic Tree")
tree_rooted <- root(tree, outgroup = "AB015094.1")
plot(tree_rooted, main = "Rooted Woolly Mammoth Phylogenetic Tree")
nodelabels()
edgelabels()
pdf("woolly_mammoth_phylogeny.pdf")
plot(tree_rooted, main = "Rooted Woolly Mammoth Phylogenetic Tree")
dev.off()

Tree saved as "woolly_mammoth_phylogeny.pdf"

RAxML-NG and R:
- Assumptions: I assumed that the GTR+G model would be appropriate for this analysis, as it accounts for nucleotide substitutions and gamma-distributed rate heterogeneity across sites. This model is commonly used for mitochondrial DNA sequences. I assumed that the outgroup I selected was suitable for rooting the phylogenetic tree.
- Limitations: The quality of the sequence alignment is paramount for constructing an accurate tree. Any errors in the alignment could lead to inaccuracies in the inferred phylogeny. I took care to ensure that the alignment was as accurate as possible, but errors can still occur, especially with large datasets. While I used the GTR+G model, there is always the possibility that a different substitution model could fit the data better. The GTR model is widely used and often produces good results, but it may not always be the optimal choice for every dataset.I also relied on bootstrap values to assess the support for the branches of the tree. However, if the bootstrap values are low, this could indicate that the tree is not well-supported and may not accurately reflect the true evolutionary relationships between the species.
- Strengths: RAxML-NG: This software is highly efficient for large datasets and produces phylogenetic trees using maximum likelihood estimation, a widely accepted and robust method for tree construction. It also supports bootstrapping, which provides an estimate of branch support. R Visualization: The ape package in R is incredibly versatile for visualizing phylogenetic trees. It allows me to customize the appearance of the tree, add labels, and export high-quality graphics for publication or presentation. This flexibility is essential for conveying my results clearly and accurately.

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
