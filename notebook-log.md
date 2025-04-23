# Siberian Woolly Mammoth Phylogenetic Analysis
- Name of Article: Intraspecific phylogenetic analysis of Siberian woolly mammoths using complete mitochondrial genomes
- DOI:  https://doi.org/10.1073/pnas.0802315105
- Link to Article: https://www.pnas.org/doi/10.1073/pnas.0802315105#supplementary-materials

## Species Information
The Siberian woolly mammoth is an extinct species of elephant that roamed parts of the Northern Hemisphere during the Pleistocene epoch. Their remains have been remarkably well-preserved in the permafrost of Siberia, allowing for easier study of their DNA. In this study, mitochondrial DNA (mtDNA) data sets collected from hair shafts will be used to complete a phylogenetic analysis. Samples from thirteen species will be used.

## Data Source
I retrieved complete mitochondrial genome sequences for various Siberian woolly mammoth (Mammuthus primigenius) specimens from the NCBI GenBank database. Each genome was downloaded individually in FASTA format based on accession numbers provided in the published article. The sequences were then compiled into a single raw file and saved as woolly_mammoth_data_raw.fasta.

## Alignment Process
To align the complete mitochondrial genomes, I used the multiple sequence alignment tool ClustalW in the terminal with the following command:
```
clustalw -ALIGN -INFILE=woolly_mammoth_data_raw.fasta -OUTFILE=woolly_mammoth_aligned_final.fasta -OUTPUT=FASTA
```

Alignment Score: 2983150
Output File: woolly_mammoth_aligned_final.fasta

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
Tree saved as "woolly_mammoth_NJ_tree.pdf"

- Assumptions: All input sequences are properly aligned and of mitochondrial origin. NJ method assumes equal evolutionary rates across lineages (no molecular clock). The raw distance model assumes equal mutation weights and does not correct for multiple substitutions.
- Limitations: File size (215 KB) made direct computation slow; splitting into parts was necessary. Less accurate than model-based methods (e.g., Maximum Likelihood or Bayesian). If sequence quality or alignment errors are present, trees may reflect noise instead of true evolutionary relationships.
- Strengths: NJ is computationally inexpensive and appropriate for large alignments. Splitting the dataset allowed analysis to proceed despite memory constraints. Side-by-side trees allow comparison of topology and potential clustering.

### Method 2: RAxML and R
I constructed a phylogenetic tree of woolly mammoth specimens using RAxML, a tool for maximum likelihood-based phylogenetic analysis. The tree was estimated using the GTR+G model, which accounts for nucleotide substitution and gamma-distributed rate heterogeneity across sites. The analysis was performed on the full dataset of aligned sequences.

RAxML in Terminal:
The sequence data was first modified to ensure unique sequence names and then processed using RAxML in the terminal. The following command was used to generate the tree:
```
/usr/local/Cellar/raxml/8.2.12/bin/raxmlHPC-AVX -s /Users/graceswenson/Desktop/563/modified_woolly_mammoth_aligned_final_unique.fasta -n my_tree_run3 -m GTRGAMMA -p 12345
```

Code: R Visualization of Phylogenetic Tree
Once the tree was generated by RAxML, it was visualized in R using the ape package. The following steps were taken:
```
install.packages("ape") 
library(ape) tree <- read.tree("/Users/graceswenson/Desktop/563/RAxML_bestTree.my_tree_run3")
plot(tree)
```

Tree saved as: "woolly_mammoth_raxml_tree.pdf"

- Assumptions: The input sequences are aligned and of mitochondrial origin. The GTR+G model was selected to account for nucleotide substitutions and site-rate heterogeneity. RAxML uses maximum likelihood estimation, which provides accurate tree topologies for large datasets.
- Limitations: The quality of the alignment is crucial. Any misalignment or sequence errors can affect the results. The RAxML analysis assumes no molecular clock, meaning it does not assume constant rates of evolution across all branches. Computationally intensive for large datasets, but the model offers a good balance of accuracy and computational feasibility.
- Strengths: RAxML provides robust maximum likelihood trees, widely used for large-scale phylogenetic analysis. The GTR+G model accounts for rate heterogeneity, making it suitable for diverse evolutionary scenarios. Visualizing the tree in R helps assess the quality and topology of the tree and provides a straightforward interpretation.

## Bayesian Inference
In this analysis, I used MrBayes to perform Bayesian inference of phylogeny for a set of mitochondrial genome sequences. This method utilizes Markov Chain Monte Carlo (MCMC) to estimate the most probable tree topologies and branch lengths, along with other model parameters. The analysis includes model selection, parameter settings, and interpretation of the resulting phylogenetic tree.

Step 1: Convert your aligned FASTA file (wolly_mammoth_aligned_final.fasta) into a NEXUS format using the following R code
```
library(ape)
fasta_file <- "wolly_mammoth_aligned_final.fasta"
alignment <- read.dna(fasta_file, format = "fasta")
nexus_file <- "wolly_mammoth_aligned_final.nex"
write.nexus(alignment, file = nexus_file)
```
Step 2: The following command block was used in the MrBayes analysis. It includes settings for the model and MCMC parameters, as well as a specified outgroup.
```
begin mrbayes;
    set autoclose=yes;
    prset brlenspr=unconstrained:exp(10.0);
    prset shapepr=exp(1.0);
    prset tratiopr=beta(1.0,1.0);
    prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
    lset nst=2 rates=gamma ngammacat=4;
    mcmcp ngen=1000000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
    outgroup EU153450.1;  # Change to a valid taxon if needed
    mcmc;
    sumt;
end;
```
To run MrBayes, open the terminal and navigate to your project directory where the NEXUS file is located. Use the following command:
```
mb wolly_mammoth_aligned_final.nex
```
Step 3: Now that you have the con.tre file, you can visualize the tree using IcyTree. Open IcyTree then upload the con.tre file.

Tree saved as: "woolly_mammoth_bayesian.jpeg"

Assumptions: Assumes that the evolutionary model you select (e.g., substitution rates, rate variation among sites) accurately reflects the underlying biological processes. Assumes sites in the sequence alignment evolve independently and identically according to the chosen model (though it can accommodate variation with models like Gamma). Assumes that the priors you specify for model parameters are appropriate—poor choices can lead to misleading results.
Strengths: Uses Markov Chain Monte Carlo (MCMC) to sample from the posterior distribution, providing posterior probabilities for tree topology, branch lengths, and model parameters. Directly estimates uncertainty in tree topology, branch lengths, and model parameters. 
Limitations: Poor or uninformed prior choices can significantly affect results. MCMC methods can be computationally intensive, especially for large datasets. Posterior probabilities are often higher than bootstrap values from ML, possibly overstating support.

## Coalescent method: BUCKy


## Software and Versions Used
- R (v4.3.1): Installed from [CRAN](https://cran.r-project.org/)
- ape package (v5.7-1): Installed via `install.packages("ape")`
- ClustalW (v2.1): Installed via Homebrew with `brew install clustal-w`
- RAxML (v8.2.12): Installed via Homebrew with `brew install brewsci/bio/raxml`
- MrBayes(v3.2.7a): Installed via Homebrew with `brew tap brewsci/bio` and `brew install mrbayes open-mpi`
- macOS Ventura 13.4 (M1 Chip)


