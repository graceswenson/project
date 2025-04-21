# Siberian Wolly Mammoth Phylogenetic Analysis
- Name of Article: Intraspecific phylogenetic analysis of Siberian woolly mammoths using complete mitochondrial genomes
- DOI:  https://doi.org/10.1073/pnas.0802315105
- Link to Article: https://www.pnas.org/doi/10.1073/pnas.0802315105#supplementary-materials

## Species Information
The Siberian Wolly Mammoth is an extinct species of elephant that roamed parts of the Northern Hemisphere during the Pleistocene epoch. Their remains have been remarkably well-preserved in the permafrost of Siberia, allowing for easier study of their DNA. In this study, mitochondrial DNA (mtDNA) data sets collected from hair shafts will be used to complete a phylogenetic analysis. Samples from __________ species will be used.

## Data Source
I retrieved complete mitochondrial genome sequences for various Siberian woolly mammoth (Mammuthus primigenius) specimens from the NCBI GenBank database. Each genome was downloaded individually in FASTA format based on accession numbers provided in the published article. The sequences were then compiled into a single raw file and saved as wolly_mammoth_data_raw.fasta.

## ALignment Process
To align the complete mitochondrial genomes, I used the multiple sequence alignment tool ClustalW in the terminal with the following command:
clustalw -ALIGN -INFILE=wolly_mammoth_data_raw.fasta -OUTFILE=wolly_mammoth_aligned_final.fasta -OUTPUT=FASTA

Alignment Score: 2983150
Output File: wolly_mammoth_aligned_final.fasta

ClustalW
- Description: ClustalW is a widely used progressive alignment algorithm that aligns biological sequences by first conducting pairwise alignments, constructing a guide tree based on sequence similarity, and then progressively aligning sequences according to the guide tree.
- Assumptions: Sequences with higher similarity are more closely related evolutionarily. Evolutionary events like substitutions, insertions, and deletions occur independently. Gap penalties approximate biological realities.
- Limitations: Errors introduced in early pairwise alignments can propagate through the final alignment. Predefined gap penalties may not always reflect true evolutionary constraints. For large datasets, ClustalW is generally slower and less accurate compared to newer tools like MUSCLE.
- Source for ClustalW: https://doi.org/10.1093/nar/22.22.4673




