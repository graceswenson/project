# Bubo virginianus Phylogenetic Analysis

## Dataset Overview
This dataset includes genome-wide nuclear and mitochondrial DNA sequence data from 27 samples of *Bubo virginianus* (Great Horned Owl) across 13 of the 16 recognized subspecies. 
https://pmc.ncbi.nlm.nih.gov/articles/PMC10422955/

### Key Points
- **Species:** *Bubo virginianus*
- **Data Type:** Nuclear and mitochondrial DNA sequences
- **Sampling Locations:** 13 subspecies across North and South America

## Data Collection
- The dataset is curated from publicly available genomic data and is provided as part of a study focusing on owl species and their evolutionary history.
- Downloaded sequence data from NCBI.
- Retrieved mitochondrial DNA and UCE data files in FASTA format for analysis.

## UCE Sequence Data Processing
- Used the Phyluce pipeline. Clean reads assembled with Trinity and mapped with Phyluce mapping workflow. Final alignments were done using MAFFT.
Description:
- This script aligns UCE sequences using MAFFT, a multiple sequence alignment tool.
- MAFFT (Multiple Alignment using Fast Fourier Transform) is an efficient and widely used method for aligning large sequence datasets.
Assumptions:
1. Input sequences are homologous and suitable for multiple sequence alignment.
2. The dataset does not have excessive gaps or sequencing errors that could negatively impact alignment quality.
3. The computational resources are sufficient to handle the dataset.
Limitations:
1. MAFFT can struggle with highly divergent sequences, leading to misalignments.
2. The alignment may introduce artificial gaps, affecting downstream analyses.
3. Computational cost increases with larger datasets.

# Load MAFFT module (if running on an HPC system with module loading)
# module load mafft

# Define input and output files
INPUT_SEQ="uce_sequences.fasta"
OUTPUT_ALIGN="uce_aligned.fasta"

# Run MAFFT with auto mode for best strategy selection
mafft --auto $INPUT_SEQ > $OUTPUT_ALIGN

# Print completion message
echo "Alignment complete. Output saved to $OUTPUT_ALIGN"

# Algorithm Description:
# The Neighbor-Joining (NJ) method constructs a phylogenetic tree by iteratively
# finding the closest pair of operational taxonomic units (OTUs) and merging them.
# It minimizes the total branch length, making it computationally efficient.
#
# Assumptions:
# - Sequences evolve under a distance-based model.
# - The true evolutionary tree is approximately additive.
#
# Limitations:
# - NJ does not model character evolution explicitly.
# - It can be sensitive to noise and long-branch attraction.

# ---- Parsimony-based tree estimation ----
# Convert DNA alignment into a phyDat object (needed for parsimony analysis)
alignment_phyDat <- phyDat(alignment, type = "DNA")

# Generate an initial tree using the NJ method
init_tree <- nj(dist_matrix)

# Perform parsimony tree search using the Fitch algorithm
pars_tree <- optim.parsimony(init_tree, alignment_phyDat)

# Plot the Parsimony tree
plot(pars_tree, main = "Parsimony-based Tree")

# Save the parsimony tree to a file
write.tree(pars_tree, file = "parsimony_tree.nwk")

# Algorithm Description:
# The maximum parsimony method finds the tree that requires the fewest evolutionary
# changes to explain the observed sequences.
#
# Assumptions:
# - Evolution follows the principle of parsimony (minimizing changes).
# - All changes are equally probable.
#
# Limitations:
# - Can be computationally intensive for large datasets.
# - Sensitive to long-branch attraction.



https://www.pnas.org/doi/10.1073/pnas.0802315105#supplementary-materials

