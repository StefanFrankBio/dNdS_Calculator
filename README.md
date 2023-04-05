# dNdS Calculator

This repository contains a Python script that calculates the dN/dS ratio for multiple sequence alignments of coding DNA sequences. The dN/dS ratio is an important measure in molecular evolution, comparing the rate of non-synonymous substitutions (dN) to the rate of synonymous substitutions (dS) in coding DNA sequences. A dN/dS ratio > 1 indicates positive selection, a ratio < 1 indicates purifying (negative) selection, and a ratio ≈ 1 suggests neutral evolution.

## Features  

The script can be run in different modes to:

- Calculate site counts for non-synonymous and synonymous sites in codons.
- Calculate substitution counts for non-synonymous and synonymous substitutions between codon pairs.
- Transpose a multiple sequence alignment, creating a JSON file with counts of each codon per site.
- Generate a consensus sequence from the transposed alignment.
- Calculate dN/dS ratios per sequence in the alignment compared to the consensus sequence.
- Calculate dN/dS ratios per site in the alignment.

## Command Line Usage
You can run the script using different commands and options:
### Site Counts
Calculate the non-synonymous and synonymous site counts for each codon:
```bash
./dNdS_calculator.py site_counts -c CODON_TABLE -o OUTFILE
```
### Substitution Counts
Calculate the non-synonymous and synonymous substitution counts for each pair of codons:
```bash
./dNdS_calculator.py sub_counts -c CODON_TABLE -o OUTFILE
```
### Transpose
Transpose a multiple sequence alignment, creating a JSON file with counts of each codon per site:
```bash
./dNdS_calculator.py transpose -i INFILE -i PREFIX
```
### Consensus
Generate a consensus sequence from the transposed alignment:
```bash
./dNdS_calculator.py consensus -i INFILE -p PREFIX
```
### Per Sequence dN/dS
Calculate dN/dS ratios per sequence in the alignment compared to the consensus sequence:
```bash
./dNdS_calculator.py per_sequence -i INFILE -p PREFIX -s SITE_COUNTS -u SUB_COUNTS
```
### Per Site dN/dS
Calculate dN/dS ratios per site in the alignment:
```bash
./dNdS_calculator.py per_site -i INFILE -p PREFIX -s SITE_COUNTS -u SUB_COUNTS
```
## Future Work
Potential future additions to this repository include:
1. Adding support for different genetic codes, such as mitochondrial or alternative nuclear codes.
2. Extending the script to accommodate for gaps and ambiguous characters in the alignment.
3. Implementing sliding window analysis for calculating dN/dS ratios over a specified window size.
4. Providing visualization options for the generated dN/dS data, such as heatmaps or line plots.
5. Allowing for parallel processing to speed up calculations on large datasets.



