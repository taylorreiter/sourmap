# sourmap

Building diffusion maps from sourmash signatures. 
A [recent paper](https://www.nature.com/articles/s41467-020-18695-z) demonstrated that using metabolic models, diffusion maps capture the metabolic niche space of 2,621 bacteria in RefSeq.
In this study, the metabolic models were constructed using CarveMe, a tool that build genome-specific reaction sets from BiGG models of biochemical reactions. 
This method requires assembly which is computationally intensive and usually fails for some fraction of metagenome or metatranscriptome reads. 
Therefore, we sought a technique to build diffusion maps from bacterial sequencing data that does not rely on assembly and annotation. 

Sourmash generates scaled MinHah sketches (*signatures*) for sequencing reads or assemblies using k-mer sequences. 
K-mers are words of length *k* in nucleotide sequences.
There are about as many k-mers in a sequence as there are nucleotides, meaning that computational operations on all k-mers can be very expensive.
Scaled MinHash sketching consistently downsamples k-mers to a representative subset, enabling fast comparisons between or against many samples. 
K-mers and scaled MinHash sketches can approximate taxonomic relationships between samples ([see here](https://msystems.asm.org/content/1/3/e00020-16)). 
A *k* size of 21 approximately captures genus-level relationships, *k* = 31 approximately captures species-level relationships, and *k* = 51 approximately captures strain-level relationships.
Given this relationship with taxonomy, we hypothesized that abundance-aware scaled MinHash sketches could be used to build a diffusion map that recapitulates taxononomic relationships between samples. 
This repository tests this idea. 
It contains implementations of diffusion maps in R and in python (although the python implementation has not be properly tested and so may not currently work properly).
While we first plan to test diffusion maps on a small set of assembled genomes (GTDB, RefSeq), if we are successful in this space and especially if diffusion maps capture not just taxonomic relationships but also functional relationships, we plan to scale this analysis to metagenomes in the Sequence Read Archive. 
Given the scale of sequencing data in the SRA, we aim to build diffusion maps in Rust to better control memory usage and performance.
