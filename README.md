# sourmap

Building diffusion maps from sourmash signatures. 
A [recent paper](https://www.nature.com/articles/s41467-020-18695-z) demonstrated that using metabolic models, diffusion maps capture the metabolic niche space of 2,621 bacteria in RefSeq.
In this study, the metabolic models were constructed using CarveMe, a tool that build genome-specific reaction sets from BiGG models of biochemical reactions. 
This method requires assembly which is computationally intensive and usually fails for some fraction of metagenome or metatranscriptome reads. 
Therefore, we sought a technique to build diffusion maps from bacterial sequencing data that does not rely on assembly and annotation. 

Sourmash generates scaled MinHah sketches (*signatures*) for sequencing reads or assemblies using k-mer sequences. 
K-mers are words of length *k* in nucleotide sequences.
There are about as many k-mers in a sequence as there are nucleotides.
With such a large number of k-mers, computational operations on all k-mers can be very expensive -- e.g. all-by-all comparisons ([see here](https://peerj.com/articles/cs-94/)). 
Scaled MinHash sketching consistently downsamples k-mers to a representative subset, enabling fast comparisons between or against many samples. 
K-mers and scaled MinHash sketches can approximate taxonomic relationships between samples ([see here](https://msystems.asm.org/content/1/3/e00020-16)). 
A *k* size of 21 approximately captures genus-level relationships, *k* = 31 approximately captures species-level relationships, and *k* = 51 approximately captures strain-level relationships.
Given this relationship with taxonomy, we hypothesized that abundance-aware scaled MinHash sketches could be used to build a diffusion map that recapitulates taxononomic relationships between samples. 
This repository tests this idea. 
It contains implementations of diffusion maps in R and in python (although the python implementation has not be properly tested and so may not currently work properly).
While we first plan to test diffusion maps on a small set of assembled genomes (GTDB, RefSeq), if we are successful in this space and especially if diffusion maps capture not just taxonomic relationships but also functional relationships, we plan to scale this analysis to metagenomes in the Sequence Read Archive. 
Given the scale of sequencing data in the SRA, we aim to build diffusion maps in Rust to better control memory usage and performance.

@taylorreiter @luizirber Nov 2020

## how to run sourmap

One day sourmap will be a snakemake pipeline, but not today! 
Run instructions are below.
From a set of signatures that live in a directory and contain a certain string, build an abundance csv. 
Read this csv into R and transform to samples x features. 
Build diffusion map. 

Install dependencies

```
conda create -n sourmap sourmash rust r-magrittr r-readr
conda activate sourmap
git clone https://github.com/taylorreiter/sourmap.git
```

Build the abundance csv. 
Note that you'll need to replace the path to the signatures with whereever sourmash signatures live on your system.
Also note that you must currently supply the `scaled` value that is in the signatures (see [issue #5](https://github.com/taylorreiter/sourmap/issues/5) for what the error message looks like when you supply the wrong scaled value).

```
cd sourmap/sourmap_rs
cargo run --release -- abundance-matrix -k 21 --scaled 2000 --from-file -o abundance.csv <(find ~/work/sourmash-bio/greyhound/data/gtdb-r95 -type f -iname "*GCA_005*")
cd ..
```

In R, run the code below.
(note that if you're on a mac and have RStudio installed you can use `open -na RStudio` from within your conda env to launch an R session within your conda env)

```
library(magrittr)
library(readr)
source("build_dm.R")

abund <- read_csv("sourmap_rs/abundance.csv", col_names = F)
abund <- t(abund)  # this is only for the current implementation. Will be removed later.
dm <- build_dm(abund)
```



