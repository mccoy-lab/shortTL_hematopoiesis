# Short telomere clonal hematopoiesis
This repository is an archive for the bioinformatic code required to reproduce the evolutionary modeling and clone-based WGS data analyses presented in *Schratz et al. (TBD)*.

[![DOI](https://zenodo.org/badge/1135736159.svg)](https://doi.org/10.5281/zenodo.24139210)

---

### [analyses](analyses/)

#### [1_clone_wgs_processing](analyses/1_clone_wgs_processing)

Snakemake pipeline used for processing of clone-based WGS data to generate somatic phylogenies. An additional [`readme`](analyses/1_clone_wgs_processing/README.md) is available for a full description of the pipeline. 

#### [2_modeling_convergence](analyses/2_modeling_convergence)

Rmarkdown notebook detailing the methods used to estimate the probability of recurrent mutation (≥ 2 mutation events) in a population of hematopoietic progenitors given empirical estimates for $\mu$ in the HSC compartment.

#### [3_stochastic_character_mapping](analyses/3_stochastic_character_mapping)

Rmarkdown notebook detailing the methods for performing stochastic character mapping to infer the evolutionary origin of a somatic mutation given a phylogenetic tree.

#### [4_hspc_simulation](analyses/4_hspc_simulation)

Rmarkdown notebook describing simulations of TERTp somatic evolution in HSPC populations, and code for running simulations and plotting output.

#### [5_phylogenies](analyses/5_phylogenies)

Directories containing phylogenetic trees scaled by [molecular time](analyses/5_phylogenies/cellphy), [mutation burden](analyses/5_phylogenies/treemut), and [age](analyses/5_phylogenies/ultrametric).
