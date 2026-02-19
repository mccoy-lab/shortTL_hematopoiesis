# Short telomere clonal hematopoiesis
This repository is an archive for the bioinformatic code required to reproduce the evolutionary modeling and clone-based WGS data analyses presented in *Schratz et al. (TBD)*.

---

### [analyses](analyses/)

#### [1_clone_wgs_processing](analyses/1_clone_wgs_processing)

Snakemake pipeline used for processing of clone-based WGS data to generate somatic phylogenies. An additional [`readme`](analyses/1_clone_wgs_processing/README.md) is available for a full description of the pipeline. 

#### [2_modeling_convergence](analyses/2_modeling_convergence)

Rmarkdown notebook detailing the methods used to estimate the probability of recurrent mutation (≥ 2 mutation events) in a population of hematopoietic progenitors given empirical estimates for $\mu$ in the HSC compartment.

#### [3_vaf_modeling](analyses/3_vaf_modeling)

Jupyter notebook for estimating the mutation rate and fitness effects of *TERTp* somatic variants using the model presented in [MacGregor et al. 2025](https://www.biorxiv.org/content/10.1101/2023.12.21.572706v2).

#### [4_stochastic_character_mapping](analyses/4_stochastic_character_mapping)

Rmarkdown notebook detailing the methods for performing stochastic character mapping to infer the evolutionary origin of a somatic mutation given a phylogenetic tree.