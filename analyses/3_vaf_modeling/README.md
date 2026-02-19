## Evolutionary Modeling of TERTp variants

This directory contains code, data, and plotting notebooks for estimating the fitness effects of TERTp variants in context of other drivers of CH. 

To run any of the code, we suggest creating a conda environment: 

```
mamba env create -f env.yaml
conda activate vaf_modeling
```

The directories are organized as: 

### `code`

Contains code used for likelihoods under the neutral Wright-Fisher model and under the model from Macgregor et al 2025. These are two models for modeling variant-allele-fraction in bulk sequencing data 

### `data` 

Contains all relevant summary data for recreating plots pertaining to TERTp VAF 

### `notebooks`

Contains notebooks used to generate supplementary plots.

## Contact

* Arjun Biddanda
* Michael Tassia
* Rajiv C. McCoy
