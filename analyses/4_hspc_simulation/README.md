# Simulations of TERTp evolution in a population of hematopoietic stem/progenitor cells (HSPCs) over a human lifetime

## Scripts in this repo

* `00_make_popTL_distribution.R`: R script for generating 1000 uniformly distributed telomere lengths between the 1st and 99th population percentiles
* `01_simulation_wrapper.sh`: Wrapper script for submitting simulation jobs with different TERTp mutation rates
* `02_simulation_populationTL_plots.R`: Analysis and plotting code
* `hsc_simulation_v10_shortTL_vaf.slim`: SLiM code for simulation
* `run_simulations.sh`: SLURM array script for simulating HSPC population dynamics for one individual's lifetime

## Description of simulation

The simulation and this description is adapted from the simulation in [Deboy & Tassia et al., 2023](https://github.com/mccoy-lab/pot1_ch_2022/blob/main/simulations/README.md).

### Overview

We used simulations, constructed with the software SLiM ([Haller and Messer, 2019](https://academic.oup.com/mbe/article/36/3/632/5229931)), to test the role of telomere length (TL) and somatic mutation rate on the acquisition and allele frequency of a somatic mutation mimicking the effect of TERTp -124C>T.

We ran 1000 simulations, with each simulation modeling HSPCs population dynamics in one individual from birth to age 90. For each simulation, the individual's constitutional telomere length (TL) at birth was drawn from a uniform distribution from the 1st to 99th percentiles of clinically validated flow-FISH measured average TLs ([Alder et al., 2018](https://www.pnas.org/doi/10.1073/pnas.1720427115)):

- 1st percentile (8.6 kb)
- 99th percentile (13.4 kb)

We ran three groups of 1000 simulations, where HSPCs in each group had a different probability of acquiring the TERTp mutation with each cell division:

- `1e-04` mutation probability per cell per division
- `1e-06` mutation probability per cell per division
- `1e-08` mutation probability per cell per division

The lowest simulated mutation rate of `1e-08` mutations per cell per division is equivalent to ~40 mutations HSPC per year, in line with upper estimates of 22 mutations per HSPC per year reported by [Osorio et al.](https://www.cell.com/cell-reports/fulltext/S2211-1247(18)31760-1) and 16-32 by [Mitchell et al.](https://www.nature.com/articles/s41586-022-04786-y)

### Setup

We begin each simulation with a population of 100,000 cells, none of which carry the TERTp mutation.

We draw every cell’s telomere length from a normal distribution with mean = the individual's assigned average telomere length at birth and standard deviation = 1.25 kb. This standard deviation was determined from the distribution of telomere lengths within an individual at birth, as estimated in [Mitchell et al.](https://www.nature.com/articles/s41586-022-04786-y)

### Events

Within each day of the simulation, the same set of events occur:

#### 1. Cell fate decision

##### Division vs. non division

We set the HSPC division rate at once every 280 days (daily division probability = 1/280), as described in [Mitchell et al.](https://www.nature.com/articles/s41586-022-04786-y)

##### Fate decision

Probabilities of fate decisions were derived from [Lahouel et al.](https://www.pnas.org/doi/10.1073/pnas.1914589117) as follows:

* 0.9 asymmetric division
* 0.0525 symmetric division
* 0.0475 differentiation (i.e., dropout from HSPC population)

#### 2. Introduction of TERTp mutation

All _dividing_ cells have a probability of acquiring the TERTp mutation, equivalent to the mutation probability per cell per division that we set for each simulation. Cells with the TERTp mutation are exempt from telomere shortening (see next section).

We did not consider the effects of a mutation at any genomic site except TERTp in our simulation. Any additional mutation "hits" at the TERTp site beyond the first mutation had no additional effects.

#### 3. Telomere shortening

All dividing cells experience a 100 bp/division telomere shortening rate, estimated from [Mitchell et al.](https://www.nature.com/articles/s41586-022-04786-y).

Cells with the TERTp mutation are exempt from telomere shortening, meaning they are effectively immune from the cell death mechanism described in the next section. _Note that cells carrying TERTp can still stochastically drop out of the HSPC population via the "differentiation" fate decision._

#### 4. Cell death caused by short telomeres

Every cell with TL <4kb experiences cell death with a probability dependent entirely on telomere length. Likelihood of death increased linearly from 0 (certain survival) to 1 (certain death) between 4 kb and 2 kb ([Alder et al., 2018](https://www.pnas.org/doi/10.1073/pnas.1720427115)).

We note that dropout from the HSPC population in this simulation can occur due to either:

- **Telomere shortening** (somewhat deterministic - although the cell division decision that leads to telomere shortening, the exact amount of telomere shortening, and the probability of death at TL<4kb are all stochastic)
- **A "differentiation" cell fate decision** that removes the cell from the HSPC pool (fully stochastic)