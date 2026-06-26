# Phylogenetic trees

This directory contains tree files for each of the phylogenies presented in this study. 
- [Maximum-likelihood trees](cellphy/) were inferred from somatic SNVs using [CellPhy v0.9.2](https://github.com/amkozlov/cellphy/tree/v0.9.2). 
- [Phylogenies with branch lengths scaled to mutation burden](treemut/) using [`treemut`](https://github.com/NickWilliamsSanger/treemut). An additional tree with branch lengths corrected for sensitivity are also provided; sensitivity correction was performed using the methods described in [*Mitchel et al. 2022*](https://www.nature.com/articles/s41586-022-04786-y)
- [Ultrametric trees](ultrametric/) were inferred from the sensitivity-correceted burden trees using the methods described in [*Mitchel et al. 2022*](https://www.nature.com/articles/s41586-022-04786-y)

All trees were rooted by setting each donor's respective matched normal as the outgroup. Matched normals were subsequently been pruned from the tree. 

The trees provided are all in newick format and are compatible with common web-based visualization tools like [iTOL](https://itol.embl.de/).