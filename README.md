[![tests](https://github.com/ACEsuit/ACEtb.jl/actions/workflows/Tests.yml/badge.svg)](https://github.com/ACEsuit/ACEtb.jl/actions/workflows/Tests.yml)
[![documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://acesuit.github.io/ACEtb.jl//dev)


# ACEtb

A data-driven scheme to construct predictive models of Hamiltonian and overlap matrices in atomic orbital representation from ab initio data as a function of local atomic and bond environments. The scheme goes beyond conventional tight binding descriptions as it represents the ab initio model to full order, rather than in two-centre or three-centre approximations. We achieve this by introducing an extension to the Atomic Cluster Expansion (ACE) descriptor that represents intraatomic onsite and interatomic offsite blocks of Hamiltonian and overlap matrices that transform equivariantly with respect to the full rotation group in 3 dimensions. 

# Installation

```
] registry add https://github.com/JuliaMolSim/MolSim.git
] add https://github.com/ACEsuit/ACEtb.jl.git
```

# Paper

http://arxiv.org/abs/2111.13736

> Zhang, L., B. Onat, G. Dusson, G. Anand, R. J. Maurer, C. Ortner, and J. R. Kermode, 2021, Equivariant analytical mapping of first principles Hamiltonians to accurate and transferable materials models: ArXiv:2111.13736 [Cond-Mat].
