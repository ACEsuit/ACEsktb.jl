
# ACEtb.jl Documentation


This package extends `ACE.jl` which provides approximation schemes for permutation and isometry invariant functions,
based on symmetric polynomial, in the case where the function to approximate satisfy permutation invariance with respect to like atoms, and rotation invariance with respect to a bond.

To ensure the symmetry, the basis functions are of the form

$f(R_{ij}; R_{env}) = g(R_{env}-R_i) + g(R_{env}-R_j)$,

where $g$ is an ACE function, i.e. satisfies rotation invariance. Here $R_{ij}$ is the vector of the bond and $R_{env}$ contains the coordinates of the atoms near the bond.

## Developers

 * Geneviève Dusson
 * Berk Onat
 * Reinhard Maurer
 * Christoph Ortner
 * James R. Kermode

## Citation

http://arxiv.org/abs/2111.13736

> Zhang, L., B. Onat, G. Dusson, G. Anand, R. J. Maurer, C. Ortner, and J. R. Kermode, 2021, Equivariant analytical mapping of first principles Hamiltonians to accurate and transferable materials models: ArXiv:2111.13736 [Cond-Mat].

## License

MIT
