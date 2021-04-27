
# Fitting

The fitting of the bond integrals is done using a least-square fit (see the module `Predictions` - `predictions.jl`).

The fitting requires the choice of a basis, defined through body-order, maximum total degree, and environment degree, this last degree defining the contribution of the interactions between atoms that do not constitute the bond -see `degreeM`. The basis also depends on different cutoff parameters.

The data consists of configurations and corresponding bond integrals. From the configurations are computed the bonds and their corresponding neighbouring atoms (see `bonds.jl`).

The ACEtb basis functions are then evaluated on each bond+environment. This gives the LSQ matrix. The coefficients for the approximations of the bond integrals as linear combinations of ACEtb basis functions are then found solving a LSQ system (see `fit_BI`).
