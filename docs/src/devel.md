
# Developer Notes

This section is dedicated to developing notes on the integration of ACEtb model to other codes such as DFTB+ as well as give detailed description of the main functions that is used in the calculations.

## Slater-Koster Functions

ACEtb package includes functions to efficiently generate and use Slater-Koster bond hopping parameters. The main functions are `sk2cart` and `cart2sk` that gather SK parameters for the given orbitals (`s`, `p`, `d` ...) and the corresponding angular momentums as well as the overlapping distances of orbitals $R$ and rotation angle, $\theta$.

To effectively store the orbital definitions, package uses `SK_Bond` definitions, which stores the orbitals in the given order and provides the combination of orbital-orbital bonds and the type of bonding such as $\sigma$, $\pi$ and so on.

Here for the example case, we select Al. A full orbital definition of Al up to `3d` orbital can be given as

$1s^2$, $2s^2$, $2p^6$, $3s^2$, $3p^1$ and $3d^0$

In above definition, the superscript numbers indicate the number of electrons in those orbitals. As ACEtb Slater-Koster module can hold any given number of orbitals for the given type of element, we can easily define the orbitals above in any order. Here, we will exactl use FHI-aims orbital definitions for the example. In FHI-aims, Tier 1 level Al is defined with `s`, `s`, `s`, `p`, `p`, `d`. These orbitals correspond to the 1s, 2s, 3s, 2p, 3p, and 3d orbitals.

In ACEtb package using `SlaterKoster` module one can define the sequence as follows:

```
orbitals = ["s", "s", "s", "p", "p", "d"]
SKorbs = SKH([SKOrbitals(o) for o in orbitals])
```

In above definition, two functions are 


