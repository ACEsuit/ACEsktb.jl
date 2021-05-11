
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

## Generating Hamiltonian `H` and Overlap Matrix `S` in ACEtb.

The main function that is used to generate Hamiltonian and overlap is `buildHS`.

```
buildHS(SKH_list, H, S, iatf, iatl, natoms, pos, species, nneigh, ineigh, ipair, i2a, norbe, onsite_terms, Bondint_table, cutoff_func, cutoff; MPIproc=MPIproc)
```

The arguements of the functions are defined as follows:

- SKH_list : SlaterKoster orbital list such as s-s, s-p, p-p but also includes sigma and pi orbital definitions in the same order as in H and S.
- H : 1D flattened Hamiltonian matrix
- S : 1D flattened overlap matrix
- iatf : index of first atom (Ex.:iatf=1)
- iatl : index of last atom (Ex.:iatl=729)
- natoms : Total number of atoms in configuration
- pos : positions of atoms (729 x 3)
- species : element type list for each atom in pos.
- nneigh : the array of the number of neighbours for each atom "i" (1D array).
- ineigh : the array of the neighbouring atom ids ("j") for each atom "i" (1D array).
- ipair : the index where we need to write the Hamiltonian values that correspond to an atom pair.
- i2a : index of atom in pos, species array to actual atom element type.
- norbe : number of orbitals for each element type.
- onsite_terms : the onsite block's diagonal hamiltonian elements.
- Bondint_table : The function for ACEtb model that holds bond integral predictions. 
- cutoff_func : The cut off function for ACEtb model.
- cutoff : The cut off distance in Angstrom.
- MPIproc : If it is 1, the running process is the leading process in DFTB+ MPI parallelisation.


