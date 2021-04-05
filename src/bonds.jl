
module Bonds


using ACE, JuLIP, NeighbourLists
using LinearAlgebra: norm, dot
using StaticArrays

using JuLIP.Potentials: neigsz!
using JuLIP: evaluate

import JuLIP: cutoff, write_dict, read_dict, evaluate!
import JuLIP.Potentials: zlist, z2i, alloc_temp, numz, z2i, i2z
import ACE: get_basis_spec, fltype, rfltype, OneParticleBasis, add_into_A!
import Base: ==, length
using ACE: z2i, i2z


# ------------------------------------------------------
# Auxiliary Objects: fcut

@doc raw"""
`BondCutoff: ` This defines the a sensible first guess at a
cutoff mechanism for environment-dependent pair bonds. The
envelope for the bond-length is simply
```math
   f_{\rm cut}(r) = (r - r_{\rm cut})^p
```
with p normally equal to 2.

The envelope for an environment bond is
less obvious to define. Here we make the choice to give a constant
radius into the z and r directions. Suppose that in cylindrical coordinates
there is an atom at position $(z, r, \theta)$ relative to the mid-point of the
bond, then the envelope will be
```math
   f_{\rm env}(r, z) =
      (z^2 - z_{\rm eff}^2)^p
      \cdot
      (r^2 - r_{\rm env}^2)^p
```
where $z_{\rm eff} = r_{\rm bond}/2 + z_{\rm env}$.
"""
struct BondCutoff
   pcut::Int
   # cutoff parameters
   rcut::Float64
   # envelope parameters
   renv::Float64
   zenv::Float64
end

BondCutoff(; r0 = nothing, pcut = 2,
             rcut = 3 * r0, renv = 2 * r0,  zenv = 2 * r0) =
   BondCutoff(pcut, rcut, renv, zenv)

# cutoff for the bond
fcut(cut::BondCutoff, R) = ((norm(R)/cut.rcut)^2 - 1)^cut.pcut*(norm(R)<=cut.rcut)

function fenv(cut::BondCutoff, R, R0)
   z, r = _get_zr(R, R0)
   # zeff = r/2 + cut.zenv
   zeff = norm(R0)/2 + cut.zenv
   return ((z/zeff)^2 - 1)^cut.pcut * (abs.(z)<=zeff) *
          ((r/cut.renv)^2 - 1)^cut.pcut * (r<=cut.renv)
end

function _get_zr(R, R0)
   R̂0 = R0/norm(R0)
   o = R0/2
   z = dot(R - o, R̂0)
   r = norm(R - o - z * R̂0)
   # r = norm(R-o)
   return z, r
end

function get_i_neigh_Rs(i, coords, nnei, inei)
   Rt = [] # holds neighbours of atom i
   offset = i == 1 ? 0 : sum(nnei[1:i-1])
   for nj = 1:nnei[i]
      jn = offset + i + nj
      j = inei[jn]
      Rij =  SVector((coords[:,j] - coords[:,i])...)
      push!(Rt,Rij)
   end
   return Rt
end

function get_all_neighs(natoms, coords, nnei, inei)
   #Rt = [ [] for ia = 1:natoms ] # holds neighs of all atoms
   Rt = []
   
   #Threads.@threads for ia = 1:natoms
   for ia = 1:natoms
      Rij = []
      offset = ia == 1 ? 0 : sum(nnei[1:ia-1])
      for nj = 1:nnei[ia]
         jn = offset + ia + nj
         ja = inei[jn]
         R0 =  SVector((coords[:,ja] - coords[:,ia])...)
         push!(Rij,R0)
      end
      push!(Rt,Rij)
   end
   return Rt
end

function get_env_neighs(Rt, R0, cut::BondCutoff)
   Renv = []
   # condition on the bond length
   rmax = sqrt((norm(R0)+abs(cut.zenv))^2 + (cut.renv)^2)
   if norm(R0) <= cut.rcut
      for R in Rt
         if (norm(R) < rmax) || norm(R-R0) > 1e-10
            # Get the length and radius of cylinder 
            #  that encloses i-j neighbours.
            z, r = _get_zr(R, R0) 
            if (z<= cut.zenv)&&(r<=cut.renv)
               push!(Renv,R)
            end
         end
      end
   end
   return Renv
end

function get_env(at::AbstractAtoms{T}, R0, i, cut::BondCutoff; nlist = nothing) where {T}
   Renv = []
   # condition on the bond length
   if norm(R0) <= cut.rcut
      if nlist == nothing
         rmax = sqrt((norm(R0)+abs(cut.zenv))^2 + (cut.renv)^2)
         nlist = neighbourlist(at, rmax)
      end
      maxR = maxneigs(nlist)
      tmpRZ = (R = zeros(JVec{T}, maxR), Z = zeros(AtomicNumber, maxR))
      j, Rt, Z = neigsz!(tmpRZ, nlist, at, i)
      for R in Rt
         if norm(R-R0) > 1e-10
            z, r = _get_zr(R, R0)
            if (z<= cut.zenv)&&(r<=cut.renv)
               push!(Renv,R)
            end
         end
      end
   end
   return Renv
end

function eval_bond(B, Rs, Zs, z0)
   r̂ = Rs[1] / norm(Rs[1])
   o = Rs[1]/2
   Rr = map( r_ -> (r = r_ - o; o + r - 2*dot(r,r̂)*r̂), Rs )
   Rr[1] = Rs[1]
   return evaluate(B, Rs, Zs, z0) + evaluate(B, Rr, Zs, z0)
end

function get_basis(order, degree, Fcut; Deg = nothing)
    # Create a basis with cylindrical symmetry start with a
    # standard ACE basis
    if Deg == nothing
       Dn = Dict( "default" => 1.,)
       Dl = Dict( "default" => 1.,)
       Dd = Dict( "default" => 1.,)
       Deg = ACE.RPI.SparsePSHDegreeM(Dn, Dl, Dd)
    end
    B = ACE.Utils.ace_basis( species = [:X, :Al], N = order,
                                 pin = 0, pcut = 0,
                                 maxdeg = degree,
                                 D = Deg)
    # convert the radial ACE basis into a cylindrical ACE basis.
    Bbonds = RPIBonds(B, Fcut)
    iX = z2i(B, AtomicNumber(:X))
    # b_index = B.pibasis.inner[iX].AAindices
    b_index = B.Bz0inds[1]
    return Bbonds, b_index
end


write_dict(basis::BondCutoff) = Dict(
         "__id__" => "ACEtb_BondCutoff",
            "pcut" => basis.pcut,
            "rcut" => basis.rcut,
            "renv" => basis.renv,
            "zenv" => basis.zenv )

read_dict(::Val{:ACEtb_BondCutoff}, D::Dict) =
   BondCutoff(D["pcut"], D["rcut"], D["renv"], D["zenv"])


# ------------ Main Object - Bond1pBasis
#  1-particle basis for bonds.

"""
`struct Bond1pBasis` : type of 1-p basis to be used for bonds rather than
site energies. Some functions will be overloaded to properly account for
the additional symmetries.

The basis must be built for the actual species occuring in the system of
interest as well as an additional artificial species `:X` which is used
to indicate the atom defining the bond.
"""
struct Bond1pBasis{T, TACE} <: OneParticleBasis{T}
   ace::TACE
   fcut::BondCutoff
   __T::Type{T}
end

# TODO: add sanity check that ace has the :X in the right place?!
Bond1pBasis(ace::OneParticleBasis{T}, fcut) where {T} =
   Bond1pBasis(ace, fcut, T)

# ------ book-keeping

zlist(basis::Bond1pBasis) = zlist(basis.ace)
cutoff(basis::Bond1pBasis) = cutoff(basis.ace)
length(basis::Bond1pBasis, args...) = length(basis.ace, args...)
get_basis_spec(basis::Bond1pBasis, args...) = get_basis_spec(basis.ace, args...)

==(P1::Bond1pBasis, P2::Bond1pBasis) = (P1.ace == P2.ace) && (P1.fcut == P2.fcut)

write_dict(basis::Bond1pBasis) = Dict(
         "__id__" => "ACEtb_Bond1pBasis",
            "ace" => write_dict(basis.ace),
           "fcut" => write_dict(basis.fcut) )

read_dict(::Val{:ACEtb_Bond1pBasis}, D::Dict) =
   Bond1pBasis(read_dict(D["ace"]), read_dict(D["fcut"]))


# ------------------------------------------------------
#  Evaluation code

fltype(basis::Bond1pBasis) = fltype(basis.ace)
rfltype(basis::Bond1pBasis) = rfltype(basis.ace)

maxlength(basis::Bond1pBasis) =
      maximum( length(basis, iz1, iz2)
               for (iz1, iz2) in Base.Iterators.product(1:numz(basis), 1:numz(basis)) )

alloc_temp(basis::Bond1pBasis, args...) =
   (
      Ptmp = zeros(fltype(basis), maxlength(basis)),
      tmpace = alloc_temp(basis.ace, args...)
   )


function evaluate!(A, tmp, basis::Bond1pBasis{TACE},
                   Rs, Zs::AbstractVector, z0) where {TACE}
   Rbond = Rs[1]
   fill!(A, 0)
   iz0 = z2i(basis, z0)
   # Should this be replaced by iz0 = z2i(basis, Zs[1]) ?
   P = tmp.Ptmp

   # center-bond
   @assert Zs[1] == 0
   fill!(P, 0)
   iz = z2i(basis, AtomicNumber(0))
   add_into_A!(P, tmp.tmpace, basis.ace, Rbond, iz, iz0)
   fc = fcut(basis.fcut, Rbond)
   Av = (@view A[basis.ace.Aindices[iz, iz0]])
   @. Av[:] = fc * P

   # environment
   for (R, Z) in zip(Rs[2:end], Zs[2:end])
      iz = z2i(basis, Z)
      fill!(P, 0)
      add_into_A!(P, tmp.tmpace, basis.ace, R, iz, iz0)
      Av = (@view A[basis.ace.Aindices[iz, iz0]])
      fenv_ = fenv(basis.fcut, R, Rbond)
      @. Av[:] += fenv_ * P
      # @. Av[:] += P
   end
   return A
end


# ------------------------------

function RPIBonds(B::RPIBasis, Fcut)
   basis1p = B.pibasis.basis1p
   bonds1p = Bond1pBasis(basis1p, Fcut)
   pibasis = PIBasis(bonds1p, B.pibasis.zlist, B.pibasis.inner, B.pibasis.evaluator)
   Bbonds = RPIBasis(pibasis, B.A2Bmaps, B.Bz0inds)
   return Bbonds
end


end
