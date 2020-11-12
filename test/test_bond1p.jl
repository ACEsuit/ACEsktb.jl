
using ACE, LinearAlgebra, StaticArrays
using ACE: evaluate
using ACEtb
using Random: shuffle

# this should possibly go into the main ACEtb module, or possibly into
# the hamiltonian assembly code
function eval_bond(B, Rs, Zs, z0)
   r̂ = Rs[1] / norm(Rs[1])
   o = Rs[1]/2
   Rr = map( r_ -> (r = r_ - o; o + r - 2*dot(r,r̂)*r̂), Rs )
   Rr[1] = Rs[1]
   return evaluate(B, Rs, Zs, z0) + evaluate(B, Rr, Zs, z0)
end

# utility function for generating random isometries
function randsym(Rs)
   Rs1 = [ [Rs[1]]; shuffle(Rs[2:end]) ]
   K = randn(3, 3)
   K = K - K'
   Q = SMatrix{3,3}(rand([-1,1]) * exp(K)...)
   return [ Q * R for R in Rs1 ]
end

#---

# to create a basis with cylindrical symmetry start with a
# standard ACE basis
B = ACE.Utils.ace_basis( species = [:X, :Al], N = 1,
                         pin = 0, pcut = 0)

# define the cylindrical cut-off function
# Fcut = ACEtb.Bonds.BondCutoff(r0 = rnn(:Al))
Fcut = ACEtb.Bonds.BondCutoff(r0 = nothing, rcut = 5., renv = 3., zenv=2., pcut = 2)

# convert the radial ACE basis into a cylindrical ACE basis.
Bbonds = ACEtb.Bonds.RPIBonds(B, Fcut)


#---

# Testing symmetry
for ntest = 1:10
   Rs, Zs, z0 = ACE.Random.rand_nhd(10, B.pibasis.basis1p.J, :Al)
   Zs[1] = 0
   b0 = eval_bond(Bbonds, Rs, Zs, z0)
   Rs1 = randsym(Rs)
   b1 = eval_bond(Bbonds, Rs1, Zs, z0)
   @show b1 ≈ b0
end

# Plotting bond basis values
x = []
y = []
v = []
for ntest = 1:1000
   Rs, Zs, z0 = ACE.Random.rand_nhd(2, B.pibasis.basis1p.J, :Al)
   # z0 = AtomicNumber(:X)
   a = Rs[2][1]
   b = Rs[2][2]
   Rs[2] = SVector([a b 0.]...)
   Rs[1] = SVector([1. 0. 0.]...)
   Zs[1] = 0
   b0 = eval_bond(Bbonds, Rs, Zs, z0)
   push!(v,b0)
   push!(x,a)
   push!(y,b)
end

v
vplot = [v[k][25] for k in 1:length(v)]
using Plots
scatter(x,y, marker_z = vplot, color = :jet)

minimum(abs.(vplot))
maximum(abs.(vplot))
