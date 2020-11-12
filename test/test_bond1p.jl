
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
r0 = rnn(:Al)
B = ACE.Utils.ace_basis( species = [:X, :Al], N = 1,
                         pin = 0, pcut = 0)

# define the cylindrical cut-off function
# Fcut = ACEtb.Bonds.BondCutoff(r0 = rnn(:Al))
Fcut = ACEtb.Bonds.BondCutoff(r0 = nothing, rcut = 2.5*r0, renv = 1.9*r0, zenv=1.9*r0, pcut = 2)

# convert the radial ACE basis into a cylindrical ACE basis.
Bbonds = ACEtb.Bonds.RPIBonds(B, Fcut)

#---



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

#---
# Plotting bond basis values
xgr = range(-2*r0, 3*r0, length = 40)
ygr = range(-2*r0, 2*r0, length = 40)
R0 = SVector([r0, 0.0, 0.0]...)
Zs = AtomicNumber.([:X, :Al])
z0 = AtomicNumber(:X)
X = []
Y = []
V = []

for x in xgr, y in ygr
   R1 = SVector([x, y, 0.0]...)
   # if norm(R1) < 0.5*r0 || norm(R1 - R0) < 0.5*r0
   #    continue
   # end
   Rs = [R0, R1]
   b0 = eval_bond(Bbonds, Rs, Zs, z0)
   # b0 = evaluate(Bbonds, Rs, Zs, z0)
   push!(V,b0)
   push!(X,x)
   push!(Y,y)
end

vplot = [V[k][12] for k in 1:length(V)]
using Plots
scatter(X,Y, marker_z = vplot, color = :jet, clims = (-0.7, 0.7))
plot!([0.0, r0], [0.0, 0.0], lw=4, m=:o, ms=6, c=:red, label = "")

# #---
# minimum(abs.(vplot))
# maximum(abs.(vplot))
# extrema(vplot)
#
# v[1]
# # display([ v[1] v[2] ][1:16,:])
