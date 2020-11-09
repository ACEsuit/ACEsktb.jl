
using ACE, LinearAlgebra, StaticArrays
using ACE: evaluate
using ACEtb

B = ACE.Utils.ace_basis( species = [:X, :Al], N = 4,
                         pin = 0, pcut = 0)

Fcut = ACEtb.Bonds.BondCutoff(r0 = rnn(:Al))
Bbonds = ACEtb.Bonds.RPIBonds(B, Fcut)

function eval_bond(B, Rs, Zs, z0)
   r̂ = Rs[1] / norm(Rs[1])
   o = Rs[1]/2
   Rr = map( r_ -> (r = r_ - o; o + r - 2*dot(r,r̂)*r̂), Rs )
   Rr[1] = Rs[1]
   return evaluate(B, Rs, Zs, z0) + evaluate(B, Rr, Zs, z0)
end

function randsym(Rs)
   Rs1 = [ [Rs[1]]; shuffle(Rs[2:end]) ]
   K = randn(3, 3)
   K = K - K'
   Q = SMatrix{3,3}(rand([-1,1]) * exp(K)...)
   return [ Q * R for R in Rs1 ]
end

Rs, Zs, z0 = ACE.Random.rand_nhd(10, B.pibasis.basis1p.J, :Al)
Zs[1] = 0
b0 = eval_bond(Bbonds, Rs, Zs, z0)


for ntest = 1:10
   Rs, Zs, z0 = ACE.Random.rand_nhd(10, B.pibasis.basis1p.J, :Al)
   Zs[1] = 0
   b0 = eval_bond(Bbonds, Rs, Zs, z0)
   Rs1 = randsym(Rs)
   b1 = eval_bond(Bbonds, Rs1, Zs, z0)
   @show b1 ≈ b0
end
