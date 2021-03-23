module Predictions

using LinearAlgebra, LowRankApprox, StaticArrays
using ACE, ACEtb
using ACEtb.Bonds: BondCutoff, get_env, eval_bond, get_basis
using ACEtb.Utils: get_data

export predict

function predict(filenames, cutoff_params, fit_params)
    rcut = cutoff_params["rcut"]
    renv = cutoff_params["renv"]
    zenv = cutoff_params["zenv"]
    pcut = cutoff_params["pcut"]
    degree = fit_params["degree"]
    order = fit_params["order"]
    cutfunc = BondCutoff(pcut, rcut, renv, zenv)
    data_train = get_data(filenames, cutfunc, get_env)
    BII = fit_BI(data_train, order, degree, cutfunc; test = nothing)
    return BII, cutfunc
end

# -------------------------------------------------------------------------
# LSQ functions
function fit_BI(train, order, degree, cutfunc; test = nothing)
   Bbonds, b_ind = get_basis(order, degree, cutfunc) #get the basis
   nbonds = length(train[1][3])
   A = zeros(ComplexF64, (length(train), length(b_ind)))
   y = zeros(ComplexF64, (length(train), nbonds))
   for (i, (R0, Renv, V)) in enumerate(train)
     Rs = [[R0]; Renv]
     Zs = [[AtomicNumber(:X)];[AtomicNumber(:Al) for _ = 1:length(Renv)]]
     z0 = AtomicNumber(:X)
     A[i, :] = eval_bond(Bbonds, Rs, Zs, z0)[b_ind]
     y[i, :] = V
   end
   c = qr(A) \ y

   function BIfunc(R0,Renv)
       Rs = [[R0]; Renv]
       Zs = [[AtomicNumber(:X)];[AtomicNumber(:Al) for _ = 1:length(Renv)]]
       z0 = AtomicNumber(:X)
      return [dot(c[:,i],eval_bond(Bbonds, Rs, Zs, z0)[b_ind]) for i in 1:nbonds]
   end
   return BIfunc
end

end # End of module
