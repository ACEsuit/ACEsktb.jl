module Predictions

using LinearAlgebra, LowRankApprox, StaticArrays
using ACE, ACEtb
using ACEtb.Bonds: BondCutoff, get_env, eval_bond, get_basis, read_dict, write_dict
using ACEtb.Utils: get_data, read_json, write_json
using JSON
using IterativeSolvers

export predict, train_and_predict

# Define the dictionaries from the parameters
function degreeM(deg,order;env_deg = deg)
   zX = AtomicNumber(:X)
   if order == 1
       Dn = Dict( "default" => deg,
                    (zX,zX) => 1.
                 )
       Dl = Dict( "default" => 1., )
       Dd = Dict( "default" => 1,
                  0 => eps(),
                  1 => deg,
               )
       return ACE.RPI.SparsePSHDegreeM(Dn, Dl, Dd)
  else
      Dn = Dict( "default" => env_deg,
                   (zX,zX) => 1.,
                )
      Dl = Dict( "default" => 1., )
      Dd = Dict( "default" => 1,
                 0 => eps(),
                 1 => deg,
                 2 => deg,
                 3 => deg,
                 4 => deg,
              )
      return ACE.RPI.SparsePSHDegreeM(Dn, Dl, Dd)
  end
end

function train_and_predict(filenames, specie_syms, cutoff_params, fit_params)
    rcut = cutoff_params["rcut"]
    renv = cutoff_params["renv"]
    zenv = cutoff_params["zenv"]
    pcut = cutoff_params["pcut"]
    degree = fit_params["degree"]
    order = fit_params["order"]
    env_deg = fit_params["env_deg"]
    cutfunc = BondCutoff(pcut, rcut, renv, zenv)
    @info "│    Getting data..."
    data_train = get_data(filenames, cutfunc, get_env)
    @info "│    fitting BI..."
    BII, train_dict = fit_BI(data_train, specie_syms, order, degree, env_deg, cutfunc; test = nothing)
    return BII, cutfunc, train_dict
end

function predict(poten_dict, cutoff_params, fit_params)
    rcut = cutoff_params["rcut"]
    renv = cutoff_params["renv"]
    zenv = cutoff_params["zenv"]
    pcut = cutoff_params["pcut"]
    order = fit_params["order"]
    degree = fit_params["degree"]
    env_deg = fit_params["env_deg"]
    cutfunc = BondCutoff(pcut, rcut, renv, zenv)
    BII = load_BI(poten_dict; test = nothing)
    return BII, cutfunc
end

function load_BI(poten_dict; test = nothing)
   basis_string = poten_dict["basis"]
   basis = read_dict(JSON.parse(basis_string))
   b_index = poten_dict["basis_index"]
   c = convert(Array{Float64,2}, poten_dict["c"])
   @info "c: ",c
   nbonds = poten_dict["nbonds"]
   specie_syms = poten_dict["elm_names"]

   @info "│    set BI func."
   function BIfunc(R0,Renv)
      Rs = [[R0]; Renv]
      Zs = [[AtomicNumber(:X)];[AtomicNumber(Symbol(specie_syms[1])) for _ = 1:length(Renv)]]
      z0 = AtomicNumber(:X)
      return [dot(c[:,i],eval_bond(basis, Rs, Zs, z0)[b_index]) for i in 1:nbonds]
   end
   @info "│    return BI func."
   return BIfunc
end

function fit_BI(train, specie_syms, order, degree, env_deg, cutfunc; test = nothing)
   @info "│    setting degreeM."
   Deg = degreeM(degree,order;env_deg = env_deg) 
   @info "│    setting basis and b_index."
   basis, b_index = get_basis(order, 1., cutfunc; Deg = Deg) 
   @info "│    basis functions set."
   nbonds = length(train[1][3])
   A = zeros(ComplexF64, (length(train), length(b_index)))
   y = zeros(ComplexF64, (length(train), nbonds))
   for (i, (R0, Renv, V)) in enumerate(train)
     Rs = [[R0]; Renv]
     Zs = [[AtomicNumber(:X)];[AtomicNumber(Symbol(specie_syms[1])) for _ = 1:length(Renv)]]
     z0 = AtomicNumber(:X)
     A[i, :] = eval_bond(basis, Rs, Zs, z0)[b_index]
     y[i, :] = V
   end
   @info "│    LSQ run."
   c = qr(A) \ y

   @info "│    set dict."
   train_dict = Dict()
   train_dict["size_A"] = size(A)
   train_dict["cond_A"] = cond(A)
   train_dict["rank_A"] = rank(A)
   train_dict["c"] = real.(c)
   train_dict["error_train"] = norm(A*c-y)
   train_dict["basis"] = JSON.json(write_dict(basis))
   train_dict["basis_index"] = b_index
   train_dict["nbonds"] = nbonds
   train_dict["elm_names"] = specie_syms
   
   @info "│    set BI func."
   function BIfunc(R0,Renv)
       Rs = [[R0]; Renv]
       Zs = [[AtomicNumber(:X)];[AtomicNumber(Symbol(specie_syms[1])) for _ = 1:length(Renv)]]
       z0 = AtomicNumber(:X)
      return [dot(c[:,i],eval_bond(basis, Rs, Zs, z0)[b_index]) for i in 1:nbonds]
   end
   @info "│    return BI funcs."
   return BIfunc, train_dict
end

end # End of module
