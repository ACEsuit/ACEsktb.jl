module Predictions

using LinearAlgebra, LowRankApprox, StaticArrays
using ACE, ACEtb
using ACEtb.Bonds: BondCutoff, get_env, eval_bond, get_basis, read_dict, write_dict
using ACEtb.Utils: get_data, read_json, write_json
using JSON
using IterativeSolvers
using LowRankApprox: pqrfact

using Infiltrator

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

function train_and_predict(filenames, specie_syms, cutoff_params, fit_params; MPIproc=1)
    rcut = cutoff_params["rcut"]
    renv = cutoff_params["renv"]
    zenv = cutoff_params["zenv"]
    pcut = cutoff_params["pcut"]
    degree = fit_params["degree"]
    order = fit_params["order"]
    env_deg = fit_params["env_deg"]
    rtol = get(fit_params, "rtol", 1e-15)
    method = Symbol(get(fit_params, "method", "QR"))
    cutfunc = BondCutoff(pcut, rcut, renv, zenv)
    if(MPIproc == 1)
       @info "│    Getting data..."
    end
    data_train = get_data(filenames, cutfunc, get_env)
    if(MPIproc == 1)
       @info "│    fitting BI..."
    end
    BII, train_dict = fit_BI(data_train, specie_syms, order, degree, env_deg, cutfunc; 
                             test = nothing, MPIproc=MPIproc, rtol=rtol, method=method)
    return BII, cutfunc, train_dict
end

function predict(poten_dict, cutoff_params, fit_params; MPIproc=1)
    rcut = cutoff_params["rcut"]
    renv = cutoff_params["renv"]
    zenv = cutoff_params["zenv"]
    pcut = cutoff_params["pcut"]
    order = fit_params["order"]
    degree = fit_params["degree"]
    env_deg = fit_params["env_deg"]
    cutfunc = BondCutoff(pcut, rcut, renv, zenv)
    BII = load_BI(poten_dict; test = nothing, MPIproc=MPIproc)
    return BII, cutfunc
end

function load_BI(poten_dict; test = nothing, MPIproc=1)
   basis_string = poten_dict["basis"]
   basis = read_dict(JSON.parse(basis_string))
   b_index = poten_dict["basis_index"]
   c = Array{Float64}(hcat(poten_dict["c"]...))
   nbonds = poten_dict["nbonds"]
   specie_syms = poten_dict["elm_names"]

   if(MPIproc == 1)
      @info "│    set BI func."
   end
   function BIfunc(R0,Renv)
      Rs = [[R0]; Renv]
      Zs = [[AtomicNumber(:X)];[AtomicNumber(Symbol(specie_syms[1])) for _ = 1:length(Renv)]]
      z0 = AtomicNumber(:X)
      B = eval_bond(basis, Rs, Zs, z0)[b_index]
      return [dot(c[:,i], B) for i in 1:nbonds]
   end
   if(MPIproc == 1)
      @info "│    return BI funcs."
   end
   return BIfunc
end

function fit_BI(train, specie_syms, order, degree, env_deg, cutfunc; 
                test = nothing, MPIproc=1, rtol=1e-15, method=:QR)
   if(MPIproc == 1)
      @info "│    setting degreeM."
   end
   Deg = degreeM(degree,order;env_deg = env_deg)
   if(MPIproc == 1)
      @info "│    setting basis and b_index."
   end
   basis, b_index = get_basis(order, 1., cutfunc; Deg = Deg)
   if(MPIproc == 1)
      @info "│    basis function is set."
   end
   nbonds = length(train[1][3])
   A = zeros((length(train), length(b_index)))
   y = zeros((length(train), nbonds))
   for (i, (R0, Renv, V)) in enumerate(train)
     Rs = [[R0]; Renv]
     Zs = [[AtomicNumber(:X)];[AtomicNumber(Symbol(specie_syms[1])) for _ = 1:length(Renv)]]
     z0 = AtomicNumber(:X)
     A[i, :] = eval_bond(basis, Rs, Zs, z0)[b_index]
     y[i, :] = V
   end
   train_dict = Dict()
   if method == :QR
      if(MPIproc == 1)
         @info "│    Solving LSQ with QR ..."
      end   
      qrA = qr(A)
      c = qrA \ y

   elseif method == :RRQR
      if(MPIproc == 1)
         @info "│    Solving LSQ with RRQR rtol=$rtol..."
      end   
      qrA = pqrfact(A, rtol=rtol)
      c = qrA \ y

   elseif method == :RegRRQR
      if(MPIproc == 1)
         @info "│    Solving LSQ with regularised RRQR rtol=$rtol..."
      end
      Γ = ACEtb.Bonds.scaling(basis, cutfunc.pcut)[b_index]
      Areg = A * Diagonal( 1 ./ Γ )
      qrAreg = pqrfact(Areg, rtol=rtol)
      c = Diagonal(1 ./ Γ) * (qrAreg \ y)

      train_dict["rank_qrAreg"] = qrAreg[:k]
      train_dict["abs_err_Areg"] = snormdiff(qrAreg, Areg)
      train_dict["rel_err_Areg"] = train_dict["abs_err_Areg"] / snorm(Areg)
   else
      error("unknown method $method")
   end
   if(MPIproc == 1)
      @info "│    set dict."
   end
   train_dict["size_A"] = size(A)
   train_dict["cond_A"] = cond(A)
   train_dict["rank_A"] = rank(A)
   train_dict["c"] = c
   train_dict["error_train"] = norm(A*c - y)
   train_dict["error_train_rel"] = train_dict["error_train"] / norm(y)
   train_dict["basis"] = JSON.json(write_dict(basis))
   train_dict["basis_index"] = b_index
   train_dict["nbonds"] = nbonds
   train_dict["elm_names"] = specie_syms

   if(MPIproc == 1)
      @info "│    setting BI function."
   end
   function BIfunc(R0,Renv)
      Rs = [[R0]; Renv]
      Zs = [[AtomicNumber(:X)];[AtomicNumber(Symbol(specie_syms[1])) for _ = 1:length(Renv)]]
      z0 = AtomicNumber(:X)
      B = eval_bond(basis, Rs, Zs, z0)[b_index]
      return [dot(c[:,i], B) for i in 1:nbonds]
   end
   if(MPIproc == 1)
      @info "│    return BI function."
   end
   return BIfunc, train_dict
end

end # End of module
