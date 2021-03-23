module Bindings

using Pkg
using LinearAlgebra, LowRankApprox, Statistics, StaticArrays
using JuLIP
using ACE, ACEtb
using ACEtb.Bonds: BondCutoff, get_env, eval_bond, get_basis
using ACEtb.SlaterKoster
import ACEtb.SlaterKoster.CodeGeneration
using ACEtb.SlaterKoster: SKH, sk2cart, cart2sk, allbonds, nbonds
using ACEtb.Utils: get_data, read_json
using ACEtb.Predictions: predict
using ACEtb.TBhelpers
using ProgressMeter

bohr2ang = 0.188972598857892E+01

export set_model, model_predict, acetb_greetings

acetb_dct = Dict()

Bondint_table = nothing
cutoff_func = nothing
saveh5_HH = nothing
saveh5_SS = nothing
saveh5_satoms = nothing
model_onsite_vals = nothing

function buildHS(SKH_list, H, S, istart, iend, coords, species, nnei, inei, ipair, norbs, onsite_terms, atoms, Bondint_table, cutoff_func, cutoff; Bpredict=true, HH=nothing, SS=nothing, supercell_atoms=nothing, prt=0)

    if(prt == 1)
       prgres = Progress(iend-istart+1, dt=0.25, desc="AITB : Calculating ... ",
                         barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▅' ,'▆', '▇'],' ','|',),
                         barlen=20)
    end

    if Bpredict == false
        invcell = inv(atoms.cell)
        mesh = [9, 9, 9]
        rcn = get_tbcells(supercell_atoms, invcell, mesh, 365)
    end

    for ia = istart:iend
       isp = species[ia]
       offset = ia == 1 ? 0 : sum(nnei[1:ia-1])
       # Onsite blocks
       io = ipair[offset + ia] + 1
       if Bpredict == false
          nno = norbs[isp] * norbs[isp]
          shift = convert(Array{Int64,1}, [0,0,0])
          cid = find_row_in_matrix(shift, rcn)
          H[io : io + nno - 1] = vcat(HH[cid,:,:]...)
          S[io : io + nno - 1] = vcat(SS[cid,:,:]...)
       else
          for ib = 1:norbs[isp]
             H[io] = onsite_terms[isp][ib]
             S[io] = 1.0
             io += norbs[isp] + 1
          end
       end

       # Offsite blocks
       for nj = 1:nnei[ia]
          jn = offset + ia + nj
          ja = inei[jn]
          jsp = species[ja]
          ix = ipair[jn]
          iy = ix + norbs[isp] * norbs[jsp]
          ix += 1
          Rij =  transpose(coords[:,ja] - coords[:,ia])
          if cutoff < norm(Rij)
             continue
          end

          if Bpredict == false
              shift = convert(Array{Int64,1}, round.(invcell * Rij'))
              shift = wrap_shift(shift, mesh)
              cid = find_row_in_matrix(shift, rcn)
              H[ix : iy] = vcat(HH[cid,:,:]...)
              S[ix : iy] = vcat(SS[cid,:,:]...)
          else
              # Predictions 
              R0 = SVector(Rij...)
              Renv = get_env(atoms, R0, ia, cutoff_func)
              VV = Bondint_table(R0,Renv)

              # Set H and S
              lnb = length(SKH_list[isp].bonds)
              VH = VV[1:lnb]
              E  = sk2cart(SKH_list[isp], Rij, VH, FHIaims=true)
              H[ix : iy] = vcat(E...)
              VS = VV[lnb+1:end]
              ES = sk2cart(SKH_list[isp], Rij, VS, FHIaims=true)
              S[ix : iy] = vcat(ES...)
          end
       end
       if(prt == 1)
          next!(prgres)
       end
    end
    if(prt == 1)
       flush(stdout)
    end
end

function get_list_str(str::Array{UInt8}, ln, stride)
   slen = Int64(ln)
   snames = []
   nstr = String( copy(str[1:slen]) )
   if ln == stride - 1
      push!(snames,strip(nstr[1:ln-1]))
   else
      si = 1
      sf = 0
      for i = 1:ln
         sf += 1
         if sf == stride
            elm = nstr[si:i-1]
            elm = strip(elm)
            push!(snames,elm)
            si = i + 1
            sf = 0
         end
      end
   end
   return snames
end

function get_specie_name(str::Array{UInt8}, ln, n)
   slen = Int64(ln / n) - 1
   snames = []
   for i = 1:n
      si = (i - 1) * slen + (i - 1) + i
      sf = i * slen - 1 + (i - 1)
      elm = String( copy(str[si:sf]) )
      elm = strip(elm)
      push!(snames,elm)
   end
   return snames
end

function set_model(natoms, nspecies, 
                   coords::Array{Float64},
                   latvecs::Array{Float64},
                   origin::Array{Float64},
                   species::Array{Int32},
                   nshells::Array{Int32},
                   norbe::Array{Int32},
                   norba::Array{Int32},
                   maxns, maxno, maxs, 
                   maxo, toto,
                   angshell::Array{Int32},
                   ishell::Array{Int32},
                   posshell::Array{Int32},
                   lstr, specienames::Array{UInt8},
                   cutoff, 
                   lfn, fnames::Array{UInt8}, 
                   prnt,
                   stat_jl::Array{Int32})
    if(prnt == 1)
        @info "AITB : Julia set_model function."
    end
    try
        global acetb_dct["natoms"] = natoms
        global acetb_dct["ntypes"] = nspecies
        pos = reshape(coords,3,:) ./ bohr2ang
        cell = reshape(latvecs,3,:) ./ bohr2ang
        global acetb_dct["latvecs"] = cell[:,:]
        global acetb_dct["origin"] = origin[:]
        global acetb_dct["norbe"] = norbe
   
        elm_names = get_specie_name(specienames, lstr, nspecies)
        global acetb_dct["elm_names"] = elm_names[:]
        atnums = [] 
        atmass = [] 
        for i=1:natoms
           sym = Symbol(elm_names[species[i]])
           am = atomic_mass(sym)
           az = atomic_number(sym)
           push!(atnums,az)
           push!(atmass,am)
        end
        
        at = Atoms(; X = pos[:,1:natoms], Z = atnums, cell = cell,
                     pbc = [true, true, true])
        global acetb_dct["julip_atoms"] = at
    
        SKH_list = set_SK_orbitals(nspecies,nshells,angshell)
        global acetb_dct["SKH_list"] = SKH_list
       
        model_files = get_list_str(fnames, lfn, 200)
        global acetb_dct["model_files"] = model_files

        filedata = read_json(model_files[1])
        if(prnt == 1)
            println("read done.")
        end
        global acetb_dct["file_data"] = filedata
        trdata = filedata["training_datasets"]
        cutoff_params = filedata["model"]["cutoff_params"]
        fit_params = filedata["model"]["fit_params"]
        predict_params = filedata["model"]["predict"]
        if predict_params == 0
           Bpredict = false
        else
           Bpredict = true
        end
        if Bpredict
            if(prnt == 1)
                println("Fitting...")
            end
            Bint_table, cutf_func = predict(trdata, cutoff_params, fit_params)
            global Bondint_table = Bint_table
            global cutoff_func = cutf_func
            if(prnt == 1)
                println("Done fitting.")
            end
        else
            if(prnt == 1)
                println("Setting Bonds table...")
            end
            HSfile = filedata["HS_datasets"][1]
            if(prnt == 1)
                println("Reading file:", HSfile)
            end
            HSdata = h5read_SK(HSfile; get_HS=true, get_atoms=true, get_metadata=true, get_energies=true)
            global saveh5_HH = permutedims(HSdata[2][1], [3, 2, 1])
            global saveh5_SS = permutedims(HSdata[2][2], [3, 2, 1])
            global saveh5_satoms = HSdata[6]
            if(prnt == 1)
               println("Done.")
            end
        end
        global onsite_vals = filedata["onsite-terms"]
        stat_jl[1] = 0
    catch
        stat_jl[1] = 1
        for (exc, bt) in Base.catch_stack()
            if(prnt == 1)
                showerror(stdout, exc, bt)
                println()
            end
        end
    end
    if(prnt == 1)
        @info "AITB : Done at Julia module."
    end
end

function model_predict(iatf, iatl, natoms, 
                       coords::Array{Float64},
                       latvecs::Array{Float64},
                       nspecies, nspecarr, species::Array{Int32},
                       nH, H::Array{Float64},
                       nS, S::Array{Float64},
                       nneigh::Array{Int32},
                       ineigh::Array{Int32},
                       ipair::Array{Int32}, cutoff,
                       prnt,
                       stat_jl::Array{Int32})
    if(prnt == 1)
        @info "AITB : Julia set_predict function."
    end
    try
        pos = reshape(coords,3,:) ./ bohr2ang
        cell = reshape(latvecs,3,:) ./ bohr2ang
        
        elm_names = acetb_dct["elm_names"]
        atnums = [] 
        atmass = [] 
        for i=1:natoms
           sym = Symbol(elm_names[species[i]])
           am = atomic_mass(sym)
           az = atomic_number(sym)
           push!(atnums,az)
           push!(atmass,am)
        end
        
        at = Atoms(; X = pos[:,1:natoms], Z = atnums, cell = cell,
                     pbc = [true, true, true])

        SKH_list = acetb_dct["SKH_list"]
       
        model_files = acetb_dct["model_files"]
        filedata = acetb_dct["file_data"]
        predict_params = filedata["model"]["predict"]
        if predict_params == 0
           Bpredict = false
        else
           Bpredict = true
        end
        onsite_vals = filedata["onsite-terms"]
        onsite_terms = [onsite_vals[elm_names[species[a]]] for a=1:natoms]
        norbe = acetb_dct["norbe"]
        if(prnt == 1)
            println("AITB : Calculating bond integrals...")
        end
        buildHS(SKH_list, H, S, iatf, iatl, pos, species, nneigh, ineigh, ipair, norbe, onsite_terms, at, Bondint_table, cutoff_func, cutoff; Bpredict=Bpredict, HH=saveh5_HH, SS=saveh5_SS, supercell_atoms=saveh5_satoms, prt=prnt)
        stat_jl[1] = 0
    catch
        stat_jl[1] = 1
        for (exc, bt) in Base.catch_stack()
            if(prnt == 1)
                showerror(stdout, exc, bt)
                println()
            end
        end
    end
    if(prnt == 1)
        @info "Done at Julia module."
    end
end

function acetb_greetings()
   ctx = Pkg.Operations.Context()
   acetb_version = string(ctx.env.manifest[ctx.env.project.deps["ACEtb"]].version)
   acetb_hash = split(string(ctx.env.project.deps["ACEtb"]),"-")[1]
   println("  \033[31m\033[1m  ___\033[38;5;208m ___\033[32m ___\033[34m _____\033[38;5;54m\033[1m _    \033[0m       \n"*
           "  \033[31m\033[1m |_  \033[38;5;208m|  _\033[32m| __\033[34m|_   _\033[38;5;54m\033[1m| |_   \033[0m\033[38;5;141m\033[1mVersion\033[0m   \n"*
           "  \033[31m\033[1m | . \033[38;5;208m| |_\033[32m| __|\033[34m | | \033[38;5;54m\033[1m| . |  \033[0m\033[38;5;200m\033[1mv",acetb_version,"\033[0m   \n"*
           "  \033[31m\033[1m |___\033[38;5;208m|___\033[32m|___|\033[34m |__|\033[38;5;54m\033[1m|___|  \033[0m\033[1m[",acetb_hash,"]\033[0m\n")
   println("   ┌───────────────────────────┐\n"*
           "   │   Developers:             │\n"*
           "   ├───────────────────────────┤\n"*
           "   │     Geneviève Dusson      │\n"*
           "   │     Berk Onat             │\n"*
           "   │     Reinhard Maurer       │\n"*
           "   │     Christoph Ortner      │\n"*
           "   │     James R. Kermode      │\n"*
           "   └───────────────────────────┘\n")
   flush(stdout)
end

end # End of module
