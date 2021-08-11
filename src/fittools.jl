module FitTools

using LinearAlgebra
using JSON
using DelimitedFiles
using Glob
using Printf

using JuLIP, ACE, ACEtb
using ACEtb.Bonds: get_env
using ACEtb.Utils: read_json, write_json, get_data
using ACEtb.Predictions: predict, train_and_predict

struct FitContext
    base_dir::String
    dftb_path::String
    output_dir::String
    param_basename::String
end

function make_params(mode::Symbol, training_datasets, test_datasets; rcut=12, 
    renv=3, zenv=1, deg=5, order=2, 
    envdeg=4, pcut=1, unitcell=false, 
    method=:QR, rtol=1e-15)

    if mode == :fit
        predict = 1
        fit = 1
    elseif mode == :predict
        predict = 1
        fit = 0
    else
        error("mode should be either :fit or :predict")
    end

    param_dict = Dict(
        "training_datasets" => training_datasets,
        "test_datasets" => test_datasets,
        "model" => Dict(
                "predict" => predict,
                "fit" => fit,
                "cutoff_params" => 
                        Dict("rcut" => rcut,
                            "renv" => renv,
                            "zenv" => zenv,
                            "pcut" => pcut),
                "fit_params" => Dict(
                    "degree" => deg,
                    "order" => order,
                    "env_deg" => envdeg,
                    "method" => method,
                    "rtol" => rtol)),
        "onsite-terms" => Dict(
        "Al" => [
            -5.5708730114137865e+01,
            -4.0432979354580754e+00,
            -4.9611320942331338e-01,
            -2.6243905672312331e+00,
            -2.6243905672300967e+00,
            -2.6241786637179012e+00,
            -4.6415677514380915e-01,
            -4.6415677508608638e-01,
            -4.6415677490900636e-01,
            -2.8478483849360947e-01,
            -2.8478483839241397e-01,
            -2.8478483828616641e-01,
            -1.8885840374342641e-01,
            -1.8885840355491695e-01
        ]))
    return param_dict, params_to_string(param_dict)
end

function params_to_string(params::Dict)
    input_dict = Dict{String}{Any}()
    for (k, v) in params
        if k in ["training_datasets", "test_datasets", "model", "onsite-terms"]
            input_dict[k] = copy(v)
        end
    end
    delete!(input_dict["model"], "fit")
    delete!(input_dict["model"], "predict")
    return @sprintf "%x" hash(input_dict)
end

function do_fit(pot_dict::Dict, elm_names=["Al"])
    trdata = pot_dict["training_datasets"]
    cutoff_params = pot_dict["model"]["cutoff_params"]
    fit_params = pot_dict["model"]["fit_params"]
    @info "│    Fitting..."
    Bint_table, cutf_func, train_dict = train_and_predict(trdata, elm_names, cutoff_params, 
                                                          fit_params; MPIproc=1)

    return train_dict, Bint_table, cutf_func
end

function load_model(fc::FitContext, param_hash::String)
    pot_file = joinpath(fc.output_dir, "pot_" * param_hash, "pot.json")
    pot_dict = read_json(pot_file)

    cutoff_params = pot_dict["model"]["cutoff_params"]
    fit_params = pot_dict["model"]["fit_params"]
    @info "│    Loading model..."
    Bint_table, cutf_func = predict(pot_dict, cutoff_params, fit_params; MPIproc=1)
    @info "│    Model is loaded."
    return pot_dict, param_hash, Bint_table, cutf_func
end

function test_error(test, Bint_table)
    errors = []
    for (i, (R0, Renv, V)) in enumerate(test)
        push!(errors, Bint_table(R0, Renv) - V)
    end
    return errors
end

function fit_model(fc::FitContext, params::Dict, training_datasets, test_datasets; band_errors=true)
    pot_dict, param_hash = make_params(:fit, training_datasets, test_datasets; params...)
    pot_dir = joinpath(fc.output_dir, "pot_" * param_hash)
    pot_file = joinpath(pot_dir, "pot.json")

    if isfile(pot_file)
        println("Skipping fit as $pot_file already present")
        pot_dict, param_hash, Bint_table, cutf_func = load_model(fc, param_hash)
    else
        isdir(pot_dir) || mkdir(pot_dir)
        println("Fitting model $param_hash in directory $pot_dir")
    
        train_dict, Bint_table, cutf_func = do_fit(pot_dict, ["Al"])
        merge!(pot_dict, train_dict)
        @info "Train error $(train_dict["error_train"])"    
    end

    pot_dict["model"]["fit"] = 0 # disable fitting for subsequent runs

    if "test_error" ∉ keys(pot_dict) || "test_error_all" ∉ keys(pot_dict)
        # compute test error on independent set
        test_dataset = pot_dict["test_datasets"]
        test_data = ACEtb.Utils.get_data(test_dataset, cutf_func, get_env)
        pot_dict["test_error_all"] = test_error(test_data, Bint_table)
        pot_dict["test_error"] = norm(test_error_all)
        @info "Test error $(pot_dict["test_error"])"
    end

    # save so we can compute Fermi level; will be overwritten later
    open(joinpath(pot_dir, "pot.json"), "w") do json_file
        JSON.print(json_file, pot_dict, 4)
    end

    # compute and save Fermi level of the model for FCC and BCC
    if "fermi_level" ∉ keys(pot_dict)
        fcc_fermi_level = fermi_level(fc, :FCC, param_hash)
        bcc_fermi_level = fermi_level(fc, :BCC, param_hash)
        pot_dict["fermi_level"] = Dict("FCC" => fcc_fermi_level,
                                       "BCC" => bcc_fermi_level)
    end

    # include calculated Fermi level in potential file
    open(joinpath(pot_dir, "pot.json"), "w") do json_file
        JSON.print(json_file, pot_dict, 4)
    end

    # use the new model to compute FCC and BCC band errors
    band_errors && compute_band_errors(fc, param_hash)

    return pot_dict, param_hash, Bint_table, cutf_func
end

function fermi_level(fc::FitContext, phase::Symbol, param_hash::String)
    @assert phase in [:FCC, :BCC]
    println("Calculating Fermi level for $phase with params $param_hash")
    pot_dir = joinpath(fc.output_dir, "pot_" * param_hash)

    cp("dftb_in.hsd-$phase", joinpath(pot_dir, "dftb_in.hsd"), force=true)
    pot_file = joinpath(pot_dir,  "pot.json")
    cp(pot_file, joinpath(pot_dir, fc.param_basename); force=true)
    open(joinpath(pot_dir, "dftb_fermi_$phase.log"), "w") do logfile
        run(pipeline(setenv(`$(fc.dftb_path)/dftb+`, dir=pot_dir), stdout=logfile, stderr=logfile))
    end
    return parse(Float64, split(readchomp(`grep "Fermi level" $pot_dir/detailed.out`))[5])
end

function predict_bands(fc::FitContext, phase::Symbol, param_hash::String)
    @assert phase in [:FCC, :BCC]
    pot_dir = joinpath(fc.output_dir, "pot_" * param_hash)
    band_file = joinpath(pot_dir, "BANDS_$(phase).dat")
    if isfile(band_file)
        println("Skipping band calculation for $phase as $band_file already present")
        return read_bands(fc, phase, param_hash)
    end
    println("Calculating bands for $phase with params $param_hash")
    phase = String(phase)
    pot_file = joinpath(pot_dir, "pot.json")
    cp(pot_file, joinpath(pot_dir, fc.param_basename); force=true)
        
    # read template with correct band path definition for this phase
    lines = readlines("dftb_in.hsd-$phase-BANDS")

    # set the Fermi level to match value from model
    pot_dict = JSON.parsefile(pot_file)
    fermi_level = pot_dict["fermi_level"][phase]
    new_lines = []
    for line in lines
        m = match(r"(^\s*)FixedFermiLevel.*$", line)
        if m !== nothing
            line = m.captures[1] * "FixedFermiLevel = $(fermi_level)"
        end
        push!(new_lines, line)
    end
    open(joinpath(pot_dir, "dftb_in.hsd"), "w") do dftb_in
        write(dftb_in, join(new_lines, "\n"))
    end

    open(joinpath(pot_dir, "dftb_bands_$phase.log"), "w") do logfile
        run(pipeline(setenv(`$(fc.dftb_path)/dftb+`, dir=pot_dir), stdout=logfile, stderr=logfile))
        run(pipeline(setenv(`$(fc.dftb_path)/dp_bands band.out band`, dir=pot_dir), 
                     stdout=logfile, stderr=logfile))
    end
    cp(joinpath(pot_dir, "band_tot.dat"), band_file)

    bands = readdlm(band_file)[:, 2:end]
    return fermi_level, bands
end

function read_bands(fc::FitContext, phase::Symbol, param_hash::String)
    @assert phase in [:FCC, :BCC]
    phase = String(phase)

    pot_dir = joinpath(fc.output_dir, "pot_" * param_hash)
    pot_file = joinpath(pot_dir, "pot.json")
    pot_dict = JSON.parsefile(pot_file)
    fermi_level = get(pot_dict, "fermi_level", missing)
    if fermi_level !== missing
        fermi_level = get(fermi_level, phase, missing)
    end

    band_file = joinpath(pot_dir, "BANDS_$(phase).dat")
    bands = readdlm(band_file)[:, 2:end]

    return fermi_level, bands
end

function read_FHI_bands(fc::FitContext, phase::Symbol)
    @assert phase in [:FCC, :BCC]
    index = phase == :FCC ? 1 : 2
    phase = lowercase(String(phase))
    dict = JSON.parsefile(joinpath(fc.base_dir, "FCC_BCC_BANDS.json"))
    fermi_level = dict["Bands"]["fermi"][phase]
    bands = hcat(dict["Bands"][phase]["TB_FHIaims"][index]...)'
    # undo shift to put back in absolute terms
    bands = bands .+ fermi_level
    return fermi_level, bands
end

kB = 8.6173303e-5
T = 100.0

fermi_dirac(x, μ, T) = 1. / (1 + exp((x - μ) / (kB * T)))

function band_error(fermi_level, bands, fermi_level_ref, bands_ref)
    error_k =  [sum(abs.((bands[n, :] .* fermi_dirac.(bands[n, :], fermi_level_ref, T) - 
                          bands_ref[n, :] .* fermi_dirac.(bands_ref[n, :], fermi_level_ref, T)))) 
                          for n = 1:size(bands, 1)]
    RMSE = sqrt(sum(error_k .^ 2)/length(error_k))
    return RMSE
end

function compute_band_errors(fc::FitContext, param_hash::String)
    pot_file = joinpath(fc.output_dir, "pot_" * param_hash, "pot.json")
    pot_dict = JSON.parsefile(pot_file)

    fcc_fermi_level, fcc_bands = predict_bands(fc, :FCC, param_hash)
    bcc_fermi_level, bcc_bands = predict_bands(fc, :BCC, param_hash)
    
    fcc_fermi_level_FHI, fcc_bands_FHI = read_FHI_bands(fc, :FCC)
    bcc_fermi_level_FHI, bcc_bands_FHI = read_FHI_bands(fc, :BCC)    

    RMSE_FCC = band_error(fcc_fermi_level, fcc_bands, fcc_fermi_level_FHI, fcc_bands_FHI)
    RMSE_BCC = band_error(bcc_fermi_level, bcc_bands, bcc_fermi_level_FHI, bcc_bands_FHI)    

    band_error_dict = Dict("FCC" => RMSE_FCC,
                           "BCC" => RMSE_BCC)

    # save results in potential dict
    pot_dict["band_error"] = band_error_dict
    open(pot_file, "w") do json_file
        JSON.print(json_file, pot_dict, 4)
    end
    
    return band_error_dict
end

export fit_model, load_model, FitContext, read_bands, params_to_string, read_FHI_bands, test_error

end