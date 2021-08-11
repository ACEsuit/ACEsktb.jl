
module AnalysisTools

using LinearAlgebra
using Plots
using Plots.PlotMeasures
using DataFrames
using Glob
using JSON

using ACEtb.FitTools

export plot_bands, plot_bands!, tabulate_results, comp_plot

function plot_bands(fc::FitContext, phase::Symbol, fermi_level::Float64, bands::Matrix; plotobj=nothing, label=nothing, color=:black)
    @assert phase in [:FCC, :BCC]
    phase = String(phase)
    
    bandpath = JSON.parsefile(joinpath(fc.base_dir, "$(phase)_tb_bandpath.json"))
    kpoints = bandpath["kpoints"]
    special_points = bandpath["special_points"]

    k0 = 0.0
    x_k = []
    tick_p = [k0]
    tick_k = phase == "FCC" ? ["W"] : tick_k = ["G"]
    push!(x_k, k0)
    for xi=2:size(kpoints, 1)
        k0 += norm(kpoints[xi] .- kpoints[xi-1])
        push!(x_k, k0)
        indict = [k for (k, v) in special_points if v==kpoints[xi] ]
        if length(indict) > 0
            push!(tick_p, k0)
            push!(tick_k, indict[1])
        end
    end

    plotobj === nothing && (plotobj = plot())
    for i = 1:size(bands, 2)
        thislabel = i == 1 ? label : false
        plot!(plotobj, x_k, bands[:, i], color=color, lw=1, label=thislabel)
    end
    for p=1:length(tick_p)
        plot!(plotobj, [tick_p[p]], lw=0.5, c=:gray, seriestype=:vline, label=false)
    end
    plot!(plotobj, [fermi_level], lw=0.5, c=:black, ls=:dash, seriestype=:hline, label=false)
    plot!(plotobj; fg_legend=:transparent, bg_legend=:transparent, 
          grid=false, framestyle = :box, legend=:bottomleft)
    xticks!(plotobj, ([t for t in tick_p],[k for k in tick_k]), fontsize=14)
    ylabel!(plotobj, "Energy (eV)")
    
    xlims!(plotobj, 0.0, maximum(x_k))
    ylims!(plotobj, -20, fermi_level + 10)
    
    plot(plotobj, dpi=120)
    return plotobj
end

function plot_bands(fc::FitContext, phase::Symbol, param_hash::String; label=nothing, kwargs...)
    fermi_level, bands = read_bands(fc, phase, param_hash)
    label === nothing && (label = param_hash)
    plot_bands(fc, phase, fermi_level, bands; kwargs...,
               label=label)
end

plot_bands!(plotobj::Plots.Plot, args...; kwargs...) = plot_bands(args...; plotobj=plotobj, kwargs...)


function tabulate_results(fc::FitContext)
    df = DataFrame(param_hash=String[],
                   order=Int[], 
                   degree=Int[],
                   env_deg=Int[],
                   rcut=Float64[],
                   renv=Float64[], 
                   pcut=Int[], 
                   N_train=Int[],
                   N_bond=Int[],
                   N_basis=Int[], 
                   error_train=Float64[],
                   error_test=Union{Missing,Float64}[],
                   band_error_FCC=Union{Missing,Float64}[], 
                   band_error_BCC=Union{Missing,Float64}[], 
                   pot_file=String[], 
                   training_datasets=String[],
                   method=String[], 
                   rtol=Union{Missing,Float64}[])
    for pot_file in glob(joinpath(splitpath(fc.output_dir)[end], "pot_*/pot.json"))
        param_hash = match(r"pot_(.*)/pot.json", pot_file).captures[1]
        pot_dict = JSON.parsefile(pot_file)
        band_error = get(pot_dict, "band_error", missing)
        band_error_FCC = band_error_BCC = missing
        if band_error !== missing
            band_error_FCC = band_error["FCC"]
            band_error_BCC = band_error["BCC"]
        end
        fit_params = pot_dict["model"]["fit_params"]
        cutoff_params = pot_dict["model"]["cutoff_params"]
        row = (param_hash,
               fit_params["order"],
               fit_params["degree"],
               fit_params["env_deg"],
               cutoff_params["rcut"],
               cutoff_params["renv"],
               cutoff_params["pcut"],
               length(pot_dict["training_datasets"]),
               pot_dict["size_A"][1],
               pot_dict["size_A"][2],
               pot_dict["error_train"],
               get(pot_dict, "test_error", missing),
               band_error_FCC,
               band_error_BCC,
               pot_file,
               join(pot_dict["training_datasets"], ","),
               get(fit_params, "method", "QR"),
               get(fit_params, "rtol", missing))
        push!(df, row)
    end

    # add a descriptive label for the types of training configs
    df.config = ifelse.(occursin.("005", df.training_datasets), 
                        ifelse.(occursin.("010", df.training_datasets), "mixT", "lowT"),
                        "highT")

    return df
end

function ref_plots(fc::FitContext)
    fcc_fermi_level_FHI, fcc_bands_FHI = read_FHI_bands(fc, :FCC)
    bcc_fermi_level_FHI, bcc_bands_FHI = read_FHI_bands(fc, :BCC)
    fcc_plot = plot_bands(fc, :FCC, fcc_fermi_level_FHI, fcc_bands_FHI; label="FHIaims ref", color=:red)
    bcc_plot = plot_bands(fc, :BCC, bcc_fermi_level_FHI, bcc_bands_FHI; label="FHIaims ref", color=:red)
    display(plot(fcc_plot, bcc_plot, title=["FCC" "BCC"]))
    return fcc_plot, bcc_plot
end

function add_plots(fc::FitContext, potential, color, fcc_plot, bcc_plot; label=nothing)
    plot_bands!(fcc_plot, fc, :FCC, potential, color=color, label=label)
    plot_bands!(bcc_plot, fc, :BCC, potential, color=color, label=label)
    p = plot(fcc_plot, bcc_plot, title=["FCC" "BCC"])
    display(p)
    return p
end

function comp_plot(fc::FitContext, potential; color=:blue, label=nothing)
    fcc_plot, bcc_plot = ref_plots(fc)
    add_plots(fc, potential, color, fcc_plot, bcc_plot, label=label)
end

end