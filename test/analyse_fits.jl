using Plots
using Plots.Measures
using DataFrames
using LinearAlgebra
using JSON

using ACEtb
using ACEtb.FitTools
using ACEtb.AnalysisTools
using LaTeXStrings

##

const DFTB_PATH = "$(@__DIR__)/../dftbplus/build/install/bin/"
const OUTPUT_DIR = joinpath(@__DIR__, "Results_fits/")
const PARAM_FILE = "Al-Al.acetb.json"

fc = FitContext(@__DIR__, DFTB_PATH, OUTPUT_DIR, PARAM_FILE)

results = tabulate_results(fc)
filter!(x -> !startswith(x.param_hash, "order"), results)

subresults = dropmissing(results)
subresults = subresults[(subresults.rcut .== 10.0) .& (subresults.N_train .<= 2) .& (subresults.rtol .== 1e-8), :]

gdf = groupby(subresults, :config)

degplots = [plot(; yscale=:log10, legend=:bottomleft, xlabel="Degree",
                ylabel="Training error")
            plot(; yscale=:log10, legend=:false, xlabel="Degree",
                ylabel="Band error FCC", ylim=(10^0.2, 10^0.4)) 
            plot(; yscale=:log10, legend=:false, xlabel="Degree",
                ylabel="Band error BCC")] #ylim=(10^0.2, 10^0.4))]
for df in gdf
    for order in (2, 3)
        nrow(df) <= 1 && continue
        degree = df[df.order .== order, :degree]
        train_err = df[df.order .== order, :error_train]
        band_err_fcc = df[df.order .== order, :band_error_FCC]
        band_err_bcc = df[df.order .== order, :band_error_BCC]
        config = df.config[df.order .== order][1]
        sortorder = sortperm(degree)
        plot!(degplots[1], degree[sortorder], 
              train_err[sortorder], linewidth=2,
              label="Order $order $config")
        plot!(degplots[2], degree[sortorder], 
              band_err_fcc[sortorder], linewidth=2,
              label="Order $order $config")
        plot!(degplots[3], degree[sortorder], 
              band_err_bcc[sortorder], linewidth=2,
              label="Order $order $config")
    end
end


rcut_plot = [plot(; xlabel="Cutoff radius", ylabel="Training Error", legend=:bottomleft, yscale=:log10),
             plot(; xlabel="Cutoff radius", ylabel="FCC Band Error", legend=false, yscale=:log10),
             plot(; xlabel="Cutoff radius", ylabel="BCC Band Error", legend=false, yscale=:log10)]

order = 2
for (color, config_type) in enumerate(["highT", "lowT", "mixT"])
    rcut_mask =  ((results.config .== config_type) .& 
                    (results.rcut .<= 11.0) .&
                    (results.order .== order) .& 
                    (results.rcut .> 6.0) .&
                    (results.N_train .<= 2) .&
                    (results.pcut .== 1))
    if config_type == "highT"
        rcut_mask .&= results.degree .== 10
    else
        rcut_mask .&= results.degree .== 10
    end
    rcut_mask = findall(rcut_mask)
    rcut_mask = rcut_mask[sortperm(results[rcut_mask, :rcut])]
    @show order config_type rcut_mask, results[rcut_mask, [:rcut, :pot_file]]

    plot!(rcut_plot[1], results[rcut_mask, :rcut], marker=:dot,
            markersize=3, markerstrokecolor=color, color=color,
            results[rcut_mask, :error_train], linewidth=2,
            label="Order $order $config_type")

    plot!(rcut_plot[2], results[rcut_mask, :rcut], marker=:dot,
            markersize=3, markerstrokecolor=color, color=color,
            results[rcut_mask, :band_error_FCC], linewidth=2,
            label="Order $order $config_type")

    plot!(rcut_plot[3], results[rcut_mask, :rcut], marker=:dot,
            markersize=3, markerstrokecolor=color, color=color,
            results[rcut_mask, :band_error_BCC], linewidth=2,
            label="Order $order $config_type")

end

rp = plot(rcut_plot..., layout=(1, 3), size=(1200, 400))
savefig("rcut.pdf")

dp = plot(degplots..., layout=(1, 3), size=(1200, 400))
savefig("degree.pdf")

# bp = comp_plot(fc, "order2_deg10_envdeg4_rcut11_renv3.0_conf15_pcut1_Ntrain1",
#                label="mixT Order 2 Deg 10 rcut 11")
# savefig("bands.pdf")

# bp


##

N_train = 100 * 321

p1 = scatter(results.N_bond, results.error_train ./ sqrt.(results.N_bond), xscale=:log10, yscale=:log10, label="Train", legend=:bottomright)
xlabel!(L"N_\mathrm{bonds}")
ylabel!("RMSE Hamiltonian")

scatter!(results.N_bond, results.error_test ./ sqrt(N_train), xscale=:log10, yscale=:log10, label="Test")
# xlabel!(L"N_\mathrm{bonds}")
# ylabel!("RMSE Hamiltonian")

p3 = scatter(results.N_bond, results.band_error_FCC, xscale=:log10, yscale=:log10, label="FCC")
scatter!(results.N_bond, results.band_error_BCC, xscale=:log10, yscale=:log10, label="BCC")

xlabel!(L"N_\mathrm{bonds}")
ylabel!(L"\mathrm{Band\; error}")

p = plot(p1, p3, size=(600, 300), layout=(1, 2), leftmargin=5Plots.mm, bottommargin=5Plots.mm)

##

println("Loading test data...")
param_hash = results[argmin(results.N_train), :param_hash] # pick smallest model
pot_dict, _, Bint_table, cutf_func = load_model(fc, param_hash)
test_dataset = pot_dict["test_datasets"]
test_data = ACEtb.Utils.get_data(test_dataset, cutf_func, ACEtb.Bonds.get_env)
println("done")

##

test_errors = Dict()
train_errors = Dict()

for param_hash in sort(results, :N_train)[1:5, :param_hash]
    @show param_hash
    pot_dict, _, Bint_table, cutf_func = load_model(fc, param_hash)
    test_errors[param_hash] = hcat(test_error(test_data, Bint_table)...)
    println("Test error: 2-norm $(norm(test_errors[param_hash])), inf-norm $(norm(test_errors[param_hash],Inf)), rms $(mean(sqrt(test_errors[param_hash].^2)))")

    train_dataset = pot_dict["training_datasets"]
    train_data = ACEtb.Utils.get_data(test_dataset, cutf_func, ACEtb.Bonds.get_env)    
    train_errors[param_hash] = hcat(test_error(train_data, Bint_table)...)
    println("Test error: 2-norm $(norm(train_errors[param_hash])), inf-norm $(norm(train_errors[param_hash],Inf)), rms $(mean(sqrt(train_errors[param_hash].^2)))")
end
