using Random
using JSON

ENV["PYTHONPATH"] = "$(@__DIR__)/../dftbplus/build/install/lib/python3.7/site-packages/"
ENV["JULIA_PROJECT"] = "$(@__DIR__)"

const DFTB_PATH = "$(@__DIR__)/../dftbplus/build/install/bin/"
const OUTPUT_DIR = joinpath(@__DIR__, "Results_fits/")
const DATABCC = "/home/eng/essswb/ownCloud/SK_AITB/BCC/without_HS/MD-20K"  #"data/BCC_MD_20K/MD-20K"
const DATAFCC = "/home/eng/essswb/ownCloud/SK_AITB/FCC/without_HS/MD-500K/"  #"data/FCC_MD_500K/MD-500K"
const PARAM_FILE = "Al-Al.acetb.json"

# always include lattice expansion data?
# FCC: /home/eng/essswb/ownCloud/SK_AITB/FCC/without_HS/LATTICE_EXTENSION/SUPERCELL-9x9x9/SK-supercell-noHS-000.h5
# BCC: /home/eng/essswb/ownCloud/SK_AITB/BCC/without_HS/LATTICE_EXPANSION/SK-supercell-noHS-000.h5

using ACEtb.FitTools
fc = FitContext(@__DIR__, DFTB_PATH, OUTPUT_DIR, PARAM_FILE)

params = Dict(
    :order => 3,
    :deg => 12,
    :rcut => 10.0,
    :renv => 3.0,
    :envdeg => 4,
    :pcut => 1,
    :method => :RegRRQR,
    :rtol => 1e-8
)

bcc_files = readlines(pipeline(`find $DATABCC -name \*.h5`))
fcc_files = readlines(pipeline(`find $DATAFCC -name \*.h5`))
bcc_data = reshape([ (file, index) for file in bcc_files, index in 1:729 ], :)
fcc_data = reshape([ (file, index) for file in fcc_files, index in 1:729 ], :)
all_data = [bcc_data; fcc_data]

if !isfile("train_test.json")
    perm = randperm(length(all_data))
    open("train_test.json", "w") do json_file
        JSON.print(json_file, perm)
    end
else
    perm = JSON.parsefile("train_test.json")
end
train = perm[1:length(perm)÷2]
test = perm[length(perm)÷2+1:end]

training_data = all_data[train]
test_data = all_data[test]

cparams = copy(params)

N_train = parse(Int, ARGS[1])
N_test = 100

@show N_train, N_test
model, param_hash, Bint_table, cutf_func = fit_model(fc, cparams, 
                                                     training_data[1:N_train], 
                                                     test_data[1:N_test],
                                                     "exact-onsite-BCC.json")

