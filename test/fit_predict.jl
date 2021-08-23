using Random
using JSON

ENV["PYTHONPATH"] = "$(homedir())/lib_acetb/dftbplus/build/install/lib/python3.7/site-packages/"
ENV["JULIA_PROJECT"] = "$(@__DIR__)"

const DFTB_PATH = "$(homedir())/lib_acetb/dftbplus/build/install/bin/"
const OUTPUT_DIR = joinpath(@__DIR__, "fits/")
const DATABCC = "data/BCC/"
const DATAFCC = "data/FCC/"
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

bcc_files = readlines(pipeline(`find -L $DATABCC -name SK-supercell-\*.h5`))
fcc_files = readlines(pipeline(`find -L $DATAFCC -name SK-supercell-\*.h5`))
bcc_data = reshape([ (file, index) for file in bcc_files, index in 1:729 ], :)
fcc_data = reshape([ (file, index) for file in fcc_files, index in 1:729 ], :)
all_data = [bcc_data; fcc_data]

@info "Found $(length(bcc_data)) FCC and $(length(fcc_data)) FCC atomic environments"

if !isfile("train_test.json")
    @info "Generating train/test split..."
    perm = randperm(length(all_data))
    open("train_test.json", "w") do json_file
        JSON.print(json_file, perm)
    end
else
    @info "Loading train/test split"
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
                                                     "data/exact-onsite-BCC.json")

