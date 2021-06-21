using ACEtb
using Test

@testset "ACEtb.jl" begin
    # Write your tests here.

    include("test_bond1p.jl")
    
    # SlaterKoster tests
    include("test_basics.jl")
    include("test_sk.jl")
    #include("test_kwon.jl")
end
