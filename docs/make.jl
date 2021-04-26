using Documenter, ACEtb

makedocs(sitename="ACEtb.jl Documentation",
         pages = [
        "Home" => "index.md",
        "Introduction" => "intro.md",
        "Getting Started" => "gettingstarted.md",
        "ACEtb Model" => [
            "What is ACEtb Model" => "acetb_intro.md",
            "Basic Usage" => "acetb_usage.md",
            "How to fir ACEtb Model" => "fitting.md",
        ],
        "ACEtb APIs" => [
            "DFTB+ Integration" => "dftbplus_acetb.md",
        ],
        "Developer Docs" => "devel.md"
        ])

# deploydocs(
#     repo = "github.com/ACEsuit/ACEtb.jl.git",
# )

