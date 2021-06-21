using Documenter, ACEtb

makedocs(sitename="ACEtb.jl Documentation",
         pages = [
        "Home" => "index.md",
        "Introduction" => "intro.md",
        "Getting Started" => "gettingstarted.md",
        "Installation" => "installation.md",
        "ACEtb Model" => [
            "What is ACEtb Model" => "acetb_intro.md",
            "How to fit ACEtb Model" => "fitting.md",
            "Basic Usage" => "acetb_usage.md",
        ],
        "ACEtb API: `lib_acetb`" => [
            "DFTB+ Integration" => "dftbplus_acetb.md",
        ],
        "Developer Docs" => "devel.md"
        ])

# deploydocs(
#     repo = "github.com/ACEsuit/ACEtb.jl.git",
# )

