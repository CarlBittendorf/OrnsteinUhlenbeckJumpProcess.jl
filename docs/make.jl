using OrnsteinUhlenbeckJumpProcess
using Documenter
using DocumenterVitepress

makedocs(;
    modules = [OrnsteinUhlenbeckJumpProcess],
    authors = "Carl Bittendorf <carl.bittendorf@gmail.com> and contributors",
    sitename = "OrnsteinUhlenbeckJumpProcess.jl",
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl",
        devbranch = "main",
        devurl = "dev"
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Models" => [
            "Overview" => "models/index.md",
            "OrnsteinUhlenbeckModel" => "models/ornstein_uhlenbeck.md",
            "OrnsteinUhlenbeckJumpModel" => "models/ou_jump.md",
            "EventProbabilityModel" => "models/event_probability.md",
        ],
        "Examples" => [
            "Analysis workflow" => "examples/analysis.md",
            "Jump model workflow" => "examples/jump_model.md",
        ],
        "Resources" => [
            "Extending the package" => "resources/extending.md",
            "API reference" => "resources/api.md",
        ],
    ]
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/CarlBittendorf/OrnsteinUhlenbeckJumpProcess.jl.git",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true
)

