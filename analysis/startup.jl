using OrnsteinUhlenbeckJumpProcess, Chain, DataToolkit, DataFrames, Turing, ReverseDiff,
      StableRNGs, StatsFuns, ParetoSmooth, AlgebraOfGraphics, CairoMakie, CSV, JLD2,
      Combinatorics, SummaryTables
using Dates, Statistics, Random

include("src/utils.jl")
include("src/fitting.jl")
include("src/simulation.jl")
include("src/analysis.jl")

set_aog_theme!()

const PALETTE = AlgebraOfGraphics.aog_theme().palette.color
const BLUE, ORANGE, GREEN, PURPLE, LIGHTBLUE, RED, YELLOW = PALETTE

const MODELS = [OrnsteinUhlenbeckModel, OrnsteinUhlenbeckJumpModel]

const DESCRIPTIONS = [
    "μ_group" => "Attractor",
    "θ_group" => "Half-Life of Deviations",
    "σ_group" => "Volatility",
    "α⁺_group" => "Reactivity (Positive Events)",
    "α⁻_group" => "Reactivity (Negative Events)"
]