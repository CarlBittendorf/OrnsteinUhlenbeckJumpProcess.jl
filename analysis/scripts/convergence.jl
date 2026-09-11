include("../startup.jl")

####################################################################################################
# REAL DATA
####################################################################################################

folder = joinpath("models", "real")

rhat_values = []
ess_values = []

for model in MODELS
    chains = @chain folder begin
        make_filename(model, :Valence)
        load("chains")
    end

    append!(rhat_values, rhat(chains).nt.rhat)
    append!(ess_values, ess(chains).nt.ess)
end

maximum(rhat_values)
median(rhat_values)
median(ess_values)
minimum(ess_values)

rhat_values = []
ess_values = []

for variable in [:PositiveEvents, :NegativeEvents]
    chains = @chain begin
        make_filename(folder, EventProbabilityModel, variable)
        load("chains")
    end

    append!(rhat_values, rhat(chains).nt.rhat)
    append!(ess_values, ess(chains).nt.ess)
end

maximum(rhat_values)
median(rhat_values)
median(ess_values)
minimum(ess_values)

####################################################################################################
# SIMULATED DATA
####################################################################################################

folder = joinpath("models", "simulation")

rhat_values = []
ess_values = []

for j in eachindex(MODELS)
    folder = joinpath("models", "simulation", string(MODELS[j]))

    for i in 1:100
        for model in MODELS
            filename = joinpath(folder, string(model) * " ($i).jld2")

            chains = load(filename, "chains")

            append!(rhat_values, rhat(chains).nt.rhat)
            append!(ess_values, ess(chains).nt.ess)
        end
    end
end

maximum(rhat_values)
median(rhat_values)
median(ess_values)
minimum(ess_values)