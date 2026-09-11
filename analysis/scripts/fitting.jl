include("../startup.jl")

####################################################################################################
# REAL DATA
####################################################################################################

folder = joinpath("models", "real")

df_data = d"Dataset"

variable = :Valence

args = prepare_args(df_data, variable)

Threads.@threads for model in MODELS
    filename = make_filename(folder, model, variable)

    fit_model(model, filename, args)
end

variables = [:PositiveEvents, :NegativeEvents]
events = [args.M⁺ .> 0, args.M⁻ .> 0]

Threads.@threads for i in eachindex(variables)
    filename = make_filename(folder, EventProbabilityModel, variables[i])

    fit_model(EventProbabilityModel, filename, merge(args, (; e = events[i])))
end

####################################################################################################
# SIMULATED DATA
####################################################################################################

include("simulation.jl")

for j in eachindex(MODELS)
    folder = joinpath("models", "simulation", string(MODELS[j]))

    Threads.@threads for i in eachindex(args)
        # fit all models using the simulated data
        for model in MODELS
            filename = joinpath(folder, string(model) * " ($i).jld2")

            fit_model(model, filename, merge(args[i], (; z_obs = params[i, j].z_obs)); iterations = 1000)
        end
    end
end