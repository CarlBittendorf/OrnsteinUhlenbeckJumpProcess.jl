include("../startup.jl")
include("simulation.jl")
include("pointwise.jl")

psisloos = Array{DataFrame, 3}(undef, length(args), length(MODELS), length(MODELS))

for j in eachindex(MODELS)
    folder = joinpath("models", "simulation", string(MODELS[j]))

    Threads.@threads for i in eachindex(args)
        for (k, model) in enumerate(MODELS)
            filename = joinpath(folder, string(model) * " ($i).jld2")
            chains = load(filename, "chains")

            psisloos[
                i, j, k] = @chain begin
                model(; merge(args[i], (; z_obs = params[i, j].z_obs))...)
                psis_loo(chains; source = "mcmc")
                getproperty(:estimates)

                DataFrame
                rename([:Statistic, :Column, :Value])

                # add columns for simulated model, fitted model and iteration
                transform(
                    All() => ByRow((x...) -> string(MODELS[j])) => :ModelSimulation,
                    All() => ByRow((x...) -> string(model)) => :ModelFit,
                    All() => ByRow((x...) -> i) => :Iteration
                )
            end
        end
    end
end

df_psisloo = @chain psisloos begin
    vcat(_...)
    select(:ModelSimulation, :ModelFit, :Iteration, :Statistic, :Column, :Value)
end

CSV.write(joinpath("data", "Model Recovery PSIS-LOO.csv"), df_psisloo)