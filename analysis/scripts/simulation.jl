
# ensure reproducibility
Random.seed!(1)

# time in hours between successive queries
Δt = begin
    Δt = ones(48)
    Δt[13] = Δt[25] = Δt[37] += 12

    return Δt
end

# create 100 different arguments
args = [simulate_args(Δt, 25, 0.9) for _ in 1:100]

# create 100 different combinations of parameter values for each model
params = [rand(model(; merge(x, (; z_obs = missing))...)) for x in args, model in MODELS]