
function simulate_args!(Δt, p)
    # determine which queries the participant responded to
    responded = convert.(Bool, rand(Binomial(1, p), length(Δt)))

    # if the participant did not respond to a query, add its time to the next
    for i in 1:findlast(responded)
        if !responded[i]
            Δt[i + 1] += Δt[i]
        end
    end

    return Δt[responded]
end

simulate_args(Δt, p) = simulate_args!(copy(Δt), p)

function simulate_args!(Δt, n, p)
    # group indices
    group_indices = ones(Int, n)

    # vector of vectors of elapsed times
    vectors = map(i -> simulate_args(Δt, p), 1:n)

    # subject indices
    subject_indices = vcat(map(i -> repeat([i], length(vectors[i])), 1:n)...)

    # hours since previous queries
    Δt = vcat(vectors...)

    # jump sizes
    M = rand([-1, 1, 1], length(Δt)) .* rand(range(0, 0.01; length = 10), length(Δt))
    M⁺ = map(x -> x > 0 ? x : 0, M)
    M⁻ = map(x -> x < 0 ? x : 0, M)

    # hours between jumps and queries
    Δτ = map(x -> rand(Uniform(0, x)), Δt)

    first_indices = findall(diff([0; subject_indices]) .!= 0)
    is_transition = trues(length(Δt))
    is_transition[first_indices] .= false

    return (; group_indices, subject_indices, Δt, M⁺, M⁻, Δτ, is_transition)
end

simulate_args(Δt, n, p) = simulate_args!(copy(Δt), n, p)