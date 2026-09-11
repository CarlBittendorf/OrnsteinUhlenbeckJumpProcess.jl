
####################################################################################################
# EXPORTS
####################################################################################################

export hours_elapsed

####################################################################################################
# IMPLEMENTATIONS
####################################################################################################

"""
    hours_elapsed(x)

Compute the time difference, in hours, between each element of a datetime-like
vector and its previous element.

This function subtracts the lagged version of `x` (with the first element used
as a default for the initial lag) and converts the resulting time deltas to
hours.

# Arguments
- `x::AbstractVector`: A vector of `DateTime`, `Date`, or other
  `Dates.Period`-compatible values.

# Returns
- `Vector{Float64}`: A vector where each entry represents the number of hours
  elapsed since the previous timestamp.

# Details
- Internally uses `ShiftedArrays.lag(x; default = first(x))` to align each
  element with its predecessor.
- `Dates.value` converts the resulting `Millisecond` differences into integers.
- Dividing by `3_600_000` converts milliseconds to hours.

# Example
```julia
julia> using Dates

julia> x = [
  DateTime("2026-01-01T12:00:00"),
  DateTime("2026-01-01T15:00:00"),
  DateTime("2026-01-02T12:00:00")
  ]

julia> hours_elapsed(x)
3-element Vector{Float64}:
  0.0
  3.0
 21.0
```
"""
function hours_elapsed(x)
    Dates.value.(x .- ShiftedArrays.lag(x; default = first(x))) ./ 3_600_000
end