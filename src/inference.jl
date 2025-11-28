"""
    estimate_efficiencies_from_curves(curves, λ0_d)

Given a mapping `curves` from some key (e.g. fiber symbol) to a NamedTuple
with fields `:eps` and `:λ0` (as produced by `λ0_vs_efficiency_all` or similar),
and a dictionary `λ0_d` mapping channel symbols to observed no-light
probabilities, compute the efficiency at which the simulated curve crosses the
observed λ₀ for each channel.

# Arguments
- `curves::AbstractDict{K,NamedTuple}`: maps a key (fiber or channel) to a
  NamedTuple `(eps, λ0)` where `eps` is a vector of efficiencies and `λ0` the
  corresponding no-light probabilities.
- `λ0_d::NamedTuple`: observed λ₀ per channel symbol.

# Returns
- `Dict{Symbol,Float64}` mapping channel symbol -> estimated efficiency.

The routine finds a bracket where the curve crosses the horizontal line
`λ0 = λ0_d[channel]` and linearly interpolates between neighbouring points.
If no sign change is found, it falls back to the nearest sampled point.
"""
function estimate_efficiencies_from_curves(
    curves::AbstractDict,
    λ0_d::NamedTuple
)
    eff_est = Dict{Symbol,Union{Float64,Nothing}}()

    λ0_map = Dict(pairs(λ0_d))

    for (ch, entry) in curves
        if !(ch in keys(λ0_map))
            eff_est[ch] = nothing
            continue
        end

        eps = entry.eps
        λvals = entry.λ0
        target = float(λ0_map[ch])

        # compute differences and look for sign change
        diffs = λvals .- target
        idx =
            findfirst(i -> i < length(diffs) && (diffs[i] == 0 || (diffs[i] * diffs[i + 1] < 0)), 1:(length(diffs) - 1))

        eps_cross = nothing
        if idx !== nothing
            i0 = idx
            # linear interpolation between i0 and i0+1
            denom = (λvals[i0 + 1] - λvals[i0])
            if denom == 0
                eps_cross = eps[i0]
            else
                eps_cross = eps[i0] + (target - λvals[i0]) * (eps[i0 + 1] - eps[i0]) / denom
            end
        else
            # fallback: nearest sampled eps
            j = findmin(abs.(λvals .- target))[2]
            eps_cross = eps[j]
        end

        eff_est[ch] = eps_cross
    end

    return eff_est
end

export estimate_efficiencies_from_curves
