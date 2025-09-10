"""
Construct the likelihood function for fitting no-light probabilities.

We model the fraction of events with no detected light (`λ₀`) as follows:

- For each channel, the expected no-light probability `λ₀_model` comes from
  the simulation (`log_p0_nominal`) combined with random coincidences (`x0_random_coin`)
  and scaled by per-channel scaling factors (`params.epsilons`).

- The observed no-light fraction in data is
  `λ₀_data = N₀ / N_data`, where `N₀` is the number of no-light events
  passing a multiplicity threshold.

- Since `N_data` is large, the binomial distribution
  `N₀ ~ Binomial(N_data, λ₀_model)`
  can be approximated by a normal distribution:
  λ₀_data ~ Normal(μ = λ₀_model, σ² = λ₀_model (1 - λ₀_model) / N_data).

The likelihood is the sum of log-probabilities across all channels.

# Arguments
- `N_data`: number of selected data events after applying the multiplicity cut.
- `λ0_data`: observed no-light fractions per channel from data.
- `log_p0_nominal`: logarithm of the probability to see no light for each
  (event, channel), typically from simulations.
- `x0_random_coin`: presence of light from random coincidences for each
  (event, channel). This is typically coming from a measurement.
- `multiplicity_thr`: discard events with multiplicity below this threshold
  (optional, defaults to 0).

# Returns
- A `DensityFunction` (from BAT.jl) representing the log-likelihood.
  It can be called with a parameter set `params` that must include
  `params.epsilons`.
  
"""
function make_likelihood(
    N_data::Int,
    λ0_data::AbstractMatrix{<:Float64},
    log_p0_nominal::AbstractMatrix{<:Float64},
    x0_random_coin::AbstractMatrix{Bool},
    ;
    multiplicity_thr::Int = 8
)    
    return logfuncdensity(params -> begin

        λ0_model = expected_no_light_fraction(params.epsilons, log_p0_nominal, x0_random_coin; multiplicity_thr=multiplicity_thr)

        logpmf = 0.0
        @inbounds @simd for i in eachindex(λ0_model)
            µ = λ0_model[i]
            σ = max(sqrt(µ * (1 - µ) / N_data), 1e-6)

            logpmf += logpdf(Normal(μ, σ), λ0_data[i])
        end

        return logpmf
    end)
end