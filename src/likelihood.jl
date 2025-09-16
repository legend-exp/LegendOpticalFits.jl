"""
    make_likelihood_λ0()

Construct the likelihood of no-light probabilities.

We model the fraction of events with no detected light (`λ0`) as follows:

- For each channel, the expected no-light probability `λ0_model` comes from
  the simulation (`log_p0_nominal`) combined with random coincidences
  (`x0_random_coin`) and scaled by per-channel scaling factors
  (`params.epsilons`, the parameters of the model).

- The observed no-light fraction in data is `λ0_data = N0 / N_data`, where `N0`
  is the number of no-light events passing a multiplicity threshold.

- Since `N_data` is large, the binomial distribution `N0 ~ Binomial(N_data,
  λ0_model)` can be approximated by a normal distribution: `λ0_data ~ Normal(μ =
  λ0_model, σ² = λ0_model (1 - λ0_model) / N_data)`.

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
- A `DensityFunction` object representing the log-likelihood.  It can be called
  with a parameter set `params` that must include `params.epsilons`
"""
function make_λ0_likelihood(
    N_data::Int,
    λ0_data::Table,
    log_p0_nominal::Table,
    x0_random_coin::Table,
    ;
    multiplicity_thr::Int = 8
)
    return logfuncdensity(
        params -> begin

            λ0_model = expected_λ0(
                params.epsilons,
                log_p0_nominal,
                x0_random_coin;
                multiplicity_thr = multiplicity_thr
            )

            logpmf = 0.0
            @inbounds @simd for i in eachindex(λ0_model)
                # Binomal statistics
                µ = λ0_model[i]
                σ = sqrt(µ * (1 - µ) / N_data)

                logpmf += logpdf(Normal(μ, σ), λ0_data[i])
            end

            return logpmf
        end
    )
end
