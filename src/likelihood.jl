"""
    make_likelihood_λ0()

Construct the likelihood of no-light probabilities.

We model the fraction of events with no detected light (`λ0`) as follows:

- For each channel, the expected no-light probability `λ0_model` comes from
  the simulation (`log_p0_nominal`) combined with random coincidences
  (`x0_random_coin`) and scaled by per-channel scaling factors
  (the parameters of the model).

- The observed no-light fraction in data is `λ0_data = N0 / N_data`, where `N0`
  is the number of no-light events passing a multiplicity threshold.

- Since `N_data` is large, the binomial distribution `N0 ~ Binomial(N_data,
  λ0_model)` can be approximated by a normal distribution: `λ0_data ~ Normal(μ =
  λ0_model, σ² = λ0_model (1 - λ0_model) / N_data)`.

The likelihood is the sum of log-probabilities across all channels.

# Arguments
- `x0`: observed no-light indicators from data events.
- `log_p0_nominal`: logarithm of the probability to see no light for each
  (event, channel), typically from simulations.
- `x0_random_coin`: observed no-light indicators from random coincidence events.
- `multiplicity_thr`: discard events with multiplicity below this threshold
  (optional, defaults to 0).
- `smear_factor`: the width of the likelihood gaussian terms is increased by a
  factor `smear_factor * mean`.

# Returns
- A `DensityFunction` object representing the log-likelihood. It can be called
  with a parameter set.
"""
function make_λ0_likelihood(
    x0::Table,
    log_p0_nominal::Table,
    x0_random_coin::Table,
    ;
    multiplicity_thr::Int = 0,
    smear_factor::Real = 0
)
    # we choose as channel order the one used in x0
    ϵ_order = columnnames(x0)

    N_data = length(x0)
    data = λ0_data(x0)

    # convert to matrix with the correct order
    log_p0, _ = _to_matrix(log_p0_nominal, order = ϵ_order)
    x0_rc, _ = _to_matrix(x0_random_coin, order = ϵ_order)

    return DensityInterface.logfuncdensity(
        # params is expected to be a NamedTuple, we just pass the values to the
        # low level routines
        params -> begin
            # make sure the order of the parameters is the correct one
            p = [params[k] for k in ϵ_order]

            # compute the forward model
            model = λ0_model(p, log_p0, x0_rc, multiplicity_thr = multiplicity_thr)

            logpmf = 0.0
            @inbounds @simd for i in eachindex(model)
                # Binomal statistics
                µ = model[i]
                σ = sqrt(µ * (1 - µ) / N_data) + smear_factor * μ

                logpmf += logpdf(Normal(μ, σ), data[i])
            end

            return logpmf
        end
    )
end

export make_λ0_likelihood
