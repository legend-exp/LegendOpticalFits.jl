using BAT
using Distributions
using ValueShapes

"""
    make_efficiencies_prior(channels)

Construct a hierarchical prior for `channels` efficiencies.

- a global hyperprior `scale ~ Beta(7, 5)` sets the typical efficiency scale.
- each channel efficiency has a Beta prior with mean `scale * shift` and
  concentration `concentration`, enforcing positivity and preference for lower
  values.

# Arguments
- `channels`: A list of channel names.
"""
function make_efficiencies_prior(channels)
    # hyper prior that controls the global scale of the priors (e.g.
    # uncertainty on the light yield). it is a beta distribution located at 0.6
    # (i.e. ~30 ph/keV with a flat-top light yield of 51 ph/keV). the beta
    # distribution works well since the domain of the efficiencies is [0, 1]
    hyper_prior = NamedTupleDist(scale = Beta(7, 5))

    # the higher the concentration, the narrower the priors
    concentration = 7
    # shift towards zero of the prior location compared to the hyperprior
    # location
    shift = 0.6

    # the priors are also beta distributions, with location controlled by the
    # hyperprior
    return HierarchicalDistribution(
        v -> begin
            μ = v.scale * shift
            α = μ * concentration
            β = (1 - μ) * concentration
            NamedTupleDist(; (ch => Beta(α, β) for ch in channels)...)
        end,
        hyper_prior
    )
end

export make_efficiencies_prior
