"""
    x0_mock(efficiencies_true, log_p0_nominal, x0_random_coin)

Generate mock x0 data (bernoulli draws folded with random-coincidence mask)
from `efficiencies_true` and `log_p0_nominal` from model.
"""
function x0_mock(
    efficiencies_true::NamedTuple,
    log_p0_nominal::Table,
    x0_random_coin::Table
)
    # make sure order is consistent with provided scaling_factors
    ϵk = keys(efficiencies_true)
    ϵv = collect(values(efficiencies_true))
    log_p0, _ = _to_matrix(log_p0_nominal, order = ϵk)
    x0_rc, _ = _to_matrix(x0_random_coin, order = ϵk)

    # get random engine depending on device (CPU/GPU)
    rng = default_device_rng(get_device(log_p0))
    # pre-allocate random numbers for forward model evaluation (2D)
    n_events, n_channels = size(log_p0)
    rands = rand(rng, n_events, n_channels)

    # avoid numerical issues
    ϵ = clamp.(ϵv, 1e-10, 1 - 1e-10)

    # calculate the probability to see no light (duplicate of bulk logic)
    p0 = exp.(log_p0 .* ϵ')

    # draw bernoulli distributed numbers and fold in random coincidences
    drawn = (rands .< p0) .* x0_rc

    # return as a Table 
    return Table(; (ϵk[i] => drawn[:, i] for i in 1:size(drawn, 2))...)
end

export x0_mock
