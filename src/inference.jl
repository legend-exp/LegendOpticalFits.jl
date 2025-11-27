# here routines to perform inference on the optical model parameters
"""
    compute_λ0_curve(log_p0_nominal, x0_random_coin, channel; eps_range=range(0.0,1.0,length=101), multiplicity_thr=0, use_mc=false, n_rands=10)

Compute the per-channel no-light probability λ₀ as a function of a single
channel's efficiency.

# Behavior modes
- Analytic (default, use_mc=false): computes the exact expectation
  per event λ₀_jk = mean(p0_jk * x0_jk) with p0_jk = exp(log_p0_jk * ε) and x0_jk 
  the rando coincidence indicators; this is exact when multiplicity_thr == 0.
- Monte‑Carlo (use_mc=true): draws n_rands Bernoulli samples per event
  from Bernoulli(p0) and averages sampled outcomes folded with x0_random_coin. 
  Use this when you want to mimic the forward-model sampling.

# Arguments
- `log_p0_nominal::Table`: table of log p0 values (events × channels).
- `x0_random_coin::Table`: table of random-coincidence indicators (events × channels).
- `channel::Symbol`: Symbol identifying the channel (must be present in
  columnnames(log_p0_nominal)).
- `eps_range`: efficiency values to evaluate (vector or range).
- `multiplicity_thr`: must be 0, otherwise channel correlations arise (not an actual parameter).
- `use_mc`: whether to use Monte‑Carlo sampling (default false).
- `n_rands`: number of Monte‑Carlo draws when use_mc=true.

# Returns
- `NamedTuple` with fields :eps (vector of efficiencies) and :λ0
  (vector of no-light probabilities evaluated at those efficiencies).
"""
function compute_λ0_curve(
    log_p0_nominal::Table,
    x0_random_coin::Table,
    channel::Symbol;
    eps_range = range(0.0, 1.0, length = 101),
    multiplicity_thr::Int = 0,
    use_mc::Bool = false,
    n_rands::Int = 10
)
    if multiplicity_thr != 0
        throw(ArgumentError("λ0_vs_efficiency assumes multiplicity_thr == 0 (channels independent)."))
    end

    # derive channel ordering from the provided `log_p0_nominal` table
    order = collect(columnnames(log_p0_nominal))
    log_p0, _ = _to_matrix(log_p0_nominal, order = order)
    x0, _ = _to_matrix(x0_random_coin, order = order)

    # resolve channel index (channel must be a Symbol and present in log_p0_nominal)
    ch_idx = findfirst(==(channel), order)
    if ch_idx === nothing
        throw(ArgumentError("channel $channel not found in efficiencies keys"))
    end

    # columns for the requested channel
    log_p0_col = view(log_p0, :, ch_idx)
    x0_col = view(x0, :, ch_idx)

    eps_vec = collect(eps_range)
    n = length(eps_vec)
    λ0s = Vector{Float64}(undef, n)
    for (i, ε) in enumerate(eps_vec)
        p0 = exp.(log_p0_col .* ε)
        if use_mc
            # Monte-Carlo sampling: draw Bernoulli samples and fold in random coincidences
            n_events = length(p0)
            rng = default_device_rng(get_device(p0))
            rands = rand(rng, n_events, n_rands)
            # (rands .< p0) broadcasts p0 over columns; x0_col broadcasts similarly
            drawn = (rands .< p0) .* x0_col
            λ0s[i] = mean(drawn)
        else
            # analytic exact expectation for multiplicity-free case
            λ0s[i] = mean(p0 .* x0_col)
        end
    end

    return (eps = eps_vec, λ0 = λ0s)
end

"""
    compute_λ0_curve_all(log_p0_nominal, x0_random_coin; eps_range=range(0.0,1.0,length=101), multiplicity_thr=0, use_mc=false, n_rands=10)

Compute λ₀(eps) curves for every channel present in `log_p0_nominal`.

Description
- For each channel (columns of `log_p0_nominal`) evaluate the no-light
  probability λ₀ across `eps_range` by calling `compute_λ0_curve`.
- The returned curves are suitable for locating the efficiency at which a
  measured λ₀ would be reproduced (e.g. by interpolation).

Arguments
- see above in description of `compute_λ0_curve`

Returns
- `Dict{Symbol, NamedTuple}` mapping each channel symbol to the NamedTuple
  returned by `compute_λ0_curve` (fields `:eps` and `:λ0`).

Notes
- Channel order is derived from `columnnames(log_p0_nominal)`. Use that
  table to control which channels are computed and their ordering.

"""
function compute_λ0_curve_all(
    log_p0_nominal::Table,
    x0_random_coin::Table;
    eps_range = range(0.0, 1.0, length = 101),
    multiplicity_thr::Int = 0,
    use_mc::Bool = false,
    n_rands::Int = 10
)
    if multiplicity_thr != 0
        throw(ArgumentError("compute_λ0_curve_all assumes multiplicity_thr == 0 (channels independent)."))
    end

    # derive channel list from the table column names
    ch_list = collect(columnnames(log_p0_nominal))
    results = Dict{Symbol,Any}()
    for ch in ch_list
        results[ch] = compute_λ0_curve(
            log_p0_nominal,
            x0_random_coin,
            ch;
            eps_range = eps_range,
            multiplicity_thr = multiplicity_thr,
            use_mc = use_mc,
            n_rands = n_rands
        )
    end

    return results
end

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
