"""
    x0_data(filename, runsel; ...) -> Table

Load `x0` values from LEGEND-200 SiPM data.

For each event and channel, this function records whether no photon
above the high photoelectron threshold was detected within the
`[-1, 5] μs` coincidence window. The output is a `Table` where each column
corresponds to a SiPM channel and each row to an event.

# Arguments
- `filename`: path to the *pygama* LH5 file containing event-tier data.
- `runsel`: `(period, run)` identifier, used to select the proper channel map.
- `max_events`: (optional) maximum number of events to read, default `10_000`.

# Returns
a `Table` of booleans with dimensions `(max_events, n_channels)`.
each entry is `true` if no qualifying photon was observed,
`false` otherwise.

# Examples
```julia
x0_data("l200-p13-r003-anp-20241217T094846Z-tier_evt.lh5", (:p13, :r003))
```
"""
function x0_data(filename::AbstractString, runsel::RunSelLike; max_events = 10_000)::Table
    period, run = runsel
    chmap = CHANNELMAPS[period][run]

    x0_cols = nothing

    lh5open(filename) do f
        # only interested in SiPM data
        events = f["evt"][1:max_events].spms

        # we want to check if there is light for each event and channel
        x0_cols = Dict(sipm => trues(first(size(events))) for sipm in keys(chmap))

        # map rawid directly to the column vector
        rawid2col = Dict(v.rawid => x0_cols[k] for (k, v) in chmap)

        # FIXME: this is quite slow?
        for i in eachindex(events)
            rawid = events[i].rawid
            # we use the high pe threshold to avoid counting afterpulses
            # this mask will also bring quality cuts and use the [-1, 5] us window
            mask = events[i].is_trig_coin_pulse_high_thr

            # loop over channels
            for j in eachindex(rawid)
                # check if there is any photon
                if any(mask[j])
                    rawid2col[rawid[j]][i] = false
                end
            end
        end
    end

    return Table(; x0_cols...)
end

export x0_data

"""
    λ0_data(x0; multiplicity_thr=0) -> NamedTuple

Compute per–channel no-light fractions from a boolean `Table` of events.

# Arguments
- `x0`: a `Table` where each column corresponds to a channel and each row
  to an event; entries are `Bool` indicating whether no qualifying photon was observed in the channel.
- `multiplicity_thr`: minimum number of channels with light per event required
  for the event to be considered. Defaults to `0`.

# Returns
A `NamedTuple` with one field per channel containing the fraction of
selected events in which that channel had no light, and the total number of events
passing the multiplicity threshold.
"""
function λ0_data(x0::Table; multiplicity_thr::Int = 0)::Tuple{NamedTuple,Integer}
    # compute multiplicity per row
    mult = zeros(Int, length(x0))
    for col in columns(x0)
        mult .+= .!col
    end

    keep = mult .>= multiplicity_thr
    nsel = count(keep)
    nsel == 0 && error("no events pass multiplicity_thr=$multiplicity_thr")

    # build NamedTuple of fractions
    return NamedTuple{columnnames(x0)}(map(col -> count(col[keep]) / nsel, columns(x0))), nsel
end

export λ0_data
