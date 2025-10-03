_chmap_p13_r001 = Dict(
    :S002 => (rawid = 1057604, usable = true),
    :S003 => (rawid = 1057605, usable = true),
    :S008 => (rawid = 1059205, usable = false),
    :S011 => (rawid = 1062404, usable = true),
    :S012 => (rawid = 1062405, usable = false),
    :S013 => (rawid = 1054400, usable = false),
    :S015 => (rawid = 1064004, usable = true),
    :S017 => (rawid = 1054401, usable = true),
    :S020 => (rawid = 1064000, usable = true),
    :S023 => (rawid = 1057601, usable = true),
    :S024 => (rawid = 1057600, usable = true),
    :S025 => (rawid = 1064003, usable = true),
    :S026 => (rawid = 1064002, usable = true),
    :S027 => (rawid = 1059203, usable = false),
    :S029 => (rawid = 1056003, usable = true),
    :S030 => (rawid = 1057602, usable = true),
    :S031 => (rawid = 1057603, usable = true),
    :S032 => (rawid = 1059200, usable = true),
    :S033 => (rawid = 1064005, usable = false),
    :S035 => (rawid = 1067205, usable = false),
    :S036 => (rawid = 1059201, usable = true),
    :S037 => (rawid = 1067204, usable = true),
    :S040 => (rawid = 1065601, usable = true),
    :S041 => (rawid = 1056005, usable = true),
    :S042 => (rawid = 1056004, usable = true),
    :S043 => (rawid = 1065600, usable = true),
    :S046 => (rawid = 1062402, usable = true),
    :S047 => (rawid = 1062403, usable = true),
    :S048 => (rawid = 1065602, usable = true),
    :S049 => (rawid = 1065603, usable = true),
    :S050 => (rawid = 1067200, usable = true),
    :S051 => (rawid = 1067201, usable = true),
    :S052 => (rawid = 1065605, usable = true),
    :S053 => (rawid = 1065604, usable = true),
    :S054 => (rawid = 1052805, usable = false),
    :S055 => (rawid = 1052804, usable = true),
    :S057 => (rawid = 1060801, usable = true),
    :S058 => (rawid = 1060800, usable = false),
    :S060 => (rawid = 1052802, usable = false),
    :S061 => (rawid = 1052803, usable = true),
    :S065 => (rawid = 1060805, usable = true),
    :S067 => (rawid = 1056000, usable = false),
    :S068 => (rawid = 1056001, usable = true),
    :S070 => (rawid = 1054405, usable = true),
    :S071 => (rawid = 1054404, usable = true),
    :S073 => (rawid = 1054403, usable = true),
    :S080 => (rawid = 1064001, usable = true),
    :S082 => (rawid = 1062401, usable = true),
    :S083 => (rawid = 1054402, usable = true),
    :S085 => (rawid = 1067202, usable = false),
    :S086 => (rawid = 1067203, usable = true),
    :S087 => (rawid = 1062400, usable = true),
    :S090 => (rawid = 1056002, usable = true),
    :S094 => (rawid = 1059202, usable = true),
    :S095 => (rawid = 1060803, usable = true),
    :S096 => (rawid = 1060802, usable = true),
    :S098 => (rawid = 1059204, usable = true),
    :S099 => (rawid = 1060804, usable = true)
)

# copy r001
_chmap_p13_r002 = copy(_chmap_p13_r001)

# update existing keys and add new ones
merge!(
    _chmap_p13_r002,
    Dict(
        # value changes
        :S008 => (rawid = 1059205, usable = true),
        :S012 => (rawid = 1062405, usable = true),
        :S024 => (rawid = 1057600, usable = false),
        :S058 => (rawid = 1060800, usable = true),
        :S080 => (rawid = 1064001, usable = false),

        # additions
        :S007 => (rawid = 1067205, usable = true),
        :S028 => (rawid = 1054400, usable = false),
        :S097 => (rawid = 1064005, usable = false)
    )
)

# remove keys that are not present in r002
for k in (:S013, :S033, :S035)
    pop!(_chmap_p13_r002, k, nothing)
end

CHANNELMAPS = Dict(
    :p13 => Dict(
        :r001 => _chmap_p13_r001,
        :r002 => _chmap_p13_r002,
        :r003 => _chmap_p13_r002,
        :r004 => _chmap_p13_r002,
        :r005 => _chmap_p13_r002,
        :r006 => _chmap_p13_r002,
        :r008 => _chmap_p13_r002,
        :r009 => _chmap_p13_r002
    )
)

# NOTE: don't use this in hot loops, it's slow!
function rawid2detname(chmap::AbstractDict{Symbol,<:NamedTuple}, rawid::Integer)
    for (k, v) in chmap
        v.rawid == rawid && return k
    end
    return nothing
end
