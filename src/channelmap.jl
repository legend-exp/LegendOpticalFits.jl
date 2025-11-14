_chmap_p13_r001 = Dict(
    :S002 => (rawid = 1057604, usable = true, fiber = :IB07),
    :S003 => (rawid = 1057605, usable = true, fiber = :IB08),
    :S008 => (rawid = 1059205, usable = false, fiber = :OB20),
    :S011 => (rawid = 1062404, usable = true, fiber = :OB15),
    :S012 => (rawid = 1062405, usable = false, fiber = :OB16),
    :S013 => (rawid = 1054400, usable = false, fiber = :OB07),
    :S015 => (rawid = 1064004, usable = true, fiber = :OB33),
    :S017 => (rawid = 1054401, usable = true, fiber = :OB08),
    :S020 => (rawid = 1064000, usable = true, fiber = :OB29),
    :S023 => (rawid = 1057601, usable = true, fiber = :IB06),
    :S024 => (rawid = 1057600, usable = true, fiber = :IB05),
    :S025 => (rawid = 1064003, usable = true, fiber = :OB32),
    :S026 => (rawid = 1064002, usable = true, fiber = :OB31),
    :S027 => (rawid = 1059203, usable = false, fiber = :OB12),
    :S029 => (rawid = 1056003, usable = true, fiber = :IB04),
    :S030 => (rawid = 1057602, usable = true, fiber = :OB09),
    :S031 => (rawid = 1057603, usable = true, fiber = :OB10),
    :S032 => (rawid = 1059200, usable = true, fiber = :IB09),
    :S033 => (rawid = 1064005, usable = false, fiber = :OB34),
    :S035 => (rawid = 1067205, usable = false, fiber = :OB04),
    :S036 => (rawid = 1059201, usable = true, fiber = :IB10),
    :S037 => (rawid = 1067204, usable = true, fiber = :OB03),
    :S040 => (rawid = 1065601, usable = true, fiber = :OB14),
    :S041 => (rawid = 1056005, usable = true, fiber = :OB24),
    :S042 => (rawid = 1056004, usable = true, fiber = :OB23),
    :S043 => (rawid = 1065600, usable = true, fiber = :OB13),
    :S046 => (rawid = 1062402, usable = true, fiber = :OB17),
    :S047 => (rawid = 1062403, usable = true, fiber = :OB18),
    :S048 => (rawid = 1065602, usable = true, fiber = :OB35),
    :S049 => (rawid = 1065603, usable = true, fiber = :OB36),
    :S050 => (rawid = 1067200, usable = true, fiber = :OB39),
    :S051 => (rawid = 1067201, usable = true, fiber = :OB40),
    :S052 => (rawid = 1065605, usable = true, fiber = :OB38),
    :S053 => (rawid = 1065604, usable = true, fiber = :OB37),
    :S054 => (rawid = 1052805, usable = false, fiber = :OB06),
    :S055 => (rawid = 1052804, usable = true, fiber = :OB05),
    :S057 => (rawid = 1060801, usable = true, fiber = :IB12),
    :S058 => (rawid = 1060800, usable = false, fiber = :IB11),
    :S060 => (rawid = 1052802, usable = false, fiber = :IB15),
    :S061 => (rawid = 1052803, usable = true, fiber = :IB16),
    :S065 => (rawid = 1060805, usable = true, fiber = :OB26),
    :S067 => (rawid = 1056000, usable = false, fiber = :OB21),
    :S068 => (rawid = 1056001, usable = true, fiber = :OB22),
    :S070 => (rawid = 1054405, usable = true, fiber = :IB02),
    :S071 => (rawid = 1054404, usable = true, fiber = :IB01),
    :S073 => (rawid = 1054403, usable = true, fiber = :IB18),
    :S080 => (rawid = 1064001, usable = true, fiber = :OB30),
    :S082 => (rawid = 1062401, usable = true, fiber = :OB28),
    :S083 => (rawid = 1054402, usable = true, fiber = :IB17),
    :S085 => (rawid = 1067202, usable = false, fiber = :OB01),
    :S086 => (rawid = 1067203, usable = true, fiber = :OB02),
    :S087 => (rawid = 1062400, usable = true, fiber = :OB27),
    :S090 => (rawid = 1056002, usable = true, fiber = :IB03),
    :S094 => (rawid = 1059202, usable = true, fiber = :OB11),
    :S095 => (rawid = 1060803, usable = true, fiber = :IB14),
    :S096 => (rawid = 1060802, usable = true, fiber = :IB13),
    :S098 => (rawid = 1059204, usable = true, fiber = :OB19),
    :S099 => (rawid = 1060804, usable = true, fiber = :OB25)
)

# copy r001
_chmap_p13_r002 = copy(_chmap_p13_r001)

# update existing keys and add new ones
merge!(
    _chmap_p13_r002,
    Dict(
        # value changes
        :S008 => (rawid = 1059205, usable = true, fiber = :OB20),
        :S012 => (rawid = 1062405, usable = true, fiber = :OB16),
        :S024 => (rawid = 1057600, usable = false, fiber = :IB05),
        :S058 => (rawid = 1060800, usable = true, fiber = :IB11),
        :S080 => (rawid = 1064001, usable = false, fiber = :OB30),

        # additions
        :S007 => (rawid = 1067205, usable = true, fiber = :OB04),
        :S028 => (rawid = 1054400, usable = false, fiber = :OB07),
        :S097 => (rawid = 1064005, usable = false, fiber = :OB34)
    )
)

# remove keys that are not present in r002
for k in (:S013, :S033, :S035)
    pop!(_chmap_p13_r002, k, nothing)
end

_chmap_p14_r006 = Dict(
    :S002 => (rawid = 1057604, usable = true, fiber = :IB07),
    :S003 => (rawid = 1057605, usable = true, fiber = :IB08),
    :S007 => (rawid = 1067205, usable = true, fiber = :OB04),
    :S008 => (rawid = 1059205, usable = true, fiber = :OB20),
    :S011 => (rawid = 1062404, usable = true, fiber = :OB15),
    :S012 => (rawid = 1062405, usable = true, fiber = :OB16),
    :S015 => (rawid = 1064004, usable = true, fiber = :OB33),
    :S017 => (rawid = 1054401, usable = true, fiber = :OB08),
    :S020 => (rawid = 1064000, usable = true, fiber = :OB29),
    :S023 => (rawid = 1057601, usable = true, fiber = :IB06),
    :S024 => (rawid = 1057600, usable = false, fiber = :IB05),
    :S025 => (rawid = 1064003, usable = true, fiber = :OB32),
    :S026 => (rawid = 1064002, usable = true, fiber = :OB31),
    :S027 => (rawid = 1059203, usable = false, fiber = :OB12),
    :S028 => (rawid = 1054400, usable = false, fiber = :OB07),
    :S029 => (rawid = 1056003, usable = true, fiber = :IB04),
    :S030 => (rawid = 1057602, usable = true, fiber = :OB09),
    :S031 => (rawid = 1057603, usable = true, fiber = :OB10),
    :S032 => (rawid = 1059200, usable = true, fiber = :IB09),
    :S036 => (rawid = 1059201, usable = true, fiber = :IB10),
    :S037 => (rawid = 1067204, usable = true, fiber = :OB03),
    :S040 => (rawid = 1065601, usable = true, fiber = :OB14),
    :S041 => (rawid = 1056005, usable = true, fiber = :OB24),
    :S042 => (rawid = 1056004, usable = true, fiber = :OB23),
    :S043 => (rawid = 1065600, usable = true, fiber = :OB13),
    :S046 => (rawid = 1062402, usable = true, fiber = :OB17),
    :S047 => (rawid = 1062403, usable = true, fiber = :OB18),
    :S048 => (rawid = 1065602, usable = true, fiber = :OB35),
    :S049 => (rawid = 1065603, usable = true, fiber = :OB36),
    :S050 => (rawid = 1067200, usable = true, fiber = :OB39),
    :S051 => (rawid = 1067201, usable = true, fiber = :OB40),
    :S052 => (rawid = 1065605, usable = true, fiber = :OB38),
    :S053 => (rawid = 1065604, usable = true, fiber = :OB37),
    :S054 => (rawid = 1052805, usable = false, fiber = :OB06),
    :S055 => (rawid = 1052804, usable = true, fiber = :OB05),
    :S057 => (rawid = 1060801, usable = true, fiber = :IB12),
    :S058 => (rawid = 1060800, usable = false, fiber = :IB11),
    :S060 => (rawid = 1052802, usable = false, fiber = :IB15),
    :S061 => (rawid = 1052803, usable = true, fiber = :IB16),
    :S065 => (rawid = 1060805, usable = true, fiber = :OB26),
    :S067 => (rawid = 1056000, usable = true, fiber = :OB21),
    :S068 => (rawid = 1056001, usable = true, fiber = :OB22),
    :S070 => (rawid = 1054405, usable = true, fiber = :IB02),
    :S071 => (rawid = 1054404, usable = true, fiber = :IB01),
    :S073 => (rawid = 1054403, usable = true, fiber = :IB18),
    :S080 => (rawid = 1064001, usable = true, fiber = :OB30),
    :S082 => (rawid = 1062401, usable = true, fiber = :OB28),
    :S083 => (rawid = 1054402, usable = true, fiber = :IB17),
    :S085 => (rawid = 1067202, usable = true, fiber = :OB01),
    :S086 => (rawid = 1067203, usable = true, fiber = :OB02),
    :S087 => (rawid = 1062400, usable = true, fiber = :OB27),
    :S090 => (rawid = 1056002, usable = false, fiber = :IB03),
    :S094 => (rawid = 1059202, usable = true, fiber = :OB11),
    :S095 => (rawid = 1060803, usable = true, fiber = :IB14),
    :S096 => (rawid = 1060802, usable = false, fiber = :IB13),
    :S097 => (rawid = 1064005, usable = false, fiber = :OB34),
    :S098 => (rawid = 1059204, usable = true, fiber = :OB19),
    :S099 => (rawid = 1060804, usable = true, fiber = :OB25)
)


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
    ),
    :p14 => Dict(
        :r006 => _chmap_p14_r006
    )
)

# NOTE: don't use this in hot loops, it's slow!
function rawid2detname(chmap::AbstractDict{Symbol,<:NamedTuple}, rawid::Integer)
    for (k, v) in chmap
        v.rawid == rawid && return k
    end
    return nothing
end
