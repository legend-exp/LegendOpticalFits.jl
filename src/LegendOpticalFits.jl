module LegendOpticalFits

using CSV
using Distributions
using HDF5
using LegendDataManagement: RunSelLike
using LegendHDF5IO
using Pkg.Artifacts
using Random
using StatsBase

include("utils.jl")
include("channelmap.jl")
include("optmap.jl")
include("models.jl")

end
