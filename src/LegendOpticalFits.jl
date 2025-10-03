module LegendOpticalFits

using CSV
using DensityInterface
using Distributions
using HDF5
using HeterogeneousComputing
using LegendDataManagement: RunSelLike
using LegendHDF5IO
using MLDataDevices
using Pkg.Artifacts
using Random
using StatsBase
using TypedTables
using Tables
using Unitful

include("utils.jl")
include("channelmap.jl")
include("data.jl")
include("optmap.jl")
include("simulation.jl")
include("models.jl")
include("likelihood.jl")
include("priors.jl")
include("inference.jl")

end
