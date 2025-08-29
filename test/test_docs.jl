# This file is a part of LegendOpticalFits.jl, licensed under the MIT License (MIT).

using Test
using LegendOpticalFits
import Documenter

Documenter.DocMeta.setdocmeta!(
    LegendOpticalFits,
    :DocTestSetup,
    :(using LegendOpticalFits);
    recursive = true
)
Documenter.doctest(LegendOpticalFits)
