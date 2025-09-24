# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using LegendOpticalFits

# Doctest setup
DocMeta.setdocmeta!(
    LegendOpticalFits,
    :DocTestSetup,
    :(using LegendOpticalFits);
    recursive = true
)

makedocs(
    sitename = "LegendOpticalFits",
    modules = [LegendOpticalFits],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/LegendOpticalFits.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md"
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = false,
    warnonly = ("nonstrict" in ARGS)
)

deploydocs(
    repo = "github.com/legend-exp/LegendOpticalFits.jl.git",
    forcepush = true,
    push_preview = true
)
