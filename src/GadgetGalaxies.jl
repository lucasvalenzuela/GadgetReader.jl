__precompile__()
module GadgetGalaxies

using GadgetIO

include("snapshot.jl")
include("galaxy.jl")

# structs
export Snapshot,
       Galaxy

end
