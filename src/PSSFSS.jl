module PSSFSS

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 3
end

include("Constants.jl")
include("PSSFSSLen.jl")
include("Rings.jl")
include("Layers.jl")
include("Sheets.jl")
include("Meshsub.jl")
include("Elements.jl")
include("RWG.jl")
include("PGF.jl")
include("Zint.jl")
include("FillZY.jl")
include("GSMs.jl")
include("Modes.jl")

using Reexport
using .Rings
@reexport using .PSSFSSLen
@reexport using .Layers
@reexport using .Sheets: RWGSheet
@reexport using .Elements: rectstrip, polyring, meander, loadedcross, jerusalemcross
@reexport using .RWG: RWGData, setup_rwg, edge_current_unit_vector
@reexport using .PGF: electric_modal_sum_funcs, magnetic_modal_sum_funcs
@reexport using .GSMs: GSM, cascade, cascade!, gsm_electric_gblock,
                      gsm_slab_interface, initialize_gsm_file, translate_gsm!
@reexport using .FillZY: fillz, filly

end
