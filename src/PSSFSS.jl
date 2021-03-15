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
include("Outputs.jl")

using Reexport
using LinearAlgebra: ×, norm, ⋅, factorize
using StaticArrays: MVector, MArray
using .Rings
@reexport using .PSSFSSLen
@reexport using .Layers: Layer, TEorTM, TE, TM
@reexport using .Elements: rectstrip, polyring, meander, loadedcross, jerusalemcross
@reexport using .Outputs
using .Sheets: Sheet, RWGSheet
using .RWG: setup_rwg, rwgbfft!, RWGData
using .GSMs: GSM, cascade, cascade, gsm_electric_gblock, gsm_magnetic_gblock,
                      gsm_slab_interface, initialize_gsm_file, translate_gsm!
using .FillZY: fillz, filly
using .Modes: zhatcross
using .Constants: twopi


const tdigits = 4 # Number of decimal places used to display elapsed time

"""
    calculate_jtype_gsm(layers, sheet::RWGSheet, u::Real, rwgdat::RWGData, s::Int, k0, k⃗inc, is_global::Int) -> gsm

Compute the generalized scattering matrix for a sheet of class `'J'`.

### Input Arguments

- `layers`: An iterable of `Layer` instances containing the layers for the `Gblock`
    associated with the `sheet` under consideration.
- `sheet`:  A sheet of class `'J'` for which the GSM is desired.
- `u`: The Green's function smoothing parameter for the `sheet`.
- `rwgdat`: The `RWGData` object associated with the `sheet` argument.
- `s`: The interface number within `layers` at which the sheet is located.
- `k0`: The free-space wavenumber in radians/meter.
- `k⃗inc`: A 2-vector containing the incident field wave vector x and y components. Note
    that by the phase match condition this vector is the same for all layers in the entire FSS
    structure.
- `is_global`: The global sheet index for the `sheet` argument within the global list of sheets.

### Return Value

- `gsm::GSM`  The full GSM for the GBlock including incident fields and scattered fields due
    to currents induced on the sheet surface.
"""
function calculate_jtype_gsm(layers, sheet::RWGSheet, u::Real,
                                 rwgdat::RWGData, s::Int, k0, k⃗inc, is_global::Int)
    one_meter = ustrip(Float64, sheet.units, 1u"m")
    area = norm(sheet.s₁ × sheet.s₂) / one_meter^2 # Unit cell area (m^2).
    acf = MVector(0.0,0.0)
    nmodesmax = max(length(layers[begin].P), length(layers[end].P)) 
    nbf = size(rwgdat.bfe,2) # Number of basis functions
    bfftstore = zeros(MArray{Tuple{2},ComplexF64,1,2}, (nbf, 2, nmodesmax))

    # Compute area correction factors for the mode normalization constants of 
    # the two end regions:
    for (i,l) in enumerate(@view layers[[begin,end]])
      area_i = twopi * twopi / norm(l.β₁ × l.β₂)
      acf[i] = √(area_i / area)
    end

    # Set up the partial GSM due to incident field
    (gsm, tlgfvi, vincs) = gsm_electric_gblock(layers, s, k0)
    if sheet.style == "NULL"
        return gsm
    end

    # Calculate the scattered field partial scattering matrix of the 
    # FSS sheet at this junction. Then add it to GSM already computed 
    # for the dielectric discontinuity...
    
    #  Fill the interaction matrix for the current sheet:
    t_temp = time()
    ψ₁ = k⃗inc ⋅ sheet.s₁ / one_meter
    ψ₂ = k⃗inc ⋅ sheet.s₂ / one_meter
    @info "  Beginning matrix fill for sheet $(is_global)"
    zmat = fillz(k0,u,layers,s,ψ₁,ψ₂,sheet,rwgdat)
    t_fill = round(time() - t_temp, digits=tdigits)
    @info "    Total matrix fill time for sheet $(is_global): $(t_fill) sec"
    # Factor the matrix:
    t_temp = time()
    zmatf = factorize(zmat)
    t_factor = round(time() - t_temp, digits=tdigits)
    @info "  Matrix factor for sheet $(is_global) used $(t_factor) sec"
    t_temp = time()
    # Compute and store the basis function Fourier transforms:
    i_ft = 0
    for (sr,l) in enumerate(@view layers[[begin,end]]) # loop over possible source regions
        for qp = 1:length(l.P) 
            kvec = l.β[qp]
            # If desired F.T. has already been computed, then copy it.
            if qp > 1 && kvec ≈ l.β[qp-1]
                bfftstore[:,sr,qp] = bfftstore[:,sr,qp-1] 
                continue
            end
            if sr == 2 && length(layers[1].β) ≥ qp && kvec ≈ layers[1].β[qp]
                bfftstore[:,2,qp] = bfftstore[:,1,qp] 
                continue
            end
            bfft = @view bfftstore[:,sr,qp]
            rwgbfft!(bfft, rwgdat, sheet, kvec, ψ₁, ψ₂) # Otherwise, compute from scratch
            i_ft += 1
        end
    end
    t_fft = round(time() - t_temp, digits=tdigits)
    @info "  Basis function Fourier Transforms at $(i_ft) points used $(t_fft) sec"
    nsolve = 0
    t_extract = 0.0
    i_extract = 0
    t_solve = 0.0
    for (sr,ls) in enumerate(@view layers[[begin,end]]) # Loop over source regions
        for qp in 1:length(ls.P) # Loop over srce reg modes
            # Incident field for source layer in absence of the FSS sheet:
            sourcevec = vincs[qp,sr] * ls.c[qp] * acf[sr] * ls.tvec[qp]
            # Compute generalized voltage vector:
            imat = [b ⋅ sourcevec for b in bfftstore[:,sr,qp]] # Eq. (7.39)
            # Solve the matrix equation
            t_solve1 = time()
            imat = zmatf \ imat
            t_solve2 = time()
            t_solve += t_solve2 - t_solve1
            nsolve += 1
            t_extract1 = time()
            for (or,lo) in enumerate(@view layers[[begin,end]]) # Loop over obs. regions
                smat = gsm[or,sr]
                for q in 1:length(lo.P)  # Loop obs. regn modes
                    # Extract partial scattering parameter due to scattered fields...
                    FTJ = sum((imat[n] * bfftstore[n,or,q] for n in 1:nbf)) # FT of total current
                    smat[q,qp] -= (lo.tvec[q] ⋅ FTJ) * (tlgfvi[q,or] / 
                                                (lo.c[q] * acf[or] * area)) # Eq (6.18)
                    i_extract += 1
                end
            end
            t_extract2 = time()
            t_extract = t_extract + (t_extract2 - t_extract1)
        end
    end
    @info "  Extracting $(i_extract) GSM entries used $(t_extract) sec"
    return gsm
end  

"""
    calculate_mtype_gsm(layers, sheet::RWGSheet, u::Real, rwgdat::RWGData, s::Int, k⃗inc, is_global::Int) -> gsm

Compute the generalized scattering matrix for a sheet of class `'M'`.

### Input Arguments

- `layers`: An iterable of `Layer` instances containing the layers for the `Gblock`
    associated with the `sheet` under consideration.
- `sheet`:  A sheet of class `'M'` for which the GSM is desired.
- `u`: The Green's function smoothing parameter for the `sheet`.
- `rwgdat`: The `RWGData` object associated with the `sheet` argument.
- `s`: The interface number within `layers` at which the sheet is located.
- `k0`: The free-space wavenumber in radians/meter.
- `k⃗inc`: A 2-vector containing the incident field wave vector x and y components. Note
    that by the phase match condition this vector is the same for all layers in the entire FSS
    structure.
- `is_global`: The global sheet index for the `sheet` argument within the global list of sheets.

### Return Value

- `gsm::GSM`  The full GSM for the GBlock including incident fields and scattered fields due
    to magnetic currents induced in the gaps on the sheet surface.
"""
function calculate_mtype_gsm(layers, sheet::RWGSheet, u::Real,
                                 rwgdat::RWGData, s::Int, k0, k⃗inc, is_global::Int)
    one_meter = ustrip(Float64, sheet.units, 1u"m")
    area = norm(sheet.s₁ × sheet.s₂) / one_meter^2 # Unit cell area (m^2).
    acf = MVector(0.0,0.0)
    nmodesmax = max(length(layers[begin].P), length(layers[end].P)) 
    nbf = size(rwgdat.bfe,2) # Number of basis functions
    bfftstore = zeros(MArray{Tuple{2},ComplexF64,1,2}, (nbf, 2, nmodesmax))

    # Compute area correction factors for the mode normalization constants of 
    # the two end regions:
    for (i,l) in enumerate(@view layers[[begin,end]])
      area_i = twopi * twopi / norm(l.β₁ × l.β₂)
      acf[i] = √(area_i / area)
    end

    # Set up the partial GSM due to incident field
    (gsm, tlgfiv, iincs) = gsm_magnetic_gblock(layers, s, k0)
    if sheet.style == "NULL"
        return gsm
    end

    # Calculate the scattered field partial scattering matrix of the 
    # FSS sheet at this junction. Then add it to GSM already computed 
    # for the dielectric discontinuity...
    
    #  Fill the interaction matrix for the current sheet:
    t_temp = time()
    ψ₁ = k⃗inc ⋅ sheet.s₁ / one_meter
    ψ₂ = k⃗inc ⋅ sheet.s₂ / one_meter
    @info "  Beginning matrix fill for sheet $(is_global)"
    ymat = filly(k0,u[is_global],layers,s,ψ₁,ψ₂,sheet,rwgdat)
    t_fill = round(time() - t_temp, digits=tdigits)
    @info "    Total matrix fill time for sheet $(is_global): $(t_fill) sec"
    # Factor the matrix:
    t_temp = time()
    ymatf = factorize(ymat)
    t_factor = round(time() - t_temp, digits=tdigits)
    @info "  Matrix factor for sheet $(is_global) used $(t_factor) sec"
    t_temp = time()
    # Compute and store the basis function Fourier transforms:
    i_ft = 0
    for (sr,l) in enumerate(@view layers[[begin,end]]) # loop over possible source regions
        for qp = 1:length(l.P) 
            kvec = l.β[qp]
            # If desired F.T. has already been computed, then copy it.
            if qp > 1 && kvec ≈ l.β[qp-1]
                bfftstore[:,sr,qp] = bfftstore[:,sr,qp-1] 
                continue
            end
            if sr == 2 && length(layers[1].β) ≥ qp && kvec ≈ layers[1].β[qp]
                bfftstore[:,2,qp] = bfftstore[:,1,qp] 
                continue
            end
            bfft = @view bfftstore[:,sr,qp]
            rwgbfft!(bfft, rwgdat, sheet, kvec, ψ₁, ψ₂) # Otherwise, compute from scratch
            i_ft += 1
        end
    end
    t_fft = round(time() - t_temp, digits=tdigits)
    @info "  Basis function Fourier Transforms at $(i_ft) points used $(t_fft) sec"
    nsolve = 0
    t_extract = 0.0
    t_solve = 0.0
    i_extract = 0
    σ = -1
    for (sr,ls) in enumerate(@view layers[[begin,end]]) # Loop over source regions
        σ *= -1 # 1 for sr == 1, and -1 for sr == 2
        for qp in 1:length(ls.P) # Loop over srce reg modes
            # Incident field for Region sr (Eq. (7.64))
            sourcevec = iincs[qp,sr] * ls.c[qp] * ls.Y[qp] * zhatcross(ls.tvec[qp])
            # Compute generalized current vector:
            vmat = [b ⋅ sourcevec for b in bfftstore[:,sr,qp]] # Eq. (7.64)
            # Solve the matrix equation
            t_solve1 = time()
            vmat = ymatf \ vmat
            t_solve2 = time()
            t_solve += t_solve2 - t_solve1
            nsolve += 1
            t_extract1 = time()
            for (or,lo) in enumerate(@view layers[[begin,end]]) # Loop over obs. regions
                smat = gsm[or,sr]
                for q in 1:length(lo.P)  # Loop obs. regn. modes
                    # Extract partial scattering parameter due to scattered fields...
                    FTM = sum((vmat[n] * bfftstore[n,or,q] for n in 1:nbf)) # FT of total mag. current
                    smat[q,qp] += (zhatcross(lo.tvec[q]) ⋅ FTM) * 
                                        (σ * tlgfiv[q,or] * lo.c[q]) # Eq. (6.37)
                    i_extract += 1
                end
            end
            t_extract2 = time()
            t_extract += t_extract2 - t_extract1
        end
    end
    @info "  Extracting $(i_extract) GSM entries used $(t_extract) sec"
    return gsm
end

end
