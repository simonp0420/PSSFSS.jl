module Modes

using Unitful: ustrip, @u_str
using LinearAlgebra: norm, ⋅, ×
using StaticArrays: MVector
using ..Constants: twopi, η₀
using ..Layers: Layer, TEorTM, TE, TM
using ..Sheets: Sheet, find_unique_periods
using ..GSMs: Gblock
using ..Rings: Ring
using ..PGF: mysqrt

const MNmax_default = 6

"""
    choose_layer_modes!(strata, gbl, k0max, dbmin)

    #
Set up the arrays of modal indices `M`, `N`, and `P` for all layers not included
in `Gblock`s given the wavenumber at the maximum operating frequency `k0max`,  
the list of layers and sheets `strata`, the list of `Gblock`s `gbl`, and the 
desired minimum attenuation `dbmin` for each neglected mode to encounter
when passing through the layer.  Also, allocate the arrays `β`, `γ`, `Y`, 
`c`, and `tvec` for each of the excluded layers.  Note that the modes are not 
necessarily defined consistently in each layer because the periodicity of the 
multiple FSS/PSS sheets may not all be identical.  Therefore, we store the β₁ 
and β₂ reciprocal lattice vectors for each layer in the `Layer` type and use 
these values to define the periodicity for the modes in a given layer.

## Arguments

- `strata`: A vector whose elements are either `Layer` or `Sheet` objects. `strata` defines the planar
geometry being analyzed.  It is assumed that for each `Layer` in `strata` not included in a `Gblock`, 
the permeability and permittivity have been correctly initialized.  On exit, these same excluded layers 
will have fields β₁ and β₂ appropriately set, and will have the arrays `β`, `γ`, `Y`, 
`c`, and `tvec` allocated for the  appropriate number of modes.

- `gbl`: An collection of `Gblock`s containing the definitions of the GSM block entities.

- `k0max`: Maximum free-space wavenumber to be analyzed (1/m).

- `dbmin`  Minimum attenuation any neglected mode must incur (dB > 0).
"""
function choose_layer_modes!(strata, gbl, k0max, dbmin)
    T = Union{Layer,Sheet}
    all(t isa T for t in strata) || error("strata elements must be of type Sheet or Layer")
    islayer = map(t -> t isa Layer, strata)
    layers = @view strata[islayer]
    nl = length(layers)
    nj = nl-1 # number of dielectric junctions
    issheet = map(x -> x isa Sheet, strata)
    sheets = @view strata[issheet]
    ns = length(sheets)
    ngbl = length(gbl)
    sint = cumsum(islayer)[issheet] # sint[k] contains dielectric interface number of k'th sheet 
    junc = zeros(Int, nj)
    junc[sint] = 1:ns #  junc[i] is the sheet number present at interface i, or 0 if no sheet is there


    # Make all layers unique:
    for i in 2:length(layers)
        for j in 1:i-1
            if layers[i] === layers[j]
                layers[i] = deepcopy(layers[j])
                @goto nexti
            end
        end
        @label nexti
    end

    # Possible mode indices m and n will range from -MNmax to MNmax,for a total of 2(2MNmax+1)^2 
    # modes in the layer (first factor of 2 is for both TE/TM).
    mset = falses(nl)  # Indicates that modes have not been determined yet for any layers
    for g in gbl  # Don't need to set modes for layers included in a Gblock
        mset[g.rng[2:end]] .= true
    end

    pure_radome = true
    for isht in junc
      if isht ≠ 0 && sheets[isht].style ≠ "NULL"
        pure_radome = false 
        break
      end
    end 

    if pure_radome
        for l in layers
            l.P = [TE, TM]
            l.M = [0, 0]
            l.N = [0, 0]
            l.β₁ = MVector(twopi, 0.0)
            l.β₂ = MVector(0.0, twopi) # cell 1m in x by 1m in y.
        end
        mset .= true

    else   # There is at least one non-NULL sheet in the structure
        last_gbl = findlast(x -> !iszero(x.j) && sheets[junc[x.j]].style ≠ "NULL", gbl)
        last_sheet = junc[gbl[last_gbl].j] # sheet index
        first_gbl = findfirst(x -> !iszero(x.j) && sheets[junc[x.j]].style ≠ "NULL", gbl)
        first_sheet = junc[gbl[first_gbl].j] # sheet index 
        # Assign equivalence class numbers to the interfaces based on the 
        # unit cell of the sheet located there:
        upa = find_unique_periods(junc, sheets) 
        # Treat all unassigned layers that are adjacent to a gbl containing an FSS screen:
        for igbl in 1:ngbl
            g = gbl[igbl]
            ijunc = g.j # Junction # where FSS is located
            ijunc == 0 && continue
            isht = junc[ijunc]  # Sheet number at junction ijunc
            sheet = sheets[isht]
            sheet.style == "NULL" && continue
            one_meter = ustrip(Float64, sheet.units, 1u"m")
            β₁ = sheet.β₁ * one_meter # radians/meter
            β₂ = sheet.β₂ * one_meter # radians/meter
            i1 = g.rng[1] # Index of layer to left of current GBLOCK.
            i2 = 1 + g.rng[end] # Index of layer to right of current GBLOCK.
            for i in [i1, i2]
                (i==1 || i==nl) && continue  # Skip boundary layers.
                # Determine neighboring gblock to the outside: 
                if i == i1
                    igbl_other = max(igbl-1, 1)
                else
                    igbl_other = min(igbl+1, ngbl)
                end
                ijunc_other = gbl[igbl_other].j
                if ijunc_other ≠ 0
                    isht_other = junc[ijunc_other]  # Sheet number at junction ijunc_other
                else
                    isht_other = 0
                end
                if !mset[i]  # This layer not yet assigned
                    layers[i].β₁ = copy(β₁)  # Save periodicity in layer i.
                    layers[i].β₂ = copy(β₂)  # Save periodicity in layer i.
                    MNmax = MNmax_default
                    # These values must be changed to 0 if there is also another GBLOCK
                    # adjacent to layer i containing a nonnull sheet with a different periodicity:
                    if ijunc_other ≠ 0 && upa[ijunc] ≠ upa[ijunc_other] &&
                       sheets[isht_other].style ≠ "NULL" 
                        error("""
                               Unequal unit cells in sheets $(isht) and $(isht_other)
                                      Try dividing layer $(i) into 2 halves.""")

                    end
                    # Select the mode indices for layer i:
                    fill_mmax_pmn!(layers[i], MNmax, k0max, dbmin)
                    mset[i] = true
                else
                    # mset is true, indicating that this layer belongs to a gblock containing an FSS.
                    # IF the periodicity of that FSS is inconsistent with that of the current gblock
                    # then we can't satisfy both requirements, so abort.
                    if ijunc_other ≠ 0 && upa[ijunc] ≠ upa[ijunc_other] && sheets[isht_other].style ≠ "NULL"
                        error("""
                        Unequal unit cells in sheets $(isht) and $(isht_other).
                               Try dividing layer $(i) into 2 halves. """)
                    end
                end
            end
        end # for
    

        # At this point, all layers adjacent to a Gblock containing an FSS sheet have 
        # been assigned modes.  Now we assign mode index values to the remaining layers 
        # by stepping inwards from the already assigned layers. Loop until all 
        # layers are assigned.
        while any(.!mset[2:nl-1])
            for i in 2:nl-1  # Step over each layer
                mset[i] && continue
                if mset[i-1] && mset[i+1]
                    # Both adjacent layers have been assigned.
                    if equal_layer_periodicity(layers[i-1],layers[i+1])
                        # The two adjacent layers have the same periodicity.  Choose
                        # the accessible modes from among the two sets of adjacent ones.
                        layers[i].β₁ = copy(layers[i+1].β₁)
                        layers[i].β₂ = copy(layers[i+1].β₂)
                        layers[i].P = vcat(layers[i-1].P, layers[i+1].P)
                        layers[i].M = vcat(layers[i-1].M, layers[i+1].M)
                        layers[i].N = vcat(layers[i-1].N, layers[i+1].N)
                        select_pmn!(layers[i], k0max, dbmin)
                        mset[i] = true
                    else
                        # Both adjacent layers have been assigned, but they have different
                        # periodicities.  If they each have only two modes (TE/TM), then
                        # all is well, since these m=n=0 modes are independent of the
                        # periodicity.  If not, we must truncate the modes to only these
                        # two, and print a warning message:
                        layers[i].β₁ = copy(layers[i+1].β₁) # Assign same periodicity.
                        layers[i].β₂ = copy(layers[i+1].β₂)
                        layers[i].P = [TE, TM]
                        layers[i].M = [0, 0]
                        layers[i].N = [0, 0]
                        mset[i] = true
                        if length(layers[i-1].P) ≠ 2 || length(layers[i+1].P) ≠ 2
                            @warn  """"
                            Setting # modes in Layer $(i) to 2 due to 
                                unequal unit cells in surrounding FSS sheets """
                        end
                    end
                elseif  mset[i-1] && !mset[i+1]
                    # Only layer to the left has been assigned.
                    i-1 == 1 && @goto nextouter # Don't rely on end layer.
                    copy_βs_pmn!(layers[i], layers[i-1])
                    select_pmn!(layers[i], k0max, dbmin)
                    mset[i] = true
                elseif !mset[i-1] && mset[i+1]
                    # Only layer to the right has been assigned.
                    i+1 == nl && @goto nextouter # Don't rely on end layer.
                    copy_βs_pmn!(layers[i], layers[i+1])
                    select_pmn!(layers[i], k0max, dbmin)
                    mset[i] = true
                end # if
            end # for
            @label nextouter
        end # while
        # Done with all interior layers.
        # Choose numbers of modes in the end layers. Any mode that can 
        # propagate at k0max will be defined. Use the periodicity of the 
        # closest FSS sheets for each end layer.
        for (l, sh) in ((layers[1], sheets[first_sheet]), (layers[nl], sheets[last_sheet]))
            one_meter = ustrip(sh.units, 1u"m")
            l.β₁ = sh.β₁ * one_meter
            l.β₂ = sh.β₂ * one_meter
            fill_mmax_pmn!(l, MNmax_default, k0max, dbmin)
            select_pmn!(l, k0max) # eliminate evanescant modes
        end
        mset[1] = mset[nl] = true

        # Ensure that any adjacent non-gblock layers either have the same 
        # periodicity or have only dominant modes:
        inagblock = reduce(union, (g.rng[2:end] for g in gbl)) # layer indices in Gblocks
        notinablock = sort(setdiff(1:nl, inagblock))
        for i in 1:length(notinablock)-1
            notinablock[i+1] - notinablock[i] == 1 || continue
            il = notinablock[i] # layers[il] and layers[il+1] are adjacent
            equal_layer_periodicity(layers[il], layers[il+1]) && continue
            for l in layers[il:il+1]
                l.P = [TE,TM]
                l.M = [0, 0]
                l.N = [0, 0]
            end
            @warn """
            Setting # modes in layers $(il) and $(il+1) to 2 due to
                unequal unit cells in surrounding FSS sheets.
            """
        end
    end # if


    # Now we allocate the other Layer field vectors
    maxmodes = 2*(2*MNmax_default+1)^2
    for i in notinablock
        layer = layers[i]
        nmodes = length(layer.P)
        for field in (:β, :γ, :Y, :c, :tvec)
            setfield!(layer, field, zeros(eltype(getfield(layer,field)), nmodes))
        end
        nmodes == maxmodes && @warn """
        Maximum number of modes reached. GSM results may not be accurate.
                 Consider increasing MNmax_default from its current value $(MNmax_default)"""
    end

    return nothing
end

"""
    setup_modes!(layer::Layer, k0::Real, kvec::AbstractVector)

Fill the modal layer fields.  Needed for layers not contained in a Gblock.
The arrays are assumed to have  been already allocated, and the index arrays 
`M`, `N`, and `P` are asssumed to have been already initialized.

### Input arguments

- `layer`: It is assumed that the `ϵᵣ`, `μᵣ`, `β₁`, and `β₂` have been 
    initialized, as have the modal index arrays `P`, `M`, and `N`.
- `k0`: Free-space wavenumber (1/meter).
- `kvec`: A real-valued 2-vector containing the x and y components of the incident
    plane wave unit vector that defines the unit cell incremental 
    phase shifts.

### Outputs

There is no explicit output, but the fields of `layer` will be modified, including
`β`, `tvec`, `γ`, `c`, and `Y`. 
"""
function setup_modes!(layer::Layer, k0::Real, kvec::AbstractVector)
    β₀₀ = MVector(kvec[1], kvec[2])
    β₁, β₂ = layer.β₁, layer.β₂
    area = twopi^2 / norm(β₁ × β₂)
    ksq = k0^2 * layer.ϵᵣ * layer.μᵣ
    for mode in 1:length(layer.M)
        m,n,p = layer.M[mode], layer.N[mode], layer.P[mode]
        β = β₀₀ + m*β₁ + n*β₂
        β² = β ⋅ β
        β̂ = β²*area < 1e-14 ? MVector(1.0, 0.0) : β/norm(β)
        layer.β[mode] .= β
        layer.γ[mode] = γ = mysqrt(β² - ksq)
        if p == TE
            Y = γ / (im * (k0*η₀) * layer.μᵣ)
            tvec = zhatcross(β̂)
        else
          Y = im * (k0/η₀) * layer.ϵᵣ / γ
          tvec = β̂
        end
        layer.Y[mode] = Y
        layer.tvec[mode] = tvec
        layer.c[mode] = mysqrt(1 / (area * Y))
    end
    return nothing
end

function zhatcross(x)
    y = similar(x)
    y[2] = x[1]
    y[1] = -x[2]
    y
end


"""
    equal_layer_periodicity(l1::Layer, l2::Layer) -> Bool

Determine if `l1` and `l2` have identical unit cells
by examining their reciprocal lattice vectors.
"""
function equal_layer_periodicity(l1::Layer, l2::Layer)
    denom = norm(l1.β₁) + norm(l1.β₂) + norm(l2.β₁) + norm(l2.β₂)
    norm(l1.β₁ - l2.β₁) / denom < 1e-5 && (norm(l1.β₂ - l2.β₂) / denom < 1e-5)
end



"""
copy_βs_pmn!(l1::Layer, l2::Layer)

Copy fields `β₁`, `β₂`, `P`, `M`, `N` from `l2` to `l1`, modifying `l1`.
"""
@inline function copy_βs_pmn!(l1::Layer, l2::Layer)
    l1.β₁ = copy(l2.β₁)
    l1.β₂ = copy(l2.β₂)
    l1.P = copy(l2.P)
    l1.M = copy(l2.M)
    l1.N = copy(l2.N)
    nothing
end



"""
    fill_mmax_pmn!(layer::Layer, MNmax::Int, k0max, dbmin=-Inf)

Fill the `layer` modal index arrays `M`, `N`, and `P` given the `layer`
material parameters and reciprocal lattice vectors, the wavenumber `k0max` 
at the maximum analysis frequency, and the minimum desired attenuation in dB 
`dbmin` for the neglected modes. `MNmax` is the maximum ring number in the 
mode lattice to consider.
"""
function fill_mmax_pmn!(layer::Layer, MNmax::Int, k0max, dbmin)
    MNmax ≥ 0 || error("MNmax must be nonnegative")
    n = 2 * (2MNmax+1)^2
    layer.M = zeros(Int, n)
    layer.N = zeros(Int, n)
    layer.P = [TE for _ in 1:n]
    mode = 1
    for ring in 0:MNmax, (m,n) in Ring(ring)
        layer.M[mode:mode+1] .= m
        layer.N[mode:mode+1] .= n
        layer.P[mode:mode+1] .= (TE,TM)
        mode += 2
    end
    select_pmn!(layer, k0max, dbmin)
    return nothing
end



"""
    select_pmn!(layer:Layer, k0max, dbmin)

Select the accessible modes from among the candidates provided
in the `layer`. A mode is considered not accessible if it undergoes more
than `dbmin` decibels of attenuation when propagating/attenuating through 
the thickness of dielectric substrate `layer` at the highest frequency of 
operation, for any possible angle of incidence of the plane-wave that excites
the FSS structure.

The special value `dbmin=-Inf` instructs this function to retain only
propagating modes
"""
function select_pmn!(layer::Layer, k0max, dbmin=-Inf)
    β₁, β₂ = layer.β₁, layer.β₂  # Reciprocal lattice vectors (1/m)
    kmaxsq = k0max^2  * (layer.ϵᵣ * layer.μᵣ)
    if isfinite(dbmin)
        αh_test = dbmin / 8.6859 # Convert from dB to nepers
        keep = map(1:length(layer.P)) do mode
            (layer.M[mode] == layer.N[mode] == 0) && (return true)
            # Obtain min value of ||βₘₙ|| for this mode for any angle of incidence
            β = norm(layer.M[mode]*β₁ + layer.N[mode]*β₂)
            βmin = β ≤ k0max ? 0.0 : β - k0max
            γ = mysqrt(βmin^2 - kmaxsq) # modal complex atten. constant
            αh = real(γ) * layer.width
            return  αh < αh_test
        end
    else
        keep = map(1:length(layer.P)) do mode
            βvec = layer.M[mode]*β₁ + layer.N[mode]*β₂
            β² = βvec ⋅ βvec
            γ = mysqrt(β² - kmaxsq) # modal complex atten. constant
            return  imag(γ) ≥ real(γ)
        end
    end
    layer.M = layer.M[keep]
    layer.N = layer.N[keep]
    layer.P = layer.P[keep]

    # Eliminate possible duplicate mode triples (p,m,n)
    MN = repeat(unique([layer.M  layer.N], dims=1), inner=(2,1))
    layer.M = MN[:,1]
    layer.N = MN[:,2]
    layer.P = repeat([TE,TM], outer=length(layer.M)÷2)
    return nothing
end

end # module
