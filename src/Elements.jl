module Elements

export diagstrip, jerusalemcross, loadedcross, manji, meander, pecsheet, pmcsheet, polyring, rectstrip, splitring

using ..PSSFSSLen: mm, cm, inch, mil, PSSFSSLength
using ..Sheets: RWGSheet, rotate!, translate!, combine, recttri, SV2
using ..Meshsub: meshsub
using StaticArrays: SA
using LinearAlgebra: norm, ⋅
using Printf: @sprintf
import LibGEOS # difference, readgeom, Polygon, MultiPolygon
import GeoInterface # nhole, ngeom, coordinates, getexterior, gethole
import PolygonOps

macro testpos(var)
    return :(all($(esc(var)) .> 0) || error($(esc(string(var))) * " must be positive!"))
end

mutable struct MeshsubData
    ρ::Vector{SV2}
    e1::Vector{Cint}
    e2::Vector{Cint}
    segmarkers::Vector{Cint}
    holes::Vector{SV2}
    boundary::Int
    node::Int
    area::Float64
end
MeshsubData() = MeshsubData(SV2[], Cint[], Cint[], Cint[], SV2[], 0, 0, 0.0)

"""
    _add_libgeos_geom!(msdata::MeshsubData, obj::LibGEOS.Polygon, ρ₀)

Update `msdata` for the new LibGEOS Polygon object. If the object has any holes,
then it is assumed that the polygon is actually a regular polygonal annular ring
centered on the point ρ₀.
"""
function _add_libgeos_geom!(msdata::MeshsubData, obj::LibGEOS.Polygon, ρ₀)
    allcoords = GeoInterface.coordinates(obj)
    ngeom = GeoInterface.ngeom(obj)
    @assert ngeom == length(allcoords)
    for kgeom in 1:ngeom # 1 for exterior, 2,3,... for holes
        coords = allcoords[kgeom]
        sgn = kgeom == 1 ? 1 : -1
        area = PolygonOps.area(coords)
        if sgn * area < 0
            reverse!(coords) # exterior CCW. holes CW
            area *= -1 # exterior area > 0, hole area < 0
        end
        msdata.area += area
        msdata.boundary += 1
        nodesave = msdata.node + 1
        for (i, ρ) in enumerate(coords)
            i == length(coords) && break # Last point is repeat of first
            msdata.node += 1
            push!(msdata.e1, msdata.node)
            push!(msdata.e2, msdata.node + 1)
            push!(msdata.segmarkers, msdata.boundary)
            push!(msdata.ρ, ρ)
        end
        msdata.e2[end] = nodesave
        if kgeom > 1
            # Add hole point.  This code assumes that the hole is a complete ring centered on ρ₀
            ρᵢₙ = ρ₀ + 0.999 * (msdata.ρ[end] - ρ₀)
            @assert PolygonOps.inpolygon(ρᵢₙ, coords, in=true, on=false, out=false)
            push!(msdata.holes, ρᵢₙ)
        end
    end
    return msdata
end

"""
    _add_libgeos_geom!(msdata::MeshsubData, obj::LibGEOS.MultiPolygon, ρ₀)

Update `msdata` for the new LibGEOS MultiPolygon object. If the object has any holes,
then it is assumed that the polygon is actually a regular polygonal annular ring
centered on the point ρ₀.
"""
function _add_libgeos_geom!(msdata::MeshsubData, obj::LibGEOS.MultiPolygon, ρ₀)
    for k in 1:GeoInterface.ngeom(obj)
        geom = LibGEOS.getgeom(obj, k)
        _add_libgeos_geom!(msdata, geom, ρ₀)
    end
    return msdata
end


@inline zhatcross(t::SV2) = SV2(-t[2], t[1])
@inline zhatcross(t::T) where {T<:AbstractVector} = [-t[2], t[1]]

"""
    s₁s₂2β₁β₂(s₁,s₂) -> (β₁, β₂)

Compute the reciprocal lattice vectors from the direct lattice vectors.
Inputs and outputs are static 2-vectors from StaticArrays.
"""
function s₁s₂2β₁β₂(s₁, s₂)
    fact = 2π / abs(s₁[1] * s₂[2] - s₁[2] * s₂[1])
    β₁ = -fact * zhatcross(s₂)
    β₂ = fact * zhatcross(s₁)
    return β₁, β₂
end

function replace_kw_arg!(kwargs, badkw, goodkw)
    if haskey(kwargs, badkw)
        kwargs[goodkw] = kwargs[badkw]
        delete!(kwargs, badkw)
    end
    kwargs
end


"""
    check_optional_kw_arguments!(kwargs:: AbstractDict{Symbol,T} where T)

Check the validity of the optional keyword arguments passed to one of the 
user-callable, specific sheet constructor functions.  If any of the arguments
were not passed, assign appropriate default values. 

Also, replace obsolete `:Rsheet` with `:Zsheet`
"""
function check_optional_kw_arguments!(kwargs::AbstractDict{Symbol,T} where {T})
    for (badkw, goodkw) in zip((:Rsheet, :sigma), (:Zsheet, :σ))
        replace_kw_arg!(kwargs, badkw, goodkw)
    end

    haskey(kwargs, :Zsheet) && haskey(kwargs, :σ) && error("Zsheet and σ cannot both be specified")
    haskey(kwargs, :Zsheet) && haskey(kwargs, :Rq) && error("Zsheet and Rq cannot both be specified")
    haskey(kwargs, :Zsheet) && haskey(kwargs, :disttype) && error("Zsheet and disttype cannot both be specified")


    defaults = Dict(:class => 'J', :dx => 0.0, :dy => 0.0, :rot => 0.0, 
    :Zsheet => 0.0, :σ => -Inf, :Rq => 0.0, :disttype => :normal, :save => "", 
    :fufp => false, :structuredtri => true)
    validkws = keys(defaults)

    badkws = setdiff(keys(kwargs), validkws)
    if !isempty(badkws)
        error("Illegal keywords: ", join(badkws, ", "))
    end
    for (key, val) in defaults
        haskey(kwargs, key) || (kwargs[key] = val)
    end

    class = kwargs[:class]
    class ≠ 'J' && class ≠ 'M' && error("Illegal value $class for class")

    for key in [:dx, :dy, :rot]
        kwargs[key] isa Real || error("$key must be a Real")
    end

    Zsheet = kwargs[:Zsheet]
    class ≠ 'J' && !iszero(Zsheet) && error("Nonzero surface impedance only allowed for J-class sheet")
    real(Zsheet) < 0 && error("real(Zsheet) must be nonnegative")

    σ = kwargs[:σ]
    if class ≠ 'J' 
        σ ≠ -Inf && error("Conductivity only allowed for J-class sheet")
    else
        -Inf < σ ≤ 0 && error("Conductivity must be nonnegative")
    end

    Rq = kwargs[:Rq]
    Rq < 0 && error("Rq must be nonnegative")

    kwargs[:save] isa AbstractString || error("save value must be an AbstractString")

    disttype = kwargs[:disttype]
    disttype == :normal || disttype == :rayleigh || error("Illegal value $disttype for disttype")

    return
end

const optional_kwargs = """
                        - `class::Char='J'`  Specify the class, either `'J'` or `'M'`.. If `'J'`,  the unknowns are electric surface 
                                   currents, as used to model a wire or metallic patch-type FSS.  If `'M'`,  the unknowns are
                                   magnetic surface currents, as used to model a slot or aperture in a perfectly conducting plane.
                        - `dx::Real=0.0`, `dy::Real=0.0`:  These specify the offsets in the x and y directions applied to the entire 
                                   unit cell and its contents.  Length units are as specified in the `units` keyword. 
                        - `rot::Real=0.0`:  Counterclockwise rotation angle in degrees applied to the entire unit cell and its contents. 
                                   This rotation is applied prior to any offsets specified in `dx` and `dy`.
                        - `Zsheet::Complex=0.0`:  The frequency-independent surface impedance of the FSS conductor in 
                          units of [Ω].  May only be specified for a sheet of class `'J'`.  If `Zsheet` is specified, then 
                          `sigma` (or `σ`) may not be specified.                          )
                        - `sigma` or `σ`: DC, bulk conductivity [S/m].  Only allowed for sheets of class `'J'`.  Cannot be 
                          simultaneously specified with `Zsheet`.  Is used with `Rq` by PSSFSS to calculate an effective 
                          sheet surface impedance at each frequency, using the Gradient Model (Grujić 2022).
                        - `Rq=0.0`: RMS surface roughness [m].  Only legal for class `'J'`. Only used if `sigma` (or `σ`) is 
                           also specified.  In that case is is used along with `sigma` to calculate a frequency-dependent
                           sheet impedance using the Gradient Model.  The default value of 0 denotes a smooth surface.
                        - `disttype::Symbol=:normal`: Probability distrubution type for surface roughness.  defaults
                          to `:normal`.  The other legal value is `:rayleigh`.
                        - `fufp::Bool`:  This keyword is not usually required. 
                          `fufp` is mnemonic for "Find Unique Face Pairs".  If true, the code will search the 
                          triangulation for classes of triangle
                          pairs that are the equivalent in the toeplitz sense.  I.e., if triangle pairs (A,B) and (C,D) belong
                          to the same equivalence class,  the six vertices in the pair (A,B) can be made to coincide 
                          with those of pair (C,D) by a simple translation. If there are many such equivalent pairs, 
                          a significant decrease in matrix fill time ensues by exploiting the equivalence.  The tradeoff
                          is the time needed to identify them.  The default value is `true` for the `strip`, `diagstrip`,  
                          `meander`, `manji`, `loadedcross`, `jerusalemcross`, and 4-sided `polyring` styles (those employing 
                          structured meshes) and `false` for the remaining styles (those employing unstructured meshes).
                        - `save::String=""` Specifies a file name to which the sheet triangulation and unit cell data is to be written,
                          typically to be plotted later.
                                
                        """



"""
    diagstrip(;P::Real, w::Real, orient::Real, Nl::Int, Nw::Int, units::PSSFSSLength, kwargs...)

Return a variable of type `RWGSheet` that contains the triangulation for a rectangular strip inside
a square unit cell, with the strip centerline coincident with one of the diagonals of the cell.  The strip
runs the full length of the diagonal.

# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `P`: The period (i.e. side length) of the square unit cell. Note that the center-center spacing of the strips
  is `P/√2`.
- `w`: The width of the strip.
- `orient`: The orientation of the strip within the unrotated unit cell in degrees.  The only valid values
  are `45` for a strip running from lower left to upper right and `-45` for a strip running from lower 
  right to upper left.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `Nl` and `Nw`:  Number of line segments along the length and width of the strip, for dividing up the strip into
  rectangles, which are  triangulated by adding a diagonal to each rectangle. These arguments are actually used for 
  triangulating the central, rectangular portion of the strip.  The ends of the strip are tapered in the form of 
  right, isosceles triangles, to conform to the boundaries of the square unit cell.  These triangular "end-caps" 
  are triangulated using an unstructured mesh.
    
## Optional arguments:
$(optional_kwargs)
"""
function diagstrip(; P::Real, w::Real, orient::Real, Nl::Int, Nw::Int, units::PSSFSSLength, kwarg...)::RWGSheet
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(P)
    @testpos(w)
    @testpos(Nl)
    @testpos(Nw)
    dxdy = [kwargs[:dx], kwargs[:dy]]
    if dxdy ≠ [0, 0]
        error("Translation not allowed for this style of sheet")
    end

    abs(orient) == 45 || error("orient must be 45 or -45")

    # Setup structured triangulation for central part:
    Lx = √2 * P - w # Length of rectangular portion
    xbl = w / 2
    orient < 0 && (xbl -= √2 * P)
    rhobl = SV2([xbl, -w / 2])
    rhotr = rhobl + SV2([√2 * P - w, w])
    sh1 = recttri(rhobl, rhotr, Nl, Nw)

    # Right triangular region:
    sh2 = translate!(rotate!(_tritri(w, Nw), 180.0), rhotr[1] + w / 2, 0)
    sh3 = combine(sh1, sh2, 'x', rhotr[1])
    # Left triangular region:
    sh2 = translate!(_tritri(w, Nw), xbl - w / 2, 0)
    sh1 = combine(sh3, sh2, 'x', xbl)
    # Final rotation:
    rotate!(sh1, orient)
    orient < 0 && translate!(sh1, P, 0)
    # Corner pieces:
    if orient > 0
        # top left patch:
        sh2 = translate!(rotate!(_tritri(w, Nw), -45.0), 0, P)
        sh3 = combine(sh1, sh2, ' ', 0.0)
        # bottom right patch:
        sh1 = translate!(rotate!(_tritri(w, Nw), 135.0), P, 0)
        sheet = combine(sh3, sh1, ' ', 0.0)
    else
        # top right patch:
        sh2 = translate!(rotate!(_tritri(w, Nw), -135.0), P, P)
        sh3 = combine(sh1, sh2, ' ', 0.0)
        # bottom left patch:
        sh1 = translate!(rotate!(_tritri(w, Nw), 45.0), 0, 0)
        sheet = combine(sh3, sh1, ' ', 0.0)
    end


    sheet.style = "diagstrip"
    sheet.units = units

    sheet.s₁ = SV2([P, 0.0])
    sheet.s₂ = SV2([0.0, P])
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])

    sheet.ξη_check = true

    return sheet

end # function

"""
    jerusalemcross(;P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that contains the triangulation for a 
"Jerusalem cross" type of geometry.
The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, 
and `fr` properly initialized.


The following "ascii art" attempts to show
the definitions of the geometrical parameters `P`, `L1`, `L2`, `A`, `B`, and `w`.
Note that the structure is supposed to be symmetrical wrt reflections
about its horizontal and vertical centerlines, and wrt reflections through a line oriented
at a 45 degree angle wrt the x-axis.


    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓ 
    ┃                                                       ┃ _______
    ┃               ┌────────────────────────┐              ┃    ↑
    ┃               │ ┌───────────────────┐  │              ┃    │
    ┃               │ └───────┐    ┌──────┘  │              ┃    │
    ┃               └──────┐  │    │ ┌───────┘              ┃    │
    ┃                      │  │    │ │                      ┃    │
    ┃  ┌───────┐           │  │    │ │            ┌──────┐  ┃    │
    ┃  │  ┌─┐  │           │  │    │ │            │ ┌──┐ │  ┃    │
    ┃  │  │ │  │           │  │   →│ │← w         │ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  │ │  └───────────┘  │    │ └────────────┘ │  │ │  ┃    │
    ┃  │  │ └─────────────────┘    └────────────────┘  │ │  ┃    
    ┃  │  │                                            │ │  ┃   L1 
    ┃  │  │ ┌─────────────────┐    ┌────────────────┐  │ │  ┃  
    ┃  │  │ │  ┌───────────┐  │    │ ┌────────────┐ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  └─┘  │          →│  │    │ │← L2     B →│ └──┘ │← ┃    │
    ┃  └───────┘           │  │    │ │            └──────┘  ┃    │
    ┃                      │  │    │ │                      ┃    │
    ┃               ┌──────┘  │    │ └───────┐              ┃    │
    ┃               │ ┌───────┘    └──────┐  │              ┃    │
    ┃               │ └───────────────────┘  │              ┃    │
    ┃               └────────────────────────┘              ┃ ___↓___
    ┃               |<───────── A ──────────>|              ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛ 
    |<─────────────────────── P ───────────────────────────>|
                        
    
    
# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `P`: The period, i.e. the side length of the square unit cell.
- `L1`,`L2`, `A`, `B`, `w`: Geometrical parameters as defined above.  Note that it is permissible
   to specify `w ≥ L2/2` and/or `w ≥ B/2` in which case the respective region will
   be filled in solidly with triangles.  If both conditions hold, then the entire structure will be
   filled in (i.e., singly-connected).  In that case the `L2` and `B` dimensions will be used 
   for the respective widths of the arms, and `w` will not be used.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
  the actual number will likely be different.
    
## Optional arguments:
$(optional_kwargs)
- `structuredtri::Bool=true`: If true, use a structured mesh for the triangulation.  If false,
  the unstructured mesh generator that was standard up to PSSFSS version 1.2 will be used. A structured 
  mesh can be analyzed more efficiently, but the number of triangles created by the unstructured
  mesh generator is usually closer to `ntri` than the number for the structured mesh generator.
"""
function jerusalemcross(; P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real,
    ntri::Int, units::PSSFSSLength, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    structuredtri = haskey(kwargs, :structuredtri) ? kwargs[:structuredtri] : true
    fufp = haskey(kwargs, :fufp) ? kwargs[:fufp] : structuredtri

    if structuredtri
        return jerusalemcross_structured(; fufp, P, L1, L2, A, B, w, ntri, units, kwargs...)
    else
        return jerusalemcross_unstructured(; fufp, P, L1, L2, A, B, w, ntri, units, kwargs...)
    end

end


"""
    jerusalemcross_unstructured(;P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that contains the triangulation for a 
"Jerusalem cross" type of geometry, using an unstructured mesh.
The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, 
and `fr` properly initialized.


The following "ascii art" attempts to show
the definitions of the geometrical parameters `P`, `L1`, `L2`, `A`, `B`, and `w`.
Note that the structure is supposed to be symmetrical wrt reflections
about its horizontal and vertical centerlines, and wrt reflections through a line oriented
at a 45 degree angle wrt the x-axis.


    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓ 
    ┃                                                       ┃ _______
    ┃               ┌────────────────────────┐              ┃    ↑
    ┃               │ ┌───────────────────┐  │              ┃    │
    ┃               │ └───────┐    ┌──────┘  │              ┃    │
    ┃               └──────┐  │    │ ┌───────┘              ┃    │
    ┃                      │  │    │ │                      ┃    │
    ┃  ┌───────┐           │  │    │ │            ┌──────┐  ┃    │
    ┃  │  ┌─┐  │           │  │    │ │            │ ┌──┐ │  ┃    │
    ┃  │  │ │  │           │  │   →│ │← w         │ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  │ │  └───────────┘  │    │ └────────────┘ │  │ │  ┃    │
    ┃  │  │ └─────────────────┘    └────────────────┘  │ │  ┃    
    ┃  │  │                                            │ │  ┃   L1 
    ┃  │  │ ┌─────────────────┐    ┌────────────────┐  │ │  ┃  
    ┃  │  │ │  ┌───────────┐  │    │ ┌────────────┐ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  └─┘  │          →│  │    │ │← L2     B →│ └──┘ │← ┃    │
    ┃  └───────┘           │  │    │ │            └──────┘  ┃    │
    ┃                      │  │    │ │                      ┃    │
    ┃               ┌──────┘  │    │ └───────┐              ┃    │
    ┃               │ ┌───────┘    └──────┐  │              ┃    │
    ┃               │ └───────────────────┘  │              ┃    │
    ┃               └────────────────────────┘              ┃ ___↓___
    ┃               |<───────── A ──────────>|              ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛ 
    |<─────────────────────── P ───────────────────────────>|
                        
    
    
# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `P`: The period, i.e. the side length of the square unit cell.
- `L1`,`L2`, `A`, `B`, `w`: Geometrical parameters as defined above.  Note that it is permissible
   to specify `w ≥ L2/2` and/or `w ≥ B/2` in which case the respective region will
   be filled in solidly with triangles.  If both conditions hold, then the entire structure will be
   filled in (i.e., singly-connected).  In that case the `L2` and `B` dimensions will be used 
   for the respective widths of the arms, and `w` will not be used.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
  the actual number will likely be different.
    
## Optional arguments:
$(optional_kwargs)
"""
function jerusalemcross_unstructured(; P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real,
    ntri::Int, units::PSSFSSLength, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    @testpos(A)
    @testpos(B)
    @testpos(L1)
    @testpos(L2)
    @testpos(w)
    @testpos(ntri)
    @testpos(P)


    areaouter = 4 * (A * B + ((L1 - L2) / 2 - B) * L2) + L2^2   # outer area for solid cross

    # Total number of vertices and holes and total area:
    if 2w < L2 && 2w < B
        nv = 28 + 28
        nholes = 1
        areat = areaouter - 4 * ((A - 2w) * (B - 2w) + (L2 - 2w) * ((L1 - L2) / 2 - B)) - (L2 - 2w)^2
    elseif 2w ≥ L2 && 2w ≥ B
        nv = 28
        nholes = 0
        areat = areaouter
    elseif 2w < L2 && 2w ≥ B
        nv = 28 + 12
        nholes = 1
        areat = areaouter - 4 * (L2 - 2w) * ((L1 - L2) / 2 - B) - (L2 - 2w)^2
    elseif 2w ≥ L2 && 2w < B
        nv = 28 + 4 * 4
        nholes = 4
        areat = areaouter - 4 * (A - 2w) * (B - 2w) - (L2 - 2w)^2
    end

    s1 = Cdouble[P, 0.0]
    s2 = Cdouble[0.0, P]
    r0 = 0.5 * (s1 + s2) # Calculate center of polygon.
    r = zeros(Cdouble, 2, nv)
    e1 = Array{Cint}(undef, nv)
    e2 = Array{Cint}(undef, nv)
    segmarkers = Array{Cint}(undef, nv)
    holes = Array{Cdouble}(undef, 2, nholes)

    # Set up the (outer) polygon geometry:
    r[:, 1] = [L1 / 2, A / 2]
    r[:, 2] = [L1 / 2 - B, A / 2]
    r[:, 3] = [L1 / 2 - B, L2 / 2]
    r[:, 4] = [L2 / 2, L2 / 2]
    r[:, 5] = [L2 / 2, L1 / 2 - B]
    r[:, 6] = [A / 2, L1 / 2 - B]
    r[:, 7] = [A / 2, L1 / 2]
    r[:, 8:14] = reverse([-1 0; 0 1] * r[:, 1:7], dims=2)
    r[:, 15:28] = reverse([1 0; 0 -1] * r[:, 1:14], dims=2)
    e1[1:28] = 1:28
    e2[1:27] = 2:28
    e2[28] = 1
    segmarkers[1:28] .= 1

    if 2w < L2 && 2w < B
        #  Set up inner boundary for annulus:
        r[:, 29] = r[:, 1] + [-w, -w]
        r[:, 30] = r[:, 2] + [w, -w]
        r[:, 31] = r[:, 3] + [w, -w]
        r[:, 32] = r[:, 4] + [-w, -w]
        r[:, 33] = r[:, 5] + [-w, w]
        r[:, 34] = r[:, 6] + [-w, w]
        r[:, 35] = r[:, 7] + [-w, -w]
        r[:, 36:42] = reverse([-1 0; 0 1] * r[:, 29:35], dims=2)
        r[:, 43:56] = reverse([1 0; 0 -1] * r[:, 29:42], dims=2)
        e1[29:56] = 29:56
        e2[29:55] = 30:56
        e2[56] = 29
        segmarkers[29:56] .= 2
        holes[:, 1] = [0, 0]
    elseif 2w < L2 && 2w ≥ B
        r[:, 29] = [L1 / 2 - B, L2 / 2 - w]
        r[:, 30] = [L2 / 2 - w, L2 / 2 - w]
        r[:, 31] = [L2 / 2 - w, L1 / 2 - B]
        r[:, 32:34] = reverse([-1 0; 0 1] * r[:, 29:31], dims=2)
        r[:, 35:40] = reverse([1 0; 0 -1] * r[:, 29:34], dims=2)
        e1[29:40] = 29:40
        e2[29:39] = 30:40
        e2[40] = 29
        segmarkers[29:40] .= 2
        holes[:, 1] = [0, 0]
    elseif 2w ≥ L2 && 2w < B
        r[:, 29] = [L1 / 2 - B + w, -A / 2 + w]
        r[:, 30] = [L1 / 2 - w, -A / 2 + w]
        r[:, 31] = [L1 / 2 - w, A / 2 - w]
        r[:, 32] = [L1 / 2 - B + w, A / 2 - w]
        e1[29:32] = 29:32
        e2[29:31] = 30:32
        e2[32] = 29
        segmarkers[29:32] .= 2
        holes[:, 1] = [(L1 - B) / 2, 0.0]

        r[:, 33:36] = [0 1; 1 0] * r[:, 29:32]
        e1[33:36] = 33:36
        e2[33:35] = 34:36
        e2[36] = 33
        segmarkers[33:36] .= 3
        holes[:, 2] = [0.0, (L1 - B) / 2]

        r[:, 37:40] = [-1 0; 0 1] * r[:, 29:32]
        e1[37:40] = 37:40
        e2[37:39] = 38:40
        e2[40] = 37
        segmarkers[33:36] .= 4
        holes[:, 3] = -holes[:, 1]

        r[:, 41:44] = [0 1; 1 0] * r[:, 37:40]
        e1[41:44] = 41:44
        e2[41:43] = 42:44
        e2[44] = 41
        segmarkers[41:44] .= 5
        holes[:, 4] = -holes[:, 2]

    end



    r .+= r0 # Center on unit cell
    isempty(holes) || (holes .+= r0)

    # Set up call to meshsub
    areatri = areat / ntri
    points = r
    segments = convert(Matrix{Cint}, transpose(hcat(e1, e2)))
    sheet = meshsub(points=points, seglist=segments, segmarkers=segmarkers,
        holes=holes, area=areatri, ntri=ntri)

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "jerusalemcross"
    sheet.ξη_check = false
    sheet.units = units
    sheet.s₁ = s1
    sheet.s₂ = s2
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function


"""
    loadedcross(;s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that
contains the triangulation for a "loaded cross" type of geometry.
The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, 
and `fr` properly initialized.


The following (very poor) "ascii art" attempts to show
the definitions of the geometrical parameters `L1`, `L2` and `w`.
Note that the structure is supposed to be symmetrical wrt reflections
about its horizontal and vertical centerlines, and wrt reflections through a line oriented
at a 45 degree angle wrt the x-axis.


     ^                 ----------------
     |                 |  _________   |
     |                 |  |       |   |
     |                 |  |       |   |
     |                 |  |    -->|   |<--- W
     |                 |  |       |   |
     |                 |  |       |   |
     |     ------------   |       |   -------------
     |     |  |-----------|       |------------|  |
     |     |  |                                |  |
     L1    |  |                                |  |
     |     |  |                                |  |
     |     |  |                                |  |
     |     |  ------------          ------------  |
     |     |-----------   |        |  ------------|
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |________|  |
     |                 |              |
     V                 ----------------
    
                       <---- L2 ------>
    
# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `s1` and `s2`:  2-vectors containing the unit cell lattice vectors.
- `L1`,`L2`,`w`: Geometrical parameters as defined above.  Note that it is permissible
   to specify `w ≥ L2/2` in which case a solid (i.e., singly-connected) cross will be 
   generated.  In that case the `L2` dimension will be used for the width of the cross pieces.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
  the actual number will likely be different.
    

## Optional arguments:
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to locate the initial
           vertex of the loaded cross.  The default is to locate the vertex on the
           positive x-axis.
$(optional_kwargs)
- `structuredtri::Bool=true`: If true, use a structured mesh for the triangulation.  If false,
  the unstructured mesh generator that was standard up to PSSFSS version 1.2 will be used. A structured 
  mesh can be analyzed more efficiently, but the number of triangles created by the unstructured
  mesh generator is usually closer to `ntri` than the number for the structured mesh generator.
"""
function loadedcross(; s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real,
    ntri::Int, orient::Real=0.0, units::PSSFSSLength, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    structuredtri = haskey(kwargs, :structuredtri) ? kwargs[:structuredtri] : true
    fufp = haskey(kwargs, :fufp) ? kwargs[:fufp] : structuredtri

    if structuredtri
        return loadedcross_structured(; fufp, s1, s2, L1, L2, w, ntri, orient, units, kwarg...)
    else
        return loadedcross_unstructured(; fufp, s1, s2, L1, L2, w, ntri, orient, units, kwarg...)
    end
end


"""
    loadedcross_unstructured(;s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that
contains the triangulation for a "loaded cross" type of geometry, using an unstructured 
triangulation.  The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, 
and `fr` properly initialized.


The following (very poor) "ascii art" attempts to show
the definitions of the geometrical parameters `L1`, `L2` and `w`.
Note that the structure is supposed to be symmetrical wrt reflections
about its horizontal and vertical centerlines, and wrt reflections through a line oriented
at a 45 degree angle wrt the x-axis.


     ^                 ----------------
     |                 |  _________   |
     |                 |  |       |   |
     |                 |  |       |   |
     |                 |  |    -->|   |<--- W
     |                 |  |       |   |
     |                 |  |       |   |
     |     ------------   |       |   -------------
     |     |  |-----------|       |------------|  |
     |     |  |                                |  |
     L1    |  |                                |  |
     |     |  |                                |  |
     |     |  |                                |  |
     |     |  ------------          ------------  |
     |     |-----------   |        |  ------------|
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |________|  |
     |                 |              |
     V                 ----------------
    
                       <---- L2 ------>
    
# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `s1` and `s2`:  2-vectors containing the unit cell lattice vectors.
- `L1`,`L2`,`w`: Geometrical parameters as defined above.  Note that it is permissible
   to specify `w ≥ L2/2` in which case a solid (i.e., singly-connected) cross will be 
   generated.  In that case the `L2` dimension will be used for the width of the cross pieces.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
  the actual number will likely be different.
    
## Optional arguments:
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to locate the initial
  vertex of the loaded cross.  The default is to locate the vertex on the
  positive x-axis.
$(optional_kwargs)           
"""
function loadedcross_unstructured(; s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real,
    ntri::Int, orient::Real=0.0, units::PSSFSSLength, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    @testpos(L1)
    @testpos(L2)
    @testpos(w)
    @testpos(ntri)
    (length(s1) == length(s2) == 2) || throw(ArgumentError("s1 and s2 must be 2-vectors"))

    # Initialization:
    nv = (2w < L2 ? 24 : 12) # Total number of vertices
    ρ₀ = 0.5 * (s1 + s2) # Calculate center of polygon.
    ρ = Array{SV2}(undef, nv)
    e1 = Array{Cint}(undef, nv)
    e2 = Array{Cint}(undef, nv)
    segmarkers = Array{Cint}(undef, nv)
    holes = Array{Cdouble}(undef, 2, 0)

    # Set up the (outer) polygon geometry:
    ρ[1] = SV2([L2 / 2, L2 / 2])
    ρ[2] = SV2([L1 / 2, ρ[1][2]])
    ρ[3] = SV2([ρ[2][1], -ρ[2][2]])
    ρ[4] = SV2([ρ[1][1], ρ[3][2]])
    ρ[5] = SV2([ρ[4][1], -L1 / 2])
    ρ[6] = SV2([-ρ[5][1], ρ[5][2]])
    ρ[7] = SV2([ρ[6][1], ρ[4][2]])
    ρ[8] = SV2([-ρ[3][1], ρ[7][2]])
    ρ[9] = SV2([ρ[8][1], ρ[1][2]])
    ρ[10] = SV2([ρ[7][1], ρ[9][2]])
    ρ[11] = SV2([ρ[10][1], -ρ[6][2]])
    ρ[12] = SV2([ρ[1][1], ρ[11][2]])
    e1[1:12] = 1:12
    e2[1:11] = 2:12
    e2[12] = 1
    segmarkers[1:12] .= 1
    areat = 2 * (L1 - L2) * L2 # total area for solid cross

    if 2w < L2
        #  Set up inner boundary for annulus:
        ρ[13] = ρ[1] .- w
        ρ[14] = ρ[2] .- w
        ρ[15] = ρ[3] + SV2([-w, w])
        ρ[16] = ρ[4] + SV2([-w, w])
        ρ[17] = ρ[5] + SV2([-w, w])
        ρ[18] = ρ[6] .+ w
        ρ[19] = ρ[7] .+ w
        ρ[20] = ρ[8] .+ w
        ρ[21] = ρ[9] + SV2([w, -w])
        ρ[22] = -ρ[16]
        ρ[23] = -ρ[17]
        ρ[24] = -ρ[18]
        e1[13:24] = 13:24
        e2[13:23] = 14:24
        e2[24] = 13
        segmarkers[13:24] .= 2
        holes = [holes ρ₀]
        areat -= 2 * (L1 - L2) * (L2 - 2w)  # Subract inner void
    end

    if orient ≠ 0
        s, c = sincosd(orient)
        rotmat = SA[c -s; s c]
        for n in eachindex(ρ)
            ρ[n] = rotmat * ρ[n]
        end
    end

    ρ = [t + ρ₀ for t in ρ] # Center on the unit cell

    # Set up call to meshsub
    areatri = areat / ntri
    points = convert(Matrix{Cdouble}, reshape(reinterpret(Cdouble, ρ), (2, length(ρ))))
    segments = convert(Matrix{Cint}, transpose(hcat(e1, e2)))
    sheet = meshsub(points=points, seglist=segments, segmarkers=segmarkers,
        holes=holes, area=areatri, ntri=ntri)

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "loadedcross"
    sheet.ξη_check = false
    sheet.units = units
    sheet.s₁ = s1
    sheet.s₂ = s2
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function

"""
    meander(;a::Real, b::Real, h::Real, w1::Real, w2::Real, ntri::Int,
                  units::PSSFSSLength, orient=0, kwarg...) --> sheet::RWGSheet

# Description:
Return a variable of type `RWGSheet` that contains the triangulation for 
a meanderline strip.  The returned `sheet` has the components `s₁`, `s₂`, 
`β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, and `fr` properly initialized.  
Geometrical parameters are shown in the following diagram:
 
      - - - - - - - - - - - - - - - - - - - - - - - - -             ^
     |                                                |             |
     |                                                |             |
     |                                                |             |
     |                                                |             |
     |                                                |             |
     |            <-------- a/2 ------->              |             |
     |               (center-to-center)               |             |
     |                                                |             |
     |          ----------------------------          |  ^    ^     b
     |          |                          |          |  w2   |     |
     |          |                          |          |  |    |     |
     |          | -----------------------  |          |  v    |     |
     |          | |                     |  |          |             |
     |       -->| |<--w1           w1-->|  |<--       |       h     |
     ------------ |                     |  ------------  ^          |
     |            |                     |             |  w2   |     |
     |            |                     |             |  |    |     |
     ------------ - - - - - - - - - - - ---------------  v    v     v
 
     <-------------------- a ------------------------->
 
 
`a` and `b` are unit cell dimensions.  `w1` and `w2` are the widths
   of the vertical and horizontal strips, resp. `h` is the total
   height of the meander. 

A nicer diagram:
![https://simonp0420.github.io/PSSFSS.jl/stable/assets/meanderdef.png](https://simonp0420.github.io/PSSFSS.jl/stable/assets/meanderdef.png)

# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `a`,`b`,`h`,`w1`, `w2`: Geometrical parameters as defined above.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles. 
  This is a guide, the actual number will likely be different.
    
## Optional arguments:
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to rotate the 
  meanderline orientation within the unrotated unit cell.  Nonzero values are
  allowed only when the unit cell is a square (i.e. `a` == `b`).  The only allowable
  values are positive or negative multiples of 90.
$(optional_kwargs)
"""
function meander(; a::Real, b::Real, h::Real, w1::Real, w2::Real, ntri::Int,
    units::PSSFSSLength, orient::Real=0.0, kwarg...)::RWGSheet

    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(a)
    @testpos(b)
    @testpos(h)
    @testpos(w1)
    @testpos(w2)
    @testpos(ntri)
    orient ≠ 0 && a ≠ b && error("Nonzero orient only allowed for square unit cell")
    morient = mod(orient, 360)
    morient ∈ [0, 90, 180, 270] || error("orient must be a multiple of 90")

    ρ = Array{SV2}(undef, 0)
    e1 = Array{Cint}(undef, 0)
    e2 = Array{Cint}(undef, 0)
    segmarkers = Array{Cint}(undef, 0)
    node = 0

    # Calculate chopping increment
    t1 = w2 * (a + 2w1) / (2 * (h - 2w2)^2)
    t2 = w1 / (h - 2w2)
    fny2 = sqrt(ntri / (4 * (t1 + t2)))
    ny2 = round(Int, fny2)
    ny2 ≤ 0 && (ny2 = 1)
    ny1 = round(Int, fny2 * w2 / (h - 2w2))
    ny1 ≤ 0 && (ny1 = 1)
    nx1 = round(Int, fny2 * (a - 2w1) / (4 * (h - 2w2)))
    nx1 ≤ 0 && (nx1 = 1)
    nx2 = round(Int, fny2 * w1 / (h - 2w2))
    nx2 ≤ 0 && (nx2 = 1)
    ntotal = 4 * (2ny1 * (nx1 + nx2) + nx2 * ny2) # Actual # of triangles

    # Triangulate first section:
    yoffset = (b - h) / 2
    Lx = (a / 2 - w1) / 2
    Ly = h - 2w2
    ρbl = SV2([0.0, yoffset])
    ρtr = ρbl + SV2([Lx, w2])
    sh1 = recttri(ρbl, ρtr, nx1, ny1)
    # Triangulate third section:
    ρbl = SV2([Lx, yoffset])
    ρtr = ρbl + SV2([w1, w2])
    sh2 = recttri(ρbl, ρtr, nx2, ny1)
    # Combine them:
    sh3 = combine(sh1, sh2, 'x', Lx)
    # Add to section 5, store result in sh2:
    ρbl = SV2([Lx, yoffset + w2])
    ρtr = ρbl + SV2([w1, Ly])
    sh1 = recttri(ρbl, ρtr, nx2, ny2)
    sh2 = combine(sh3, sh1, 'y', ρbl[2])
    #  Add to section 7, store result in sh3
    ρbl = SV2([Lx, yoffset + w2 + Ly])
    ρtr = ρbl + SV2([w1, w2])
    sh1 = recttri(ρbl, ρtr, nx2, ny1)
    sh3 = combine(sh2, sh1, 'y', ρbl[2])
    #  Add to sections 9 and 10, store result in sh2
    ρbl = SV2([Lx + w1, yoffset + w2 + Ly])
    ρtr = ρbl + SV2([2 * Lx, w2])
    sh1 = recttri(ρbl, ρtr, 2nx1, ny1)
    sh2 = combine(sh3, sh1, 'x', ρbl[1])
    #  Add to section 8, store result in sh3
    ρbl = SV2([ρtr[1], ρbl[2]])
    ρtr = ρbl + SV2([w1, w2])
    sh1 = recttri(ρbl, ρtr, nx2, ny1)
    sh3 = combine(sh2, sh1, 'x', ρbl[1])
    # Add to section 6, store result in sh2:
    ρbl = SV2([ρbl[1], yoffset + w2])
    ρtr = ρbl + SV2([w1, Ly])
    sh1 = recttri(ρbl, ρtr, nx2, ny2)
    sh2 = combine(sh3, sh1, 'y', ρtr[2])
    #  Add to section 4, store result in sh3
    ρbl = SV2([ρbl[1], yoffset])
    ρtr = ρbl + SV2([w1, w2])
    sh1 = recttri(ρbl, ρtr, nx2, ny1)
    sh3 = combine(sh2, sh1, 'y', ρtr[2])
    # Add to section 2, store result in sh2:
    ρbl = SV2([ρtr[1], yoffset])
    ρtr = ρbl + SV2([Lx, w2])
    sh1 = recttri(ρbl, ρtr, nx1, ny1)
    sheet = combine(sh3, sh1, 'x', ρbl[1])
    if morient == 90
        xform = (x, y) -> (b - y, x)
    elseif morient == 180
        xform = (x, y) -> (a - x, b - y)
    elseif morient == 270
        xform = (x, y) -> (y, a - x)
    end
    if morient ≠ 0
        sheet.ρ = [SV2(xform(ρ[1], ρ[2])) for ρ in sheet.ρ]
    end

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    sheet.style = "meander"
    sheet.units = units
    sheet.s₁ = SV2([a, 0.0])
    sheet.s₂ = SV2([0.00, b])
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.ξη_check = true

    return sheet

end # function

"""
    pecsheet()

Return a variable of type `RWGSheet` that contains a perfect electric conducting sheet (i.e. an "E-wall").

"""
function pecsheet()::RWGSheet
    sheet = RWGSheet()
    sheet.style = "NULL"
    sheet.class = 'E'
    return sheet
end # function

"""
    pmcsheet()

Return a variable of type `RWGSheet` that contains a perfect magnetic conducting sheet (i.e. an "H-wall").

"""
function pmcsheet()::RWGSheet
    sheet = RWGSheet()
    sheet.style = "NULL"
    sheet.class = 'H'
    return sheet
end # function

"""
    polyring(;s1::Vector, s2::Vector, a::Vector, b::Vector, sides::Int ,ntri::Int ,orient::Real, units::PSSFSSLength, kwargs...) --> RWGSheet

Return a variable of type `RWGSheet` that contains the triangulation for one or more concentric annular regions bounded by polygons.

# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `s1` and `s2`:  2-vectors containing the unit cell lattice vectors.
- `a` and `b`:  n-vectors (n>=1) of the same length providing the inner and outer radii, respectively of the polygonal rings.
  Entries in `a` and `b` must be strictly increasing, except for possibly `b[end]` as discussed 
  below. `b[i] > a[i]` ∀ `i ∈ 1:n`, except possibly `b[end]` as discussed below. 
  `a[1]` may be zero to denote a solid (non-annular) polygon as the first "ring".
  It is possible to let the outermost ring to extend completely to the unit cell boundary.  
  This is specified by setting `b[end]` < 0, in which case for unstructured meshes,
  `-b[end]` is interpreted as the number of edges along the shorter of the `s1` and `s2` lattice vectors.
- `sides`:  The number (>= 3) of polygon sides.
- `ntri`:  The desired total number of triangles distributed among all the annular regions. This is a guide, the actual number 
  will likely be different.
    
## Optional arguments:
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to locate the initial
           vertex of the polygonal rings.  The default is to locate the vertex on the
           positive x-axis.
- `structuredtri::Bool`: Defaults to `true` when `sides==4` and false otherwise. A `true` value is only
  allowed when `sides==4` and `s1` ⟂ `s2`.  If true, use a structured mesh for the triangulation.  If false,
  the unstructured mesh generator that was standard up to PSSFSS version 1.2 will be used. A structured 
  mesh can be analyzed more efficiently, but the number of triangles created by the unstructured
  mesh generator is usually closer to `ntri` than the number for the structured mesh generator.
$(optional_kwargs)
"""
function polyring(; s1::Vector, s2::Vector, a::Vector{<:Real}, b::Vector{<:Real},
    sides::Int, ntri::Int, units::PSSFSSLength,
    orient::Real=0.0, kwarg...)::RWGSheet
    kwargs = Dict{Symbol,Any}(kwarg)

    if sides == 4
        structuredtri = haskey(kwargs, :structuredtri) ? kwargs[:structuredtri] : true
    else
        structuredtri = haskey(kwargs, :structuredtri) ? kwargs[:structuredtri] : false
        structuredtri && b[end] < 0 &&
            throw(ArgumentError("structuredtri=true not not compatible for polyring with sides≠4 and b[end]<0"))
    end

    structuredtri && b[end] < 0 && !iszero(s1 ⋅ s2) && 
        throw(ArgumentError("structuredtri=true not not possible for nonrectangular unit cell with b[end]<0"))
    fufp = haskey(kwargs, :fufp) ? kwargs[:fufp] : structuredtri

    if structuredtri
        return polyring_structured(; fufp, s1, s2, a, b, sides, ntri, units, orient, kwarg...)
    else
        return polyring_unstructured(; fufp, s1, s2, a, b, sides, ntri, units, orient, kwarg...)
    end
end

function polyring_unstructured(; s1::Vector, s2::Vector, a::Vector{<:Real}, b::Vector{<:Real},
    sides::Int, ntri::Int, units::PSSFSSLength,
    orient::Real=0.0, kwarg...)::RWGSheet
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    @testpos(sides)
    @testpos(ntri)
    (length(s1) == length(s2) == 2) || throw(ArgumentError("s1 and s2 must have length 2"))


    length(a) ≠ length(b) && throw(ArgumentError("length(a) !== length(b)"))
    nring = length(a)
    for i in 1:nring
        if i < nring || (i == nring && b[nring] > 0)
            a[i] ≥ b[i] && throw(ArgumentError("a[$i] ≥ b[$i]"))
        end
    end
    for i in 1:nring-1
        b[i] ≥ a[i+1] && throw(ArgumentError("b[$i] ≥ a[$(i+1)]"))
        a[i+1] - a[i] ≤ 0 && throw(ArgumentError("Elements of a must be strictly increasing"))
        b[i+1] - b[i] ≤ 0 && i < nring - 1 &&
            throw(ArgumentError("All but final element of b must be strictly increasing"))
    end


    if b[nring] < 0
        fillcell = true  # outer ring extends all the way to unit cell boundaries.
        ns1 = Int(-b[nring])  # number of edges along s1 direction.
        norms1, norms2 = (norm(s1), norm(s2))
        if norms1 ≤ norms2
            ns2 = round(Int, ns1 * norms2 / norms1) # number of edges along s2 direction
        else
            ns2 = ns1  # given value is to be associated with shorter edge.
            ns1 = round(Int, ns2 * norms1 / norms2)
        end
    else
        fillcell = false # outer ring has finite width.
    end

    ρ₀ = 0.5 * (s1 + s2) # calculate center of polygon.
    α = 360 / sides

    # compute area of each ring and total area of all rings:
    area_factor = sides / 2 * sind(α)
    area = [(b[i]^2 - a[i]^2) * area_factor for i in 1:nring]
    if fillcell  # need to recompute outer ring's area:
        area[nring] = zhatcross(s1) ⋅ s2 - area_factor * a[nring]^2
    end
    areat = sum(area) # total area of all rings.
    areatri = areat / ntri # Desired area of a single triangle

    ρ = Array{SV2}(undef, 0)
    e1 = Array{Cint}(undef, 0)
    e2 = Array{Cint}(undef, 0)
    segmarkers = Array{Cint}(undef, 0)
    holes = Array{SV2}(undef, 0)
    boundary = 0
    node = 0
    for iring in 1:nring
        if a[iring] == 0
            #  solid polygon:
            boundary += 1
            for i = 1:sides
                node += 1
                push!(e1, node)
                push!(e2, node + 1)
                push!(segmarkers, boundary)
                push!(ρ, ρ₀ + b[iring] * SV2(reverse(sincosd(orient + (i - 1) * α))))
            end
            e2[end] -= sides
        else
            if iring == nring && fillcell
                # annulus with regular polygon as inner boundary and unit cell as outer:
                # outer boundary first:
                boundary += 1
                i1 = 1
                i2 = ns1 + 1
                nodesave = node + 1 # initial node of outer boundary
                for i in i1:i2  # bottom edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node + 1)
                    push!(segmarkers, boundary)
                    push!(ρ, (i - i1) / ns1 * s1)
                end
                i1 = ns1 + 2
                i2 = ns1 + ns2 + 1
                for i in i1:i2 # right edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node + 1)
                    push!(segmarkers, boundary)
                    push!(ρ, s1 + (i - i1 + 1) / ns2 * s2)
                end
                i1 = ns1 + ns2 + 2
                i2 = 2 * ns1 + ns2 + 1
                for i in i1:i2 # top edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node + 1)
                    push!(segmarkers, boundary)
                    push!(ρ, s1 + s2 - (i - i1 + 1) / ns1 * s1)
                end
                i1 = 2 * ns1 + ns2 + 2
                i2 = 2 * ns1 + 2 * ns2
                for i in i1:i2 # left edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node + 1)
                    push!(segmarkers, boundary)
                    push!(ρ, s2 - (i - i1 + 1) / ns2 * s2)
                end
                e2[end] = nodesave

                # now inner boundary:
                boundary += 1
                i1 = 2 * (ns1 + ns2) + 1
                i2 = i1 + sides - 1
                nodesave = node + 1 # first node of inner boundary
                for i in i1:i2
                    node += 1
                    push!(e1, node)
                    push!(e2, node + 1)
                    push!(segmarkers, boundary)
                    push!(ρ, ρ₀ + a[iring] * SV2(reverse(sincosd(orient + (i - i1) * α))))
                end
                e2[end] = nodesave
            else
                # regular polygonal annulus:
                for r in (b[iring], a[iring])
                    boundary += 1
                    nodesave = node + 1
                    for i in 1:sides
                        node += 1
                        push!(e1, node)
                        push!(e2, node + 1)
                        push!(segmarkers, boundary)
                        push!(ρ, ρ₀ + r * SV2(reverse(sincosd(orient + (i - 1) * α))))
                    end
                    e2[end] = nodesave
                end
            end
        end
    end

    # Calculation coordinates of "hole" points
    ρhole = Array{SV2}(undef, 0)
    a[1] > 0 && push!(ρhole, ρ₀)

    for i in 1:nring-1
        unitvector = SV2(reverse(sincosd(orient)))
        r = 0.5 * (b[i] + a[i+1])
        push!(ρhole, ρ₀ + r * unitvector)
    end

    # Set up call to meshsub
    points = convert(Matrix{Cdouble}, reshape(reinterpret(Cdouble, ρ), (2, length(ρ))))
    segments = convert(Matrix{Cint}, transpose(hcat(e1, e2)))
    if isempty(ρhole)
        holes = Array{Cdouble}(undef, 2, 0)
    else
        holes = convert(Matrix{Cdouble}, reshape(reinterpret(Cdouble, ρhole), (2, length(ρhole))))
    end
    sheet = meshsub(points=points, seglist=segments, segmarkers=segmarkers,
        holes=holes, area=areatri, ntri=ntri)

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "polyring"
    sheet.ξη_check = fillcell
    sheet.units = units
    sheet.s₁ = SV2(s1)
    sheet.s₂ = SV2(s2)
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function polyring_unstructured

"""
    rectstrip(;Lx::Real, Ly::Real, Nx::Int, Ny::Int, Px::Real, Py::Real, units::PSSFSSLength, kwargs...)

Return a variable of type `RWGSheet` that contains the triangulation for a rectangular strip.

# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `Lx` and `Ly`:  Lengths of the strip in the x and y directions.
- `Px` and `Py`:  Lengths (periods) of the rectangular unit cell in the x and y directions.
- `Nx` and `Ny`:  Number of line segments in the x and y directions, for dividing up the strip into
  rectangles, which are  triangulated by adding a diagonal to each rectangle.
    
## Optional arguments:
$(optional_kwargs)
"""
function rectstrip(; Lx::Real, Ly::Real, Nx::Int, Ny::Int, Px::Real, Py::Real, units::PSSFSSLength,
    kwarg...)::RWGSheet
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(Lx)
    @testpos(Ly)
    @testpos(Nx)
    @testpos(Ny)
    @testpos(Px)
    @testpos(Py)

    # Setup triangulation:
    x0 = 0.5 * (Px - Lx)  # Center strip in unit cell
    y0 = 0.5 * (Py - Ly)  # Center strip in unit cell
    rhobl = SV2([x0, y0])
    rhotr = SV2([x0 + Lx, y0 + Ly])
    sheet = recttri(rhobl, rhotr, Nx, Ny)
    sheet.style = "rectstrip"
    sheet.units = units

    sheet.s₁ = SV2([Px, 0.0])
    sheet.s₂ = SV2([0.0, Py])
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.ξη_check = (Lx == Px || Ly == Py)

    return sheet

end # function

"""
_makeregpoly(radius, sides, center = [0.,0.], orient = 0.0)

Create a LibGEOS regular polygon.

## Arguments:
* `radius`: Radius of vertices on polygon boundary, measured from `center`.
* `sides`: Number of sides of the boundary regular polygons.
* `center`: A 2-vector containing the coordinates of the polygon center.
* `orient`: Orientation angle of the first vertex wrt the center location.
"""
function _makeregpoly(radius, sides, center=[0.0, 0.0], orient=0.0)
    io = IOBuffer()
    write(io, "POLYGON((")
    for i in 0:sides
        s, c = sincosd(orient + i * 360 / sides)
        print(io, center[1] + radius * c, " ", center[2] + radius * s)
        i < sides && write(io, ",")
    end
    write(io, "))")
    s = String(take!(io))
    poly = LibGEOS.readgeom(s)
    return poly
end

"""
    _makering(a, b, sides, center = [0.,0.], orient = 0.0)

Create a LibGEOS annular regular polygon.

## Arguments:
* `a`, `b`: Radii of vertices on the ring's inner and outer boundaries, respectively.
* `sides`: Number of sides of the boundary regular polygons.
* `center`: A 2-vector containing the coordinates of the annulus center.
* `orient`: Orientation angle of the first vertex wrt the center location.
"""
function _makering(a, b, sides, center=SV2[0.0, 0.0], orient=0.0)
    return LibGEOS.difference(_makeregpoly(b, sides, center, orient), _makeregpoly(a, sides, center, orient))
end


"""
    _makewedge(center, radius, centerangle, wedgeangle; sides=20)

Make a LibGeos pie-shaped wedge centered on angle `centerangle` of specified wedge angle.
Angles are in degrees.
"""
function _makewedge(center, radius, centerangle, wedgeangle; sides=20)
    io = IOBuffer()
    print(io, "POLYGON((", center[1], " ", center[2], ",")
    for θ in range(start=centerangle - wedgeangle / 2, stop=centerangle + wedgeangle / 2, length=sides)
        s, c = sincosd(θ)
        print(io, center[1] + radius * c, " ", center[2] + radius * s, ",")
    end
    print(io, center[1], " ", center[2], "))")
    s = String(take!(io))
    poly = LibGEOS.readgeom(s)
    return poly
end

"""
    _makerect(center, len, centerangle, width)

Make a LibGeos rectangle centered on angle `centerangle` of specified width and length.
Angles are in degrees.
"""
function _makerect(center, len, centerangle, width)
    x1, y1 = 0.0, -width / 2
    x2, y3 = len, y1 + width
    prevertices = SV2[[x1, y1], [x2, y1], [x2, y3], [x1, y3], [x1, y1]]
    s, c = sincosd(centerangle)
    rotmat = SA[c -s; s c]
    io = IOBuffer()
    print(io, "POLYGON((")
    for (i, prevertex) in enumerate(prevertices)
        x, y = center + rotmat * prevertex
        print(io, x, " ", y)
        i < length(prevertices) && write(io, ",")
    end
    write(io, "))")
    s = String(take!(io))
    rect = LibGEOS.readgeom(s)
    return rect
end

function _getcoordinates(geosgeom)
    return GeoInterface.coordinates(geosgeom)
end

"""
    splitring(; s1, s2, a, b, sides, ntri, gapwidth, gapcenter, gapangle, units, kwargs...) --> RWGSheet

Return a variable of type `RWGSheet` similar to a `polyring` but with zero or more gaps in each concentric annular region.

# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `s1` and `s2`:  2-vectors containing the unit cell lattice vectors.
- `a` and `b`:  n-vectors (n>=1) of the same length providing the inner and outer radii, respectively of the polygonal rings.
  Entries in `a` and `b` must be positive and strictly increasing. `b[i] > a[i]` ∀ `i ∈ 1:n`.
- `sides`:  The number (>= 4) of polygon sides for the background regular annular polygon(s) from which the gaps are removed.
- `gapcenter`: A scalar or vector of angles in degrees that define the gap center angular location(s), measured counterclockwise.  
  A scalar implies that all rings have a gap in that same angular location.  If a vector, then it must have the 
  same length as `a` and `b`, with `gapcenter[m]` denoting the gap center location for the `m`th ring.
  However `gapcenter[m]` can be either a scalar (denoting a single gap) or an n-tuple (denoting n gaps 
  in the `m`th ring).
- `gapwidth`: A scalar or a vector of the same length as `a` and `b` containing the gap width(s) for each ring.
  A width of zero implies that the ring is not split (i.e. there is no gap).  If the `gapwidth` of all rings
  is zero, then the resulting geometry is similar to a `polyring`. If a ring is to have multiple gaps, then
  the widths of the gaps for that ring should be passed as a tuple.  For example, suppose there are three
  rings and the second ring has 2 gaps, with the others having a single gap.  Then `gapwidth = [0.5, (0.4, 0.6), 0.3]`
  would be an appropriately formatted input in this case. When `gapwidth` is specified, the gaps are
  implemented as if a rectangular region is removed from the annular polygonal rings. Note that only 
  one of `gapwidth` and `gapangle` can be specified.
- `gapangle`: A scalar or vector of the same length as `a` and `b` containing the angular widths of the gaps in degrees.
  As with `gapwidth`, for any rings with multiple gaps, the corresponding entry in `gapangle` should be a 
  tuple of the same length as the number of gaps for that ring. When `gapangle` is specified, the gap(s) 
  in the `m`th ring is/are formed as if pie-shaped wedge(s) with wedge angle(s) `gapangle[m]`, are 
  removed from the ring(s). The locations and sizes of the tuples in `gapangle` must agree with those 
  in `gapcenter`.  Note that only one of `gapangle` and `gapwidth` can be specified.
- `ntri`:  The desired total number of triangles distributed among all the annular regions. This is a guide, the actual number 
  will likely be different.

## Optional arguments:
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to locate the initial
  vertex of the polygonal rings.  The default is to locate the vertex on the
  ray from the center parallel to the positive x-axis.
$(optional_kwargs)
"""
function splitring(;
    s1::Vector{<:Real},
    s2::Vector{<:Real},
    a::Vector{<:Real},
    b::Vector{<:Real},
    sides::Int,
    ntri::Int,
    units::PSSFSSLength,
    gapcenter::Union{Real,Vector,Nothing}=nothing,
    gapwidth::Union{Real,Vector,Nothing}=nothing,
    gapangle::Union{Real,Vector,Nothing}=nothing,
    orient::Real=0.0,
    kwarg...)::RWGSheet

    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    nring = length(a)
    isngc, isngw, isnga = isnothing.((gapcenter, gapwidth, gapangle))
    if !isngc
        if (isngw && isnga) || (!isngw && !isnga)
            throw(ArgumentError("Exactly one of gapwidth or gapangle must be specified"))
        end
        gcenter = isa(gapcenter, Real) ? fill(gapcenter, length(a)) : gapcenter
        if isngw
            gaorw = isa(gapangle, Real) ? fill(gapangle, length(a)) : gapangle
            usewedge = true
        else
            gaorw = isa(gapwidth, Real) ? fill(gapwidth, length(a)) : gapwidth
            usewedge = false
        end
        userect = !usewedge
        gaps = true
    else
        !isngw && throw(ArgumentError("gapwidth cannot be specified without gapcenter"))
        !isnga && throw(ArgumentError("gapangle cannot be specified without gapcenter"))
        gaps = false
    end

    all(x -> length(x) == nring, (b, gcenter, gaorw)) ||
        throw(ArgumentError("Incompatible lengths for a, b, gapcenter, gapwidth, or gapangle"))
    for (gc, gaw) in zip(gcenter, gaorw)
        length(gc) == length(gaw) ||
            throw(ArgumentError("Incompatible gapcenter and gap width/angle"))
    end
    sides ≥ 3 || throw(ArgumentError("Number of sides must be 3 or more"))
    @testpos(ntri)
    @testpos(a)
    @testpos(b)
    @testpos(b - a)
    (length(s1) == length(s2) == 2) || throw(ArgumentError("s1 and s2 must have length 2"))

    for i in 1:nring
        a[i] ≥ b[i] && throw(ArgumentError("a[$i] ≥ b[$i]"))
    end
    for i in 1:(nring-1)
        ip1 = i + 1
        b[i] ≥ a[ip1] && throw(ArgumentError("b[$i] ≥ a[$(ip1)]"))
        a[ip1] - a[i] ≤ 0 && throw(ArgumentError("Elements of a must be strictly increasing"))
        b[ip1] - b[i] ≤ 0 && throw(ArgumentError("Elements of b must be strictly increasing"))
    end

    ρ₀ = 0.5 * (s1 + s2) # calculate center of polygon.
    msdata = MeshsubData()
    for iring in 1:nring
        ring = _makering(a[iring], b[iring], sides, ρ₀, orient)
        for (gapcen, gapspec) in zip(gcenter[iring], gaorw[iring])
            gapspec == 0 && continue
            if usewedge
                poly = _makewedge(ρ₀, 1.2 * b[iring], gapcen, gapspec)
            else
                poly = _makerect(ρ₀, 1.2 * b[iring], gapcen, gapspec)
            end
            ring = LibGEOS.difference(ring, poly)
        end
        _add_libgeos_geom!(msdata, ring, ρ₀)
    end # for iring 

    # Set up call to meshsub
    areatri = msdata.area / ntri # Desired area of a single triangle
    ρ, e1, e2, ρhole, segmarkers = msdata.ρ, msdata.e1, msdata.e2, msdata.holes, msdata.segmarkers
    points = convert(Matrix{Cdouble}, reshape(reinterpret(Cdouble, ρ), (2, length(ρ))))
    segments = convert(Matrix{Cint}, transpose(hcat(e1, e2)))
    if isempty(ρhole)
        holes = Array{Cdouble}(undef, 2, 0)
    else
        holes = convert(Matrix{Cdouble}, reshape(reinterpret(Cdouble, ρhole), (2, length(ρhole))))
    end
    sheet = meshsub(points=points, seglist=segments, segmarkers=segmarkers,
        holes=holes, area=areatri, ntri=ntri)

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "splitring"
    sheet.ξη_check = false
    sheet.units = units
    sheet.s₁ = SV2(s1)
    sheet.s₂ = SV2(s2)
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function splitring

"""
    _tritri(w::Real, nw::Int) -> sh::RWGSheet

Create a variable of type `RWGSheet` that contains the triangulation for 
a isosceles right triangle.  The apex is located at the origin and the base is bisected by the positive x-axis.
The base length is w and the triangulation is generated by dividing up the triangular region into squares and adding 
a diagonal edge across each square (an exception is the triangle adjacent to the vertex of the large triangle in the case
where nw is odd). The fields `ρ`, `e1`, `e2`, `fv`, and `fe` are properly initialized upon return.
"""
function _tritri(w::Real, nw::Int)::RWGSheet
    even = mod(nw, 2) == 0
    if even
        k = nw ÷ 2
        nodecount = (k + 1)^2
        facecount = k * (k + 1)
    else
        k = (nw + 1) ÷ 2
        nodecount = k * (k + 1) + 1
        facecount = k * k
    end
    sh = RWGSheet()
    dxy = w / nw
    points = zeros(2, nodecount)
    # Set the node coordinates:
    node = 0
    if even
        for (i, x) in enumerate(range(0.0, w / 2, step=dxy))
            for y in range(-(i - 1) * dxy, (i - 1) * dxy, step=dxy)
                node += 1
                points[1:2, node] .= x, y
            end
        end
    else
        for mx in 0:k
            x = max(0.0, (2mx - 1) / 2 * dxy)
            for y in range(-x, x, step=dxy)
                node += 1
                points[1:2, node] .= x, y
            end
        end
    end
    node == nodecount || error("Node miscount")
    area = dxy^2
    astr = @sprintf("%.14f", area)
    switches = "Da$(astr)q30.0QeYY"
    switches = "Da$(astr)q30.0QeYY"
    seglist = Array{Int,2}(undef, 2, 0)
    sh = meshsub(; points, seglist, area, ntri=facecount, switches)::RWGSheet
    return sh
end

function rotationmat(θ)
    s, c = sincosd(θ)
    return SA[c -s; s c]
end

include("structuredtri.jl") # Code for structured meshes for loadedcross, jerusalemcross

end # module
