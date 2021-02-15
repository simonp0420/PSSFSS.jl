""" 
    Substrate

 Defines the derived type Layer for a substrate layer and
 Gblock for a generalized scattering matrix block entity.
"""
module Substrate

using ..PSSFSSLen  # For length units and ustrip only
using StaticArrays: MVector

export Layer, Gblock, TEorTM

@enum TEorTM TE=1 TM=2

"""
    Layer <: Any

An instance of the `Layer` type represents a single dielectric layer of the physical structure.
It contains the electrical properties of the dielectric layer. For layers not included
in a Gblock, an instance of Layer also specifies the periodicity (via the
reciprocal lattice vectors) and stores the mode constants for the Floquet modes present 
in the layer. 

    Layer(width::0u"mm", ϵᵣ=1.0, tanδ=0.0, μᵣ=1.0, mtanδ=0.0)

Create a Layer instance with the specified electrical properties.  All arguments
are optional keyword arguments with default values as shown above. They can be 
supplied in any order.

# Arguments
- `width`: The layer width (i.e. thickness) expressed as a 
        [`Unitful`](https://github.com/timholy/Unitful.jl) 
        length quantity. For convenience the following unit
        suffixes are exported by this module: `m`, `cm`, `mil`, `inch`,
        so one can specify, e.g., `width=20mil`.  Note that `width` can 
        be negative for the first and/or final layer of the composite structure,
        which has the effect of shifting the phase reference plane towards
        the interior of the composite structure.  This is sometimes needed 
        when interfacing with other programs, such as TEP file generation for
        Ticra's `GRASP` program.
- `ϵᵣ<:Real`: Relative permittivity of the dielectric.
- `tanδ<:Real`: Loss tangent (electrical) of the dielectric.
- `μᵣ<:Real`: Relative permeability of the dielectric.
- `mtanδ<:Real`: Loss tangent (magnetic) of the dielectric.

"""
mutable struct Layer
    name::String
    ϵᵣ::ComplexF64  # Relative permittivity
    μᵣ::ComplexF64  # Relative permeability
    user_width::Unitful.Length # Layer thickness in Unitful units
    width::Float64  # (thickness in meters without attached units)
    #  Mode indices, sorted in order of increasing real part, or decreasing 
    #  imaginary part, of γ, the mode attenuation constant.  M is the x 
    #  index, N is the y index, and P = TE or TM:
    P::Vector{TEorTM}
    M::Vector{Int}
    N::Vector{Int}
    β₁::MVector{2,Float64}  # Reciprocal lattice vector (1/m)
    β₂::MVector{2,Float64}  # Reciprocal lattice vector (1/m)
    # Transverse portion of modal propagation vectors (radians/meter):
    β::Vector{MVector{2,Float64}}
    γ::Vector{ComplexF64} # Modal attenuation constants (nepers/meter)
    Y::Vector{ComplexF64} # Modal admittances (Siemens)
    c::Vector{ComplexF64} # Modal normalization constants (volts/meter)
    tvec::Vector{MVector{2,Float64}} # Modal unit electric field vectors (unitless)
    
    # Interior constructor:
    function Layer(;name="Layer", width::Unitful.Length=0u"mm", ϵᵣ::Real=1.0,
                   tanδ::Real=0.0, μᵣ::Real=1.0, mtanδ::Real=0.0)
                    cϵᵣ = ϵᵣ * complex(1.0, -tanδ)
                    cμᵣ = μᵣ * complex(1.0, -mtanδ) 
                    new(name, cϵᵣ, cμᵣ, width, float(ustrip(u"m", width)),
                        TEorTM[], Int[], Int[], MVector{2}([0.,0.]), MVector{2}([0.,0.]), 
                        MVector{2,Float64}[], ComplexF64[], ComplexF64[], ComplexF64[], 
                        MVector{2,Float64}[])
    end # function
end # struct

Base.:(==)(l1::Layer, l2::Layer) = 
             all(f -> getfield(l1, f) == getfield(l2, f), 1:nfields(l1))

Base.show(::IO, ::MIME"text/plain", l::Layer) =
    println(l.name, ": width=", l.user_width, ", ϵᵣ=", real(l.ϵᵣ), ", tanδ=", -imag(l.ϵᵣ)/real(l.ϵᵣ),
            ", μᵣ=", real(l.μᵣ), ", mtanδ=", -imag(l.μᵣ)/real(l.μᵣ), )


"""
    Gblock <: Any

 
     Gblock(rng::UnitRange{Int}, j::Int)

`Gblock` constructor.
The `Gblock` type is used to represent GSMblocks. A GSMblock is a contiguous portion
of the composite FSS structure for which a single GSM (generalized scattering
matrix) is defined and computed. A GSMblock may consist of a single dielectric 
interface plane (with or without FSS sheet present), or it may consist
of multiple, adjacent interface planes and the intervening dielectric 
layers.  In the latter case, there must be an FSS sheet present at exactly
one of the interface planes within the GSMblock. Note
that interface plane k is located between layers k and k+1, so that
there are N-1 interface planes in a FSS structure consisting of N
dielectric layers.
   
# Arguments
- `rng`: Specifies the interface planes contained in the `Gblock`. `rng` must
         not be empty, and its initial and final values must be between `1`
         and `N-1` inclusive, where `N` is the number of dielectric layers.
- `j`:  Interface plane number containing the FSS sheet.  If the Gblock does
        not contain a sheet then `j` should be set to `0`.
"""
mutable struct Gblock
    rng::UnitRange{Int}  # Indices of the interface planes in this Gblock.
    j::Int   # Index of interface plane containing FSS sheet or 0 IF no sheet.

    function Gblock(rng::UnitRange{Int}, j::Int)
        length(rng) > 0 || error("rng=$rng must have positive length")
        if j < 0
            error("negative integer j = $j")
        elseif j == 0 
            length(rng) ≠ 1 && error("length(rng) must be 1 when j==0")
        else 
            j ∈ rng || error("j=$j is not a member of rng=$rng")
        end
        new(rng, j)
    end # function

end # struct

end # module
