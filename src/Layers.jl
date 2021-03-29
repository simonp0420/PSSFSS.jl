""" 
    Layers

 Defines the `Layer` type for a substrate layer.
"""
module Layers

using ..PSSFSSLen  # For length units and ustrip only
using StaticArrays: SVector

export Layer, Gblock, TEorTM, TE, TM

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
    β₁::SVector{2,Float64}  # Reciprocal lattice vector (1/m)
    β₂::SVector{2,Float64}  # Reciprocal lattice vector (1/m)
    # Transverse portion of modal propagation vectors (radians/meter):
    β::Vector{SVector{2,Float64}}
    γ::Vector{ComplexF64} # Modal attenuation constants (nepers/meter)
    Y::Vector{ComplexF64} # Modal admittances (Siemens)
    c::Vector{ComplexF64} # Modal normalization constants (volts/meter)
    tvec::Vector{SVector{2,Float64}} # Modal unit electric field vectors (unitless)
    
    # Interior constructor:
    function Layer(;name="Layer", width::Unitful.Length=0u"mm", ϵᵣ::Real=1.0,
                   tanδ::Real=0.0, μᵣ::Real=1.0, mtanδ::Real=0.0)
                    cϵᵣ = ϵᵣ * complex(1.0, -tanδ)
                    cμᵣ = μᵣ * complex(1.0, -mtanδ) 
                    new(name, cϵᵣ, cμᵣ, width, float(ustrip(u"m", width)),
                        TEorTM[], Int[], Int[], SVector{2}([0.,0.]), SVector{2}([0.,0.]), 
                        SVector{2,Float64}[], ComplexF64[], ComplexF64[], ComplexF64[], 
                        SVector{2,Float64}[])
    end # function
end # struct

Base.:(==)(l1::Layer, l2::Layer) = 
             all(f -> getfield(l1, f) == getfield(l2, f), 1:nfields(l1))

Base.show(::IO, ::MIME"text/plain", l::Layer) =
    println(l.name, ": width=", l.user_width, ", ϵᵣ=", real(l.ϵᵣ), ", tanδ=", -imag(l.ϵᵣ)/real(l.ϵᵣ),
            ", μᵣ=", real(l.μᵣ), ", mtanδ=", -imag(l.μᵣ)/real(l.μᵣ), ", ", length(l.P), " modes" )




end # module
