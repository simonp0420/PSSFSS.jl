module ZhatCross
export ẑ

using StaticArrays: SVector
import LinearAlgebra.×

struct UnitVector{N} end # N = 1, 2, or 3 denotes x, y, or z unit vector
const ẑ = UnitVector{3}()

×(::UnitVector{3}, t::SVector{2,T}) where {T} = SVector{2,T}(-t[2], t[1])
×(t::SVector{2,T}, ::UnitVector{3}) where {T} = SVector{2,T}(t[2], -t[1])
×(::UnitVector{3}, t::SVector{3,T}) where {T} = SVector{3,T}(-t[2], t[1], zero(T))
×(t::SVector{3,T}, ::UnitVector{3}) where {T} = SVector{3,T}(t[2], -t[1], zero(T))

@inline function ×(::UnitVector{3}, t::AbstractVector)
    length(t) == 2 && return [-t[2], t[1]]
    length(t) == 3 && return [-t[2], t[1], zero(eltype(t))]
    throw(ArgumentError("Vector length not 2 or 3"))
end

@inline function ×(t::AbstractVector, ::UnitVector{3})
    length(t) == 2 && return [t[2], -t[1]]
    length(t) == 3 && return [t[2], -t[1], zero(eltype(t))]
    throw(ArgumentError("Vector length not 2 or 3"))
end

end # module