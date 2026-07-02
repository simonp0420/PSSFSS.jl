module UnitVectors
export x̂, ŷ, ẑ

using StaticArrays: SVector
import LinearAlgebra: ×, dot

struct UnitVector{N} end # N = 1, 2, or 3 denotes x, y, or z unit vector
const x̂ = UnitVector{1}()
const ŷ = UnitVector{2}()
const ẑ = UnitVector{3}()

dot(::UnitVector{T}, v::AbstractVector) where {T} = v[T]
dot(v::AbstractVector, ::UnitVector{T}) where {T} = conj(v[T])


×(::UnitVector{3}, t::SVector{2,T}) where {T} = SVector{2,T}(-t[2], t[1])

×(::UnitVector{1}, t::SVector{3,T}) where {T} = SVector{3,T}(zero(T), -t[3], t[2])
×(::UnitVector{2}, t::SVector{3,T}) where {T} = SVector{3,T}(t[3], zero(T), -t[1])
×(::UnitVector{3}, t::SVector{3,T}) where {T} = SVector{3,T}(-t[2], t[1], zero(T))

@inline function ×(::UnitVector{1}, t::AbstractVector)
    length(t) == 3 && return [zero(eltype(t)), -t[3], t[2]]
    throw(ArgumentError("Vector length not 3"))
end

@inline function ×(::UnitVector{2}, t::AbstractVector)
    length(t) == 3 && return [t[3], zero(eltype(t)), -t[1]]
    throw(ArgumentError("Vector length not 3"))
end

@inline function ×(::UnitVector{3}, t::AbstractVector)
    length(t) == 2 && return [-t[2], t[1]]
    length(t) == 3 && return [-t[2], t[1], zero(eltype(t))]
    throw(ArgumentError("Vector length not 2 or 3"))
end

×(v::AbstractVector, û::UnitVector) = ×(û, -v)


end # module
