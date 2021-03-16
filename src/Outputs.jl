module Outputs
export @outputs, Output

using LinearAlgebra: ⋅, norm
using ..Constants: c₀, twopi
using ..GSMs: GSM
using ..Layers: TEorTM, TE, TM
using ..Elements: s₁s₂2β₁β₂
using ..Modes: zhatcross
using Unitful
using StaticArrays: @SVector

@enum HorV H=1 V=2
@enum RorL R=1 L=2

using ..Sheets: SV2
using ..PSSFSSLen: PSSFSSLength

SteerType = Union{NamedTuple{(:ψ₁, :ψ₂), Tuple{Float64, Float64}}, 
                  NamedTuple{(:θ, :ϕ), Tuple{Float64, Float64}}}

struct Output
    gsm::GSM
    steering::SteerType
    β⃗₀₀::SV2 # radians/meter
    FGHz::Float64
    ϵᵣin::ComplexF64
    μᵣin::ComplexF64
    unitsin::PSSFSSLength
    s₁in::SV2  # in units of unitsin
    s₂in::SV2  # in units of unitsin
    ϵᵣout::ComplexF64
    μᵣout::ComplexF64
    unitsout::PSSFSSLength
    s₁out::SV2 # in units of unitsout
    s₂out::SV2 # in units of unitsout
end



struct Outout{F <: Function} 
    f::F
    label::String
end
(oo::Outout)(o::Output) = oo.f(o)  # Make it a functor


"""
    getsijmn(i::Int,j::Int,m,n,o::Output)

Obtain the `(m,n)` entry of the `(i,j)` partition of `o.gsm`.  Note that 
`m` and `n` can be either integers or `enums` of type `TEorTEM`, `RorL`, or `HorV`.
In either of the last two cases, the GSM is modified appropriately as described in 
Chapter 8 of the theory documentation.
"""
@inline function getsijmn(i::Int, j::Int, m::Int, n::Int, o::Output)
    o.gsm[i,j][m,n]
end


@inline function getsijmn(i::Int, j::Int, m::Integer, n, o::Output)
    (view(o.gsm[i,j], 1:2, 1:2) * sourcemat(j,n,o))[m, Int(n)]
end


@inline function getsijmn(i::Int, j::Int, m, n::Integer, o::Output)
    (obsmat(i,m,o) * view(o.gsm[i,j], 1:2, 1:2))[Int(m), n]
end

@inline function getsijmn(i::Int, j::Int, m, n, o::Output)
    (obsmat(i,m,o) * view(o.gsm[i,j], 1:2, 1:2) * sourcemat(j,n,o))[Int(m), Int(n)]
end


"""
sourcemat(j, n::union{HorV,RorL}, o::Output)

Compute a 2×2 transformation matrix which when used to right-multiply `o.gsm[i,j]` performs 
a basis change for the polarization basis vectors from TE/TM to either CP (circular polarization)
or horizontal/vertical, as determined by the type of `n`.
"""
function sourcemat(j::Int, n::HorV, o::Output)
    (θ1inc, ϕ1inc) = θϕ(o)
    if j == 1 # Region 1 incidence
        (θ,ϕ) = (θ1inc, ϕ1inc) # Eqs. (8.5a) and (8.5b)
    else # Region N incidence
        n1 = sqrt(real(o.ϵᵣin * o.μᵣin))
        n2 = sqrt(real(o.ϵᵣout * o.μᵣout))
        θ = asind(n1/n2*sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc + 180
    end 
    (ĥ3, v̂3) = ĥv̂(θ, ϕ)
    ĥ = @view ĥ3[1:2]  # Only need x and y components due to dot product later
    v̂ = @view v̂3[1:2]  # Only need x and y components due to dot product later
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([1.0, 0.0]) : o.β⃗₀₀ / β₀₀
    t̂₁ = zhatcross(β̂₀₀)
    t̂₂ = β̂₀₀
    ct = cosd(θ)
    mat = [ĥ⋅t̂₁     v̂⋅t̂₁
           ĥ⋅t̂₂/ct  v̂⋅t̂₂/ct]
    return mat
end
function sourcemat(j::Int, n::RorL, o::Output)
    (θ1inc, ϕ1inc) = θϕ(o)
    if j == 1 # Region 1 incidence
        (θ,ϕ) = (θ1inc, ϕ1inc) # Eqs. (8.5a) and (8.5b)
        sgn = 1
    else # Region N incidence
        n1 = sqrt(real(o.ϵᵣin * o.μᵣin))
        n2 = sqrt(real(o.ϵᵣout * o.μᵣout))
        θ = asind(n1/n2*sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc + 180
        sgn = -1
    end 
    (ĥ, v̂) = ĥv̂(θ, ϕ)
    L̂ = view((ĥ + sgn*im*v̂)/√2, 1:2) # Only need x and y components due to dot product later
    R̂ = view((ĥ - sgn*im*v̂)/√2, 1:2)
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([1.0, 0.0]) : o.β⃗₀₀ / β₀₀
    t̂₁ = zhatcross(β̂₀₀)
    t̂₂ = β̂₀₀
    ct = cosd(θ)
    mat = [t̂₁⋅R̂     t̂₁⋅L̂
           t̂₂⋅R̂/ct  t̂₂⋅L̂/ct] # Dot products reversed to avoid conjugation
    return mat
end


"""
obsmat(i::Int, n::union{HorV,RorL}, o::Output)

Compute a 2×2 transformation matrix which when used to left-multiply `o.gsm[i,j]` performs 
a basis change for the polarization basis vectors from TE/TM to either CP (circular polarization)
or horizontal/vertical, as determined by the type of `n`.
"""
function obsmat(i::Int, n::HorV, o::Output)
    (θ1inc, ϕ1inc) = θϕ(o)
    if i == 1 # Region 1 reflection
        (θ, ϕ) = (θ1inc, ϕ1inc+180) # Eq. (8.5c)
        sgn = 1
    else # Region N reflection
        n1 = sqrt(real(o.ϵᵣin * o.μᵣin))
        n2 = sqrt(real(o.ϵᵣout * o.μᵣout))
        θ = asind(n1/n2*sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc  # Eq. (8.5e)
        sgn = -1
    end 
    (ĥ, v̂) = ĥv̂(θ, ϕ)
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([1.0, 0.0]) : o.β⃗₀₀ / β₀₀
    t̂₁2 = zhatcross(β̂₀₀)
    t̂₁ = @SVector([t̂₁2[1], t̂₁2[2], 0.0])
    t̂₂ = @SVector([β̂₀₀[1], β̂₀₀[2], sgn*tand(θ)]) # term from Eqs. (8.20)
    mat = [ĥ⋅t̂₁  ĥ⋅t̂₂
           v̂⋅t̂₁  v̂⋅t̂₂]
    return mat
end

function obsmat(i::Int, n::RorL, o::Output)
    (θ1inc, ϕ1inc) = θϕ(o)
    if i == 1 # Region 1 reflection
        (θ,ϕ) = (θ1inc, ϕ1inc+180) # Eq. (8.5c)
        sgn = 1
    else # Region N reflection
        n1 = sqrt(o.ϵᵣin * o.μᵣin)
        n2 = sqrt(o.ϵᵣout * o.μᵣout)
        θ = asind(n1/n2*sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc  # Eq. (8.5e)
        sgn = -1
    end 
    (ĥ, v̂) = ĥv̂(θ, ϕ)
    L̂ = (ĥ - sgn*im*v̂)/√2
    R̂ = (ĥ + sgn*im*v̂)/√2
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([1.0, 0.0]) : o.β⃗₀₀ / β₀₀
    t̂₁2 = zhatcross(β̂₀₀)
    t̂₁ = @SVector([t̂₁2[1], t̂₁2[2], 0.0])
    t̂₂ = @SVector([β̂₀₀[1], β̂₀₀[2], sgn*tand(θ)]) # term from Eqs. (8.20)
    mat = [R̂⋅t̂₁  R̂⋅t̂₂
           L̂⋅t̂₁  L̂⋅t̂₂]
    return mat
end


getSIJMN(i,j,m,n) = Outout("S$(i)$(j)($(m),$n)") do o
    getsijmn(i, j, m, n, o)
end
S11(m,n) = getSIJMN(1,1,m,n)
S12(m,n) = getSIJMN(1,2,m,n)
S21(m,n) = getSIJMN(2,1,m,n)
S22(m,n) = getSIJMN(2,2,m,n)

smag(i,j,m,n) = Outout("S$(i)$(j)MAG($m,$n)") do o
    abs(getsijmn(i,j,m,n,o))
end
S11MAG(m,n) = smag(1,1,m,n)
S12MAG(m,n) = smag(1,2,m,n)
S21MAG(m,n) = smag(2,1,m,n)
S22MAG(m,n) = smag(2,2,m,n)

sdb(i,j,m,n) =  Outout("S$(i)$(j)DB($m,$n)") do o
    10*log10(abs2(getsijmn(i,j,m,n,o)))
end
S11DB(m,n) = sdb(1,1,m,n)
S12DB(m,n) = sdb(1,2,m,n)
S21DB(m,n) = sdb(2,1,m,n)
S22DB(m,n) = sdb(2,2,m,n)

sang(i,j,m,n) = Outout("S$(i)$(j)ANG($m,$n)") do o
    rad2deg(angle(getsijmn(i,j,m,n,o)))
end
S11ANG(m,n) = sang(1,1,m,n)
S12ANG(m,n) = sang(1,2,m,n)
S21ANG(m,n) = sang(2,1,m,n)
S22ANG(m,n) = sang(2,2,m,n)

sreal(i,j,m,n) = Outout("S$(i)$(j)REAL($m,$n)") do o
    real(getsijmn(i,j,m,n,o))
end
S11REAL(m,n) = sreal(1,1,m,n)
S12REAL(m,n) = sreal(1,2,m,n)
S21REAL(m,n) = sreal(2,1,m,n)
S22REAL(m,n) = sreal(2,2,m,n)

simag(i,j,m,n) = Outout("S$(i)$(j)IMAG($m,$n)") do o
    imag(getsijmn(i,j,m,n,o))
end
S11IMAG(m,n) = simag(1,1,m,n)
S12IMAG(m,n) = simag(1,2,m,n)
S21IMAG(m,n) = simag(2,1,m,n)
S22IMAG(m,n) = simag(2,2,m,n)

ΔIPD21 = Outout("ΔIPD") do o
    rad2deg(angle(getsijmn(2,1,1,1,o)/getsijmn(2,1,2,2,o)))
end
ΔIPD12 = Outout("ΔIPD") do o
    rad2deg(angle(getsijmn(1,2,1,1,o)/getsijmn(1,2,2,2,o)))
end

ΔIL21 = Outout("ΔIL") do o
    10*log10(abs2(getsijmn(2,1,1,1,o)/getsijmn(2,1,2,2,o)))
end
ΔIL12 = Outout("ΔIL") do o
    10*log10(abs2(getsijmn(1,2,1,1,o)/getsijmn(1,2,2,2,o)))
end

ardb(i,j,n) = Outout("ARdB$i$j($n)") do o
    jP = im * getsijmn(i,j,1,n,o)/getsijmn(i,j,2,n,o) # Modified Linear Pol. ratio
    Q = (1-jP)/(1+jP) # Circular polarization ratio
    absQ = abs(Q)
    absQ > 1 && (absQ = 1/absQ)
    ardb = 20*log10((1+absQ)/(1-absQ))
end
ARDB11(n) = ardb(1,1,n)
ARDB12(n) = ardb(1,2,n)
ARDB21(n) = ardb(2,1,n)
ARDB22(n) = ardb(2,2,n)


FGHZ = Outout("FGHZ") do o
    o.FGHz
end

FMHZ = Outout("FMHZ") do o
    o.FGHz * 1000
end

THETA = Outout("THETA") do o
    get(o.steering, :θ, NaN)
end
Θ = THETA

PHI = Outout("PHI") do o
    get(o.steering, :ϕ, NaN)
end
Φ = PHI

PSI1 = Outout("PSI1") do o
    get(o.steering, :ψ₁, NaN)
end
Ψ₁ = PSI1

PSI2 = Outout("PSI2") do o
    get(o.steering, :ψ₂, NaN)
end
Ψ₂ = PSI2


"""
    θϕ(o::Output) -> (θ, ϕ)

Return steering angles in degrees from `o`.  If `o` specifies `ψ₁` and `ψ₂` instead of angles,
then the latter are computed, using the Region 1 (input) periodicity and electrical parameters.
"""
function θϕ(o::Output)
    #=
    haskey(o.steering, :θ) && return (o.steering.θ, o.steering.ϕ)
    ψ₁, ψ₂ = o.steering
    units_per_meter = ustrip(Float64, o.unitsin, 1u"m")
    s₁, s₂ = [o.s₁in, o.s₂in] / units_per_meter
    β₁, β₂ = s₁s₂2β₁β₂(s₁,s₂)
    β₀₀ = (ψ₁ * β₁ + ψ₂ * β₂) / (2π)
    =#
    β₀₀² = o.β₀₀ ⋅ o.β₀₀
    β₀₀² == 0 && return (0.0, get(o.steering, ϕ, 0.0))
    k² = (twopi * FGHz * 1e9 / c₀)^2 * real(o.ϵᵣin * o.μᵣin)
    β₀₀² > k² && error("Cut-off dominant mode")
    kz = √(k² - β₀₀²)  # for out-going wave vector in Layer 1
    θ = acosd(kz/sqrt(k²))
    ϕ = atand(β₀₀[2], β₀₀[1])
    return (θ,ϕ)
end

"""
    ĥv̂(θ, ϕ)  

Compute Ludwig 3 unit vectors from spherical location vectors.  
"""
function ĥv̂(θ, ϕ)
    st,ct = sincosd(θ)
    sp,cp = sincosd(ϕ)
    x,y,z = st*cp, st*sp, ct
    r̂ = [st*cp, st*sp, ct]
    θ̂ = [ct*cp, ct*sp, -st]
    ϕ̂ = [-sp, cp, 0.0]
    ĥ = θ̂*cp - ϕ̂*sp
    v̂ = θ̂*cp + ϕ̂*sp
    ĥ, v̂
end


"""
    @outputs(args...)

Convert list of user output requests to a vector of functors that generate the requested
outputs when applied to an `Output` instance.  In the conversion process, replace
lower case letters with upper case.

### Examples

    julia> output = @outputs FGHz θ ϕ s11db(te,te) S11ang(Te,te)
    julia> output = @outputs FGHz theta phi s21db(R,H) ARdB21(H) ARdB11(v)
"""
macro outputs(args...)
    newargs = Any[]
    for (iarg,arg) in pairs(args)
        if arg isa Symbol
            push!(newargs,Symbol(uppercase(string(arg))))
        elseif arg isa Expr && arg.head == :call
            for (iarg2,arg2) in pairs(arg.args)
                arg2 isa Symbol && (arg.args[iarg2] = Symbol(uppercase(string(arg2))))
            end
            push!(newargs, arg)
        else
            error("Illegal @outputs construction")
        end
    end
    tuple([eval(a) for a in newargs]...)
end

end # module