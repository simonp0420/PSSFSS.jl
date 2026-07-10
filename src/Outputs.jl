module Outputs
export @outputs, Result, append_result_data, read_result_file, extract_result_file, extract_result, res2fresnel, res2tep

using LinearAlgebra: ⋅, norm, ×
using ..UnitVectors: ẑ
using ..Constants: c₀, twopi
using ..GSMs: GSM
using ..Layers: TEorTM, TE, TM
using ..Elements: s₁s₂2β₁β₂
using Unitful
using StaticArrays: @SVector, @SMatrix, SMatrix
using JLD2: JLD2, jldopen
using FileIO: load
using TicraUtilities: TicraUtilities
using Printf: @printf
using Dates: now

@enum HorV H = 1 V = 2
@enum RorL R = 1 L = 2

using ..Sheets: SV2
using ..PSSFSSLen: PSSFSSLength

SteerType = Union{NamedTuple{(:ψ₁, :ψ₂),Tuple{Float64,Float64}},
    NamedTuple{(:θ, :ϕ),Tuple{Float64,Float64}}}

struct Result
    gsm::GSM
    steering::SteerType
    β⃗₀₀::SV2 # radians/meter
    FGHz::Float64
    ϵᵣin::ComplexF64
    μᵣin::ComplexF64
    β₁in::SV2  # radians/meter
    β₂in::SV2  # radians/meter
    ϵᵣout::ComplexF64
    μᵣout::ComplexF64
    β₁out::SV2 # radians/meter
    β₂out::SV2 # radians/meter
end

Base.show(io::IO, ::MIME"text/plain", r::Result) =
    print(io, "Result: ", r.FGHz, " GHz, ", r.steering, ", GSM", size(r.gsm.s12))



struct Outfun{F<:Function}
    f::F
    label::String
end
(o::Outfun)(r::Result) = o.f(r)  # Make it a functor

Base.show(io::IO, ::MIME"text/plain", o::Outfun) =
    print(io, "Outfun: ", o.label)

function Base.show(io::IO, ::MIME"text/plain", t::NTuple{N,Outfun} where {N})
    print("Outfun NTuple: (")
    for (i, o) in pairs(t)
        if i < length(t)
            print(io, o.label, ", ")
        else
            print(io, o.label, ")")
        end
    end
end


"""
    getsijmn(i::Int,j::Int,m,n,o::Result)

Obtain the `(m,n)` entry of the `(i,j)` partition of `o.gsm`.  Note that
`m` and `n` can be either integers or `enums` of type `TEorTEM`, `RorL`, or `HorV`.
In either of the last two cases, the GSM is modified appropriately as described in
Chapter 8 of the theory documentation.
"""
@inline function getsijmn(i::Int, j::Int, m::Union{Int,TEorTM}, n::Union{Int,TEorTM}, o::Result)
    o.gsm[i, j][Int(m), Int(n)]
end


@inline function getsijmn(i::Int, j::Int, m::Integer, n::Union{HorV,RorL}, o::Result)
    (SMatrix{2,2}(view(o.gsm[i, j], 1:2, 1:2)) * sourcemat(j, n, o))[m, Int(n)]
end


@inline function getsijmn(i::Int, j::Int, m::Union{HorV,RorL}, n::Integer, o::Result)
    (obsmat(i, m, o) * SMatrix{2,2}(view(o.gsm[i, j], 1:2, 1:2)))[Int(m), n]
end

@inline function getsijmn(i::Int, j::Int, m::Union{HorV,RorL}, n::Union{HorV,RorL}, o::Result)
    (obsmat(i, m, o) * SMatrix{2,2}(view(o.gsm[i, j], 1:2, 1:2)) * sourcemat(j, n, o))[Int(m), Int(n)]
end


"""
sourcemat(j, n::union{HorV,RorL}, o::Result)

Compute a 2×2 transformation matrix which when used to right-multiply `o.gsm[i,j]` performs
a basis change for the polarization basis vectors from TE/TM to either CP (circular polarization)
or horizontal/vertical, as determined by the type of `n`.
"""
function sourcemat(j::Int, n::HorV, o::Result)
    (θ1inc, ϕ1inc) = θϕ(o)
    if j == 1 # Region 1 incidence
        (θ, ϕ) = (θ1inc, ϕ1inc) # Eqs. (8.5a) and (8.5b)
    else # Region N incidence
        n1 = sqrt(real(o.ϵᵣin) * real(o.μᵣin))
        n2 = sqrt(real(o.ϵᵣout) * real(o.μᵣout))
        θ = asind(n1 / n2 * sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc + 180
    end
    (ĥ3, v̂3) = ĥv̂(θ, ϕ)
    ĥ = @view ĥ3[1:2]  # Only need x and y components due to dot product later
    v̂ = @view v̂3[1:2]  # Only need x and y components due to dot product later
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([cosd(ϕ1inc), sind(ϕ1inc)]) : o.β⃗₀₀ / β₀₀
    t̂₁ = ẑ × β̂₀₀
    t̂₂ = β̂₀₀
    ct = cosd(θ)
    mat = @SMatrix [ĥ⋅t̂₁    v̂⋅t̂₁
                    ĥ⋅t̂₂/ct v̂⋅t̂₂/ct]
    return mat
end
function sourcemat(j::Int, n::RorL, o::Result)
    (θ1inc, ϕ1inc) = θϕ(o)
    if j == 1 # Region 1 incidence
        (θ, ϕ) = (θ1inc, ϕ1inc) # Eqs. (8.5a) and (8.5b)
        sgn = 1
    else # Region N incidence
        n1 = sqrt(real(o.ϵᵣin) * real(o.μᵣin))
        n2 = sqrt(real(o.ϵᵣout) * real(o.μᵣout))
        θ = asind(n1 / n2 * sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc + 180
        sgn = -1
    end
    (ĥ, v̂) = ĥv̂(θ, ϕ)
    L̂ = view((ĥ + sgn * im * v̂) / √2, 1:2) # Only need x and y components due to dot product later
    R̂ = view((ĥ - sgn * im * v̂) / √2, 1:2)
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([cosd(ϕ1inc), sind(ϕ1inc)]) : o.β⃗₀₀ / β₀₀
    t̂₁ = ẑ × β̂₀₀
    t̂₂ = β̂₀₀
    ct = cosd(θ)
    mat = @SMatrix [t̂₁⋅R̂    t̂₁⋅L̂
                    t̂₂⋅R̂/ct t̂₂⋅L̂/ct] # Dot products reversed to avoid conjugation
    return mat
end


"""
obsmat(i::Int, n::union{HorV,RorL}, o::Result)

Compute a 2×2 transformation matrix which when used to left-multiply `o.gsm[i,j]` performs
a basis change for the polarization basis vectors from TE/TM to either CP (circular polarization)
or horizontal/vertical, as determined by the type of `n`.
"""
function obsmat(i::Int, n::HorV, o::Result)
    (θ1inc, ϕ1inc) = θϕ(o)
    if i == 1 # Region 1 reflection
        (θ, ϕ) = (θ1inc, ϕ1inc + 180) # Eq. (8.5c)
        sgn = 1
    else # Region N reflection
        n1 = sqrt(real(o.ϵᵣin) * real(o.μᵣin))
        n2 = sqrt(real(o.ϵᵣout) * real(o.μᵣout))
        θ = asind(n1 / n2 * sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc  # Eq. (8.5e)
        sgn = -1
    end
    (ĥ, v̂) = ĥv̂(θ, ϕ)
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([cosd(ϕ1inc), sind(ϕ1inc)]) : o.β⃗₀₀ / β₀₀
    t̂₁2 = ẑ × β̂₀₀
    t̂₁ = @SVector([t̂₁2[1], t̂₁2[2], 0.0])
    t̂₂ = @SVector([β̂₀₀[1], β̂₀₀[2], sgn * tand(θ)]) # term from Eqs. (8.20)
    mat = @SMatrix [ĥ⋅t̂₁ ĥ⋅t̂₂
                    v̂⋅t̂₁ v̂⋅t̂₂]
    return mat
end

function obsmat(i::Int, n::RorL, o::Result)
    (θ1inc, ϕ1inc) = θϕ(o)
    if i == 1 # Region 1 reflection
        (θ, ϕ) = (θ1inc, ϕ1inc + 180) # Eq. (8.5c)
        sgn = 1
    else # Region N reflection
        n1 = sqrt(real(o.ϵᵣin) * real(o.μᵣin))
        n2 = sqrt(real(o.ϵᵣout) * real(o.μᵣout))
        θ = asind(n1 / n2 * sind(θ1inc)) # Snell's law
        ϕ = ϕ1inc  # Eq. (8.5e)
        sgn = -1
    end
    (ĥ, v̂) = ĥv̂(θ, ϕ)
    L̂ = (ĥ - sgn * im * v̂) / √2
    R̂ = (ĥ + sgn * im * v̂) / √2
    β₀₀ = norm(o.β⃗₀₀)
    β̂₀₀ = (β₀₀ == 0) ? @SVector([cosd(ϕ1inc), sind(ϕ1inc)]) : o.β⃗₀₀ / β₀₀
    t̂₁2 = ẑ × β̂₀₀
    t̂₁ = @SVector([t̂₁2[1], t̂₁2[2], 0.0])
    t̂₂ = @SVector([β̂₀₀[1], β̂₀₀[2], sgn * tand(θ)]) # term from Eqs. (8.20)
    mat = @SMatrix [R̂⋅t̂₁ R̂⋅t̂₂
                    L̂⋅t̂₁ L̂⋅t̂₂]
    return mat
end


getSIJMN(i, j, m, n) =
    Outfun("S$(i)$(j)($(m),$n)") do o
        getsijmn(i, j, m, n, o)
    end
S11(m, n) = getSIJMN(1, 1, m, n)
S12(m, n) = getSIJMN(1, 2, m, n)
S21(m, n) = getSIJMN(2, 1, m, n)
S22(m, n) = getSIJMN(2, 2, m, n)

smag(i, j, m, n) =
    Outfun("S$(i)$(j)MAG($m,$n)") do o
        abs(getsijmn(i, j, m, n, o))
    end
S11MAG(m, n) = smag(1, 1, m, n)
S12MAG(m, n) = smag(1, 2, m, n)
S21MAG(m, n) = smag(2, 1, m, n)
S22MAG(m, n) = smag(2, 2, m, n)

sdb(i, j, m, n) =
    Outfun("S$(i)$(j)DB($m,$n)") do o
        10 * log10(abs2(getsijmn(i, j, m, n, o)))
    end
S11DB(m, n) = sdb(1, 1, m, n)
S12DB(m, n) = sdb(1, 2, m, n)
S21DB(m, n) = sdb(2, 1, m, n)
S22DB(m, n) = sdb(2, 2, m, n)

sang(i, j, m, n) =
    Outfun("S$(i)$(j)ANG($m,$n)") do o
        rad2deg(angle(getsijmn(i, j, m, n, o)))
    end
S11ANG(m, n) = sang(1, 1, m, n)
S12ANG(m, n) = sang(1, 2, m, n)
S21ANG(m, n) = sang(2, 1, m, n)
S22ANG(m, n) = sang(2, 2, m, n)

sreal(i, j, m, n) =
    Outfun("S$(i)$(j)REAL($m,$n)") do o
        real(getsijmn(i, j, m, n, o))
    end
S11REAL(m, n) = sreal(1, 1, m, n)
S12REAL(m, n) = sreal(1, 2, m, n)
S21REAL(m, n) = sreal(2, 1, m, n)
S22REAL(m, n) = sreal(2, 2, m, n)

simag(i, j, m, n) =
    Outfun("S$(i)$(j)IMAG($m,$n)") do o
        imag(getsijmn(i, j, m, n, o))
    end
S11IMAG(m, n) = simag(1, 1, m, n)
S12IMAG(m, n) = simag(1, 2, m, n)
S21IMAG(m, n) = simag(2, 1, m, n)
S22IMAG(m, n) = simag(2, 2, m, n)

ΔIPD21 = Outfun("ΔIPD") do o
    rad2deg(angle(getsijmn(2, 1, 1, 1, o) / getsijmn(2, 1, 2, 2, o)))
end
DIPD21 = ΔIPD21

ΔIPD12 = Outfun("ΔIPD") do o
    rad2deg(angle(getsijmn(1, 2, 1, 1, o) / getsijmn(1, 2, 2, 2, o)))
end
DIPD12 = ΔIPD12

ΔIL21 = Outfun("ΔIL") do o
    10 * log10(abs2(getsijmn(2, 1, 1, 1, o) / getsijmn(2, 1, 2, 2, o)))
end
DIL21 = ΔIL21

ΔIL12 = Outfun("ΔIL") do o
    10 * log10(abs2(getsijmn(1, 2, 1, 1, o) / getsijmn(1, 2, 2, 2, o)))
end
DIL12 = ΔIL12

ardb(i, j, n) =
    Outfun("AR$i$j($n)dB") do o
        jP = im * getsijmn(i, j, 1, n, o) / getsijmn(i, j, 2, n, o) # Modified Linear Pol. ratio
        Q = (1 - jP) / (1 + jP) # Circular polarization ratio
        absQ = abs(Q)
        absQ > 1 && (absQ = 1 / absQ)
        ardb = 20 * log10((1 + absQ) / (1 - absQ))
    end
AR11DB(n) = ardb(1, 1, n)
AR12DB(n) = ardb(1, 2, n)
AR21DB(n) = ardb(2, 1, n)
AR22DB(n) = ardb(2, 2, n)


FTHZ = Outfun("FTHZ") do o
    o.FGHz * 1e-3
end

FGHZ = Outfun("FGHZ") do o
    o.FGHz
end

FMHZ = Outfun("FMHZ") do o
    o.FGHz * 1e3
end

FKHZ = Outfun("FKHZ") do o
    o.FGHz * 1e6
end

FHZ = Outfun("FHZ") do o
    o.FGHz * 1e9
end

THETA = Outfun("THETA") do o
    get(o.steering, :θ, NaN)
end
Θ = THETA

PHI = Outfun("PHI") do o
    get(o.steering, :ϕ, NaN)
end
Φ = PHI

PSI1 = Outfun("PSI1") do o
    get(o.steering, :ψ₁, NaN)
end
Ψ₁ = PSI1

PSI2 = Outfun("PSI2") do o
    get(o.steering, :ψ₂, NaN)
end
Ψ₂ = PSI2


"""
    θϕ(o::Result) -> (θ, ϕ)

Return steering angles in degrees from `o`.  If `o` specifies `ψ₁` and `ψ₂` instead of angles,
then the latter are computed, using the Region 1 (input) periodicity and electrical parameters.
"""
function θϕ(o::Result)
    #=
    haskey(o.steering, :θ) && return (o.steering.θ, o.steering.ϕ)
    ψ₁, ψ₂ = o.steering
    units_per_meter = ustrip(Float64, o.unitsin, 1u"m")
    s₁, s₂ = [o.s₁in, o.s₂in] / units_per_meter
    β₁, β₂ = s₁s₂2β₁β₂(s₁,s₂)
    β₀₀ = (ψ₁ * β₁ + ψ₂ * β₂) / (2π)
    =#
    β₀₀² = o.β⃗₀₀ ⋅ o.β⃗₀₀
    β₀₀² == 0 && return (0.0, get(o.steering, :ϕ, 0.0))
    k² = (twopi * o.FGHz * 1e9 / c₀)^2 * real(o.ϵᵣin * o.μᵣin)
    β₀₀² > k² && error("Cut-off dominant mode")
    haskey(o.steering, :θ) && return (o.steering.θ, o.steering.ϕ)
    kz = √(k² - β₀₀²)  # for out-going wave vector in Layer 1
    θ = acosd(kz / sqrt(k²))
    ϕ = atand(o.β⃗₀₀[2], o.β⃗₀₀[1])
    return (θ, ϕ)
end

"""
    ĥv̂(θ, ϕ)

Compute Ludwig 3 unit vectors from spherical location vectors.
"""
function ĥv̂(θ, ϕ)
    st, ct = sincosd(θ)
    sp, cp = sincosd(ϕ)
    θ̂ = @SVector [ct * cp, ct * sp, -st]
    ϕ̂ = @SVector [-sp, cp, 0.0]
    ĥ = θ̂ * cp - ϕ̂ * sp
    v̂ = θ̂ * sp + ϕ̂ * cp
    ĥ, v̂
end


"""
    @outputs(args...)

Convert list of user output requests to a vector of functors that generate the requested
outputs when applied to a `Result` instance.  In the conversion process, replace
lower case letters with upper case.

### Examples

    julia> output = @outputs FGHz θ ϕ s11db(te,te) S11ang(Te,te)
    julia> output = @outputs FGHz theta phi s21db(R,H) ARdB21(H) ARdB11(v)
"""
macro outputs(args...)
    newargs = Any[]
    for (iarg, arg) in pairs(args)
        if arg isa Symbol
            push!(newargs, Symbol(uppercase(string(arg))))
        elseif arg isa Expr && arg.head == :call
            for (iarg2, arg2) in pairs(arg.args)
                arg2 isa Symbol && (arg.args[iarg2] = Symbol(uppercase(string(arg2))))
            end
            push!(newargs, arg)
        else
            error("Illegal @outputs construction")
        end
    end
    Tuple(eval(a) for a in newargs)
end


"""
    append_result_data(fname::AbstractString, gname::String, result::Result)

Append a `Result` instance to a result file for a particular frequency and pair of scan parameters.

## Arguments

- `fname`: The name of the result file to be appended to.
- `gname`: The unique `JLD2` group name to be used in the file for grouping the data
  associated with this frequency/scan case.
- `result`:  The `Result` data to be written to the file.
"""
function append_result_data(fname::AbstractString, gname::String, result::Result)
    jldopen(fname, "a") do fid
        group = JLD2.Group(fid, gname)
        group["result"] = result
    end
    return
end

append_result_data(::Base.DevNull, ::String, ::Result) = nothing

"""
    read_result_file(fname::AbstractString) --> Vector{Result}

Read a result file (in JLD2 format) and return a vector of results.
"""
function read_result_file(fname::AbstractString)::Vector{Result}
    dat = load(fname) # a Dict
    ks = collect(keys(dat))
    sort!(ks, by=x -> parse(Int, split(x, '/')[1]))
    Result[dat[k] for k in ks]
end


function extract_result_file(fname::AbstractString, ops::Tuple)
    results = read_result_file(fname)
    [o(r) for r in results, o in ops]
end

"""
    extract_result(r::Result, ops::NTuple{N,Outfun}) --> Row Matrix
    extract_result(r::AbstractVector{Result}, ops::NTuple{N,Outfun}) --> Matrix

Return a matrix of outputs extracted from a `Result` instance or vector.  `ops` is a
`NTuple` as returned by the `@outputs` macro.

### Example
    results = analyze(...)
    ops = @outputs FGHz s11dB(h,h) s11ang(h,h)
    data = extract_result(results, ops)
    # or data = extract_result(results[1], ops) # returns a single row
"""
function extract_result(results::AbstractVector{Result}, ops::Tuple)
    [o(r) for r in results, o in ops]
end

function extract_result(results::Result, ops::Tuple)
    permutedims([o(results) for o in ops])
end

"""
    extract_result(fname::AbstractString, ops::Tuple) --> Matrix

Return a matrix of outputs extracted from a results file.  `ops` is a
Tuple returned by the `@outputs` macro.

### Example
    ops = @outputs FGHz S11DB(H,H) S11ANG(H,H)
    data = extract_result("pssfss.res", ops)
"""
extract_result(fname::AbstractString, ops::Tuple) = extract_result_file(fname, ops)


"""
    _check_results_for_tep!(results::Vector{Result})

Verify that the input vector of `Result` objects is suitable for conversion to a `TEPperiodic` object.

The following checks are performed:
1. Ensure that incidence angles rather than incremental phasings are used.
1. If more than one ϕ value is used, ensure that all ϕ values in the range `0:Δϕ:(360-Δϕ)` are present.
1. Ensure that the complete Cartesian product of angles and frequencies is present.

Finally, the vector is sorted into the appropriate order for storage into a `TEPperiodic` object.
"""
function _check_results_for_tep!(results::Vector{Result})
    kys = keys(results[1].steering)
    kys != (:θ, :ϕ) && error("results do not use θ and ϕ for steering")
    # Rearrange into proper order for storage in a TEPperiodic:
    sort!(results, by = r -> (r.FGHz, r.steering.ϕ, r.steering.θ))
    freqs = unique!([r.FGHz for r in results])
    thetas = unique!([r.steering.θ for r in results])
    phis = unique!([r.steering.ϕ for r in results])
    nf = length(freqs)
    nt = length(thetas)
    np = length(phis)
    dphi = isone(np) ? 360 : phis[2] - phis[1]
    phi_range = 0:dphi:(360-dphi)
    phis == phi_range || error("ϕ values must be equivalent to 0:Δϕ:(360-Δϕ)")
    dtheta = thetas[2] - thetas[1]
    theta_range = 0:dtheta:last(thetas)
    thetas == theta_range || error("θ values must be equivalent to 0:Δθ:θmax")
    freqs_vec = [f*u"GHz" for f in freqs]
    # Verify that full cartesian product of angles and frequencies is present:
    i = 0
    for f in freqs, p in phis, t in thetas
        i += 1
        r = results[i]
        if (r.steering.ϕ != p) || (r.steering.θ != t) || (r.FGHz != f)
            error("Missing case (FGHz, ϕ, θ) = ($f, $p, $t) in input vector")
        end
    end
    return (theta_range, phi_range, freqs_vec)
end # function

"""
    res2tep(results::Vector{Result}; name="tep", class="res2tep") -> t::TEPperiodic
    res2tep(resultfile::AbstractString; name="tep", class="res2tep") -> t::TEPperiodic
    res2tep(results::Vector{Result}, tepfile::AbstractString; name="tep", class="res2tep") -> t::TEPperiodic
    res2tep(resultfile::AbstractString, tepfile::AbstractString; name="tep", class="res2tep") -> t::TEPperiodic

Convert a vector of `Result` elements into a `TEPperiodic` object, as defined in the
[TicraUtilities](https://github.com/simonp0420/TicraUtilities.jl) package.  If positional argument
`tepfile` is provided, the `TEPperiodic` object will be saved to this file name as a TICRA-compatible
TEP (tabulated electrical properties) file. If the first positional argument is an `AbstractString`, it is
assumed to be the name of a PSSFSS results file, from which the vector of results will be read.
The keyword arguments are used to provide values for the same-named fields in the TEP structure.

## Requirements for TEP File Compatibility
Because a TEP file contains all of the information of the full 4×4 scattering matrix computed by PSSFSS, there are no
limitations on the type of unit cell geometry that can be used for creating TEP files.

TEP files use the concept of "front" and "rear" incidence.  When converting a PSSFSS analysis result to TEP format, Region 1
(the first layer in the `strata` vector) is taken as the "front" incidence region, and Region `n` (the last layer) is
taken to be the "rear" region.  Both of these layers should have zero width and assume vacuum electrical parameters.  I.e.,
they should be specified as `Layer()` in the `strata` stackup.

`results` (or `resultfile`) must contain the results of a PSSFSS analysis sweep over θ and ϕ (and possibly frequency) such that
1. Incidence angles θ and ϕ rather than incremental phasings ψ₁ and ψ₂ were used.
1. If more than one ϕ value is used, then all ϕ values in the range `0:Δϕ:(360-Δϕ)` must be present.
1. The entire 3-dimensional Cartesian product of all angles and frequencies must be present.
"""
function res2tep(results::Vector{Result}, tepfile::AbstractString; name="tep", class="res2tep")
    t = res2tep(results; name, class)
    TicraUtilities.write_tepfile(tepfile, t)
    return t
end

function res2tep(resultfile::AbstractString, tepfile::AbstractString; name="tep", class="res2tep")
    results = read_result_file(resultfile)
    t = res2tep(results; name, class)
    TicraUtilities.write_tepfile(tepfile, t)
    return t
end

function res2tep(resultfile::AbstractString; name="tep", class="res2tep")
    results = read_result_file(resultfile)
    return res2tep(results; name, class)
end

function res2tep(results::Vector{Result}; name="tep", class="res2tep")
    theta, phi, freqs = _check_results_for_tep!(results)
    mff = @SMatrix [1 -1; -1 1]
    nt, np, nf = length.((theta, phi, freqs))
    sff, sfr, srf, srr = (zeros(ComplexF64, (2,2,nt,np,nf)) for _ in 1:4)
    i = 0
    for ifr in 1:nf, ip in 1:np, it in 1:nt
        i += 1
        gsm = results[i].gsm
        sff[:,:,it,ip,ifr] .= mff .* gsm[1,1]
        sfr[:,:,it,ip,ifr] .= mff .* gsm[1,2]
        srf[:,:,it,ip,ifr] .= mff .* gsm[2,1]
        srr[:,:,it,ip,ifr] .= mff .* gsm[2,2]
    end

    tep = TicraUtilities.TEPperiodic(; name, class, theta, phi, freqs, sff, sfr, srf, srr)
    return tep
end # function


"""
    _prepare_results_for_fresnel(results::Vector{Result}) -> Vector{Result} -> newresults::Vector{Result}

Verify that the input vector of `Result` objects is suitable for creation of an HFSS-compatible Fresnel table.

The following checks are performed:
1. Ensure that incidence angles rather than incremental phasings are used.
2. Ensure that the θ angles begin at 0 and are uniformly spaced up to the maximum θ value present.
3. Ensure that the increment in θ values divides evenly into 90.
4. Ensure that if multiple frequencies are present, then they have a uniform spacing.

The output vector is modified in the following ways:
1. Results for all ϕ angles other than the minimum ϕ in absolute value are discarded.
2. The remaining results are sorted into order of increasing θ, then increasing frequency.
3. The vector is sorted into the appropriate order for creating a Fresnel table: frequency
   varies most rapidly, then θ.
4. If the input θ values do not extend all the way to 90°, then the results for the maximum supplied θ are copied and appended
   a sufficient number of times to simulate data all the way to 90°.
"""
function _prepare_results_for_fresnel(results::Vector{Result})
    kys = keys(results[1].steering)
    kys != (:θ, :ϕ) && error("results do not use θ and ϕ for steering")

    ϕmin = minimum(abs, (r.steering.ϕ for r in results))
    results = filter(r -> r.steering.ϕ == ϕmin, results)

    sort!(results, by = r -> (r.steering.θ, r.FGHz))

    freqs = unique!([r.FGHz for r in results])
    thetas = unique!([r.steering.θ for r in results])
    iszero(first(thetas)) || error("First θ value is not equal to zero")
    nf = length(freqs)
    if nf > 1
        df = freqs[2] - freqs[1]
        frange = first(freqs):df:last(freqs)
        freqs == frange || error("frequencies are not uniformly spaced")
    end
    nt = length(thetas)
    nt ≤ 1 && error("Too few θ value")
    dtheta = thetas[2] - thetas[1]
    isinteger(90 / dtheta) || error("θ increment does not divide evenly into 90")
    theta_range = 0:dtheta:last(thetas)
    thetas == theta_range || error("θ values are not uniformly spaced")

     # Verify that full cartesian product of angles and frequencies is present:
    i = 0
    for t in thetas, f in freqs
        i += 1
        r = results[i]
        if (r.steering.θ != t) || (r.FGHz != f)
            error("Missing case (FGHz, ϕ, θ) = ($f, $phimin, $t) in input vector")
        end
    end

    # Fill in any missing theta values
    if last(thetas) < 90
        newthetas = float((dtheta + last(thetas)):dtheta:90)
        irng = (length(results) - nf + 1):length(results) # Indices in results for all frequencies of last θ input value
        for newtheta in newthetas, i in irng
            r = results[i]
            steer = (θ = newtheta, ϕ = r.steering.ϕ)
            rnew = Result(r.gsm, steer, r.β⃗₀₀, r.FGHz, r.ϵᵣin, r.μᵣin, r.β₁in, r.β₂in, r.ϵᵣout, r.μᵣout, r.β₁out, r.β₂out)
            push!(results, rnew)
        end
    end

    return results, union(thetas, newthetas), freqs
end # function



"""
    res2fresnel(results::Vector{Result}, fresnelfile::AbstractString)
    res2fresnel(resultfile::AbstractString, fresnelfile::AbstractString)

Create an HFSS-compatible "Fresnel table" file from `results`, the vector of `Result` objects returned by
the `analyze` function.  If the first positional argument is an `AbstractString`, it is
assumed to be the name of a PSSFSS results file, from which the vector of results will be read.

Since Fresnel tables contain data for only a single ϕ value, if the input `result` vector contains data for multiple
ϕ values, only the value with minimum magnitude will be used.

Fresnel tables may be formatted to contain only reflection coefficients (for a so-called "opaque" structure), or they
may contain both reflection and transmission coefficients (a "non-opaque" structure).
An opaque structure is one for which the s21 partition of the generalized scattering matrix is identically zero
for all frequencies and scan angles.  The correct format to be written will be selected automatically by `res2fresnel`.

## Requirements for Fresnel Table Compatibility
The data in `results` must satisfy the following requirements:
1. Incidence angles rather than incremental phasings must be used.
2. θ angles must begin at 0 and be uniformly spaced up to the maximum θ value present.
3. The increment in θ values must divide evenly into 90.
4. If multiple frequencies are present, then they must have a uniform spacing.

A Fresnel table must contain θ values equally spaced between 0 and 90, inclusive.
If the `results` vector provided as input does not contain θ values all the way to 90, then the scattering matrix values
corresponding to the maximum provided θ value will be copied into the remaining angular "slots" as necessary to provide
a complete Fresnel table.

There are some limitations on the type of unit cell geometry that should be used for creating Fresnel tables.  First, a Fresnel
table contains data for only a single ϕ value.  This means that the geometry being analyzed must be such that the scattering
matrix of the structure is essentially independent of ϕ.  As a counterexample, a strip grid is not a suitable structure, since
its scattering properties are strongly dependent on ϕ.  Second, a Fresnel table records only co-polarized
(TE → TE and TM → TM) transmission and reflection coefficients.  This means that the structure being analyzed must not
generate cross-polarized (TE → TM or TM → TE) transmission or reflection coefficients of significant amplitude.

Fresnel tables consider only incidence from a single "front" region. When creating the Fresnel table, the front region is taken
to be Region 1 of the PSSFSS model (i.e. the first layer present in the PSSFSS `strata` vector).

### Additional Requirements for Non-Opaque Structures
When used in an HFSS SBR+ model, the scattering properties read from the Fresnel table are applied to a zero-thickness surface,
so that the transmitted ray is launched from the same "hit" point of the surface that was encountered by the incident
ray. Because of this, the phase reference plane for both input and output ports of the PSSFSS model should be located
at this front surface (i.e. the first interface plane in the `strata` vector).  This is accomplished by specifying zero
width for the first `Layer` object (i.e. using `Layer()` for the first layer), and then specifying the final layer's width
to be the negative of the sum of all the other layer widths in the `strata` vector. The negative width value shifts
the output port reference plane to coincide with that of the input port.  As an example:
```julia
strata = [Layer(), Layer(width=2mm, ϵᵣ=2.2) Layer(width=3.3mm, ϵᵣ=3.0), Layer(width=2mm, ϵᵣ=2.2), Layer(width=-7.3mm)]
```
"""
function res2fresnel(results::Vector{Result}, fresnelfile::AbstractString)

    results, thetas, freqs = _prepare_results_for_fresnel(results)
    nt = length(thetas)
    nf = length(freqs)
    ronly = all(r -> iszero(r.gsm[2,1]), results) # reflection only
    open(fresnelfile, "w") do fid
        date, clock = split(string(now()), 'T')
        if ronly
            println(fid, "# HFSS-compatible Fresnel reflection table created by PSSFSS")
            println(fid, "# Created on ", date, " at ", clock)
            println(fid, "ReflTab1e")
        else
            println(fid, "# HFSS-compatible Fresnel reflection/transmission table created by PSSFSS")
            println(fid, "# Created on ", date, " at ", clock)
            println(fid, "RTTable")
        end

        println(fid, "# <num_theta_step> = <number_of_points> - 1")
        println(fid,  nt - 1)

        if nf == 1
            println(fid, "# Mono freq. table for ", only(freqs), " GHz")
            println(fid, "MonoFreq")
            println(fid, "# Data section follows.")
        else
            println(fid, "# MultiFreq <freq_start_ghz> <freq_stop_ghz> <num_freq_steps>")
            println(fid, "MultiFreq ", first(freqs), " ", last(freqs), " ", nf-1)
            println(fid, "# Data section follows. Frequency loops within theta")
        end

        if ronly
            println(fid, "#<rte_rl> <rte_im> <rtm_rl> <rtm_im>")
        else
            println(fid, "#<rte_rl> <rte_im> <rtm_rl> <rtm_im> <tte_rl> <tte_im> <ttm_rl> <ttm_im>")
        end

        for res in results
            r = res.gsm[1,1]; rte = r[1,1]; rtm = -r[2,2]
            t = res.gsm[2,1]; tte = t[1,1]; ttm = t[2,2]
            @printf(fid, "%8.5f %8.5f %8.5f %8.5f", real(rte), imag(rte), real(rtm), imag(rtm))
            if ronly
                println(fid)
            else
                @printf(fid, " %8.5f %8.5f %8.5f %8.5f\n", real(tte), imag(tte), real(ttm), imag(ttm))
            end
        end
    end
end # function

function res2fresnel(resultfile::AbstractString, fresnelfile::AbstractString)
    results = read_result_file(resultfile)
    res2fresnel(results, fresnelfile)
end


end # module