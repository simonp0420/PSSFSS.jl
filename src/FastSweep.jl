module FastSweep

using OffsetArrays: Origin
import LinearAlgebra: norm

@inline zerobased(x) = first(eachindex(x)) == 0

"""
    (Sinterp, errest) = interp_path2(x0, x0j, Sj0)

Rational function interpolation using a Path II Neville lozenge, as defined in the reference.

## Positional Input Arguments

- `x0`: Independent variable at which the rational function is to be evaluated.
- `x0j`: Vector of independent variable sample points. `type(x0)` must equal `eltype(x0j)`.
- `Sj0`: Vector of function values evaluated at `x0j`. `length(x0j)` must equal `length(Sj0)`.

## Optional Keyword Arguments

- `store1`, `store2`, `store3`: Optional storage vectors, which will be mutated.  If not supplied,
  they will be allocated on each call to this function. They should be of the same element type and
  length as `Sj0`.

## Return Values

 - `Sinterp`: The rational function interpolation at `x0`.
 - `errest`: An error estimate for `norm(Sinterp - S(x0))`.

## Reference
Ma, X., Wan, G. and Wan, W., 2012. "A Multi-Dimensional Adaptive Sampling Method for Analysis
and Design of Frequency Selective Surface with Arbitrary Element".
Progress In Electromagnetics Research B, 41, pp.213-230.
"""
function interp_path2(
    x0::T1,
    x0jin::AbstractVector{T1},
    Sj0::AbstractVector{T2};
    store1::AbstractVector{T2} = zeros(T2, length(Sj0)),
    store2::AbstractVector{T2} = zeros(T2, length(Sj0)),
    store3::AbstractVector{T2} = zeros(T2, length(Sj0)),
) where {T1<:Real,T2}

    # trivial case where x0 is one of the sample points
    for (x,S) in zip(x0jin,Sj0)
        iszero(x-x0) && (return (S, abs(zero(T2))))
    end

    # Input checking
    K1 = length(x0jin)
    K = K1 - 1
    K1 == length(Sj0) || error("length mismatch for x0j and Sj0")
    K1 ≤ length(store1) && K1 ≤ length(store2) && K1 ≤ length(store3) ||
        error("Storage vector lengths not all long enough")

    store1n = @view store1[begin:begin+K]
    store2n = @view store2[begin:begin+K]
    store3n = @view store3[begin:begin+K]

    for i in eachindex(store1n, store2n, store3n)
        store1n[i] = zero(T2)
        store2n[i] = zero(T2)
        store3n[i] = zero(T2)
    end

    # Use zero-based arrays for convenience:
    x0j = zerobased(x0jin) ? x0jin : Origin(0)(x0jin)
    Sjk = zerobased(store1n) ? store1n : Origin(0)(store1n)
    Sjkm1 = zerobased(store2n) ? store2n : Origin(0)(store2n)
    Sjkm2 = zerobased(store3n) ? store3n : Origin(0)(store3n)


    @inbounds for j = 0:K
        Sjk[j] = Sj0[begin+j] # Works for zero or one-based Sj0
    end
    @inbounds for k = 1:K
        Sjkm2 .= Sjkm1
        Sjkm1 .= Sjk
        @inbounds for j = 0:K-k
            num1 = x0 - x0j[j]
            den1 = Sjkm1[j+1] - Sjkm2[j+1]
            num2 = x0j[j+k] - x0
            den2 = Sjkm1[j] - Sjkm2[j+1]
            bigden = fixbigden(num1 * den2 + num2 * den1)
            Sjk[j] = Sjkm2[j+1]
            Sjk[j] += (x0j[j+k] - x0j[j]) * den1 .* den2 ./ bigden
        end
    end
    return (Sjk[0], norm(Sjk[0] - Sjkm1[0]))
end

fixbigden(x::T) where {T <: Number} = ifelse(iszero(x), one(T), x)

function fixbigden(x::AbstractArray{T}) where {T}
    for i in eachindex(x)
        x[i] = ifelse(iszero(x[i]), one(T), x[i])
    end
    return x
end


"""
    interpolate_band(f::F, x::AbstractVector; max_err_lim_db=-80, nrepeat=3, xlabel="", showprogress=false) where {F<:Function}

Rational function interpolation of the function `f` evaluated at each element of `x`.

`f` is a function that returns either a number, vector, or matrix such that `LinearAlgebra.norm(f[x[i]])`
 and `zero(f[x[i]])` are defined for any `i ∈ eachindex(x)`.
`max_err_lim_db` is the maximum allowed estimated error in dB
for any of the interpolated points, and `nrepeat` is the number of times this error criterion must be
consecutively met during the interpolation procedure before it has been considered to be satisfied.  Note
that the default values are very strict.  The return value is a vector of the same length as `x` containing
the interpolated function values.
"""
function interpolate_band(
    f::F,
    x::AbstractVector{T};
    showprogress = false,
    xlabel = "",
    max_err_lim_db = -80,
    nrepeat = 3) where {F<:Function, T<:Real}

    Base.require_one_based_indexing(x)
    len = length(x)
    if len < 6
        knots = collect(eachindex(x)) # Don't bother interpolating
    else
        knots = round.(Int, collect(LinRange(firstindex(x), lastindex(x), 5)))
    end

    crclear = "\r\u1b[K" # carriage return and clear rest of line

    showprogress && println("")
    f1 = f(x[knots[1]])
    showprogress && print(crclear, "1 knot at ", x[knots[1]], " ", xlabel, ", maxerrdB = Inf")
    fknots = Array{typeof(f1), 1}(undef, length(knots))
    fknots[1] = f1
    for k in (1+firstindex(knots)):lastindex(knots)
        fknots[k] = f(x[knots[k]])
        showprogress && print(crclear, k, " knots. Added ", x[knots[k]], " ", xlabel, ", maxerrdB = Inf")
    end
    errs = ones(len)
    errs[knots] .= 0.0
    mapreduce(iszero, &, errs) && (return fknots)
    fvalues = Array{eltype(fknots), 1}(undef, len)
    fvalues[knots] .= fknots
    store1, store2, store3 = (similar(fvalues) for _ in 1:3)
    max_err_lim = 10.0^(max_err_lim_db / 20)
    max_err = 1.0
    repeats = nextknot = 0

    while repeats < nrepeat || max_err > max_err_lim
        if nextknot > 0
            fvalues[nextknot] = f(x[nextknot])
            errs[nextknot] = 0.0
            push!(knots, nextknot)
        end
        for i in eachindex(x, fvalues)
            i ∈ knots && continue
            xview = view(x, knots)
            fview = view(fvalues, knots)
            fvalues[i], errs[i] = interp_path2(x[i], xview, fview; store1, store2, store3)
        end
        (max_err, nextknot) = findmax(errs)
        if showprogress
            maxerrdB = round(20*log10(max_err), digits=2)
            print(crclear, length(knots), " knots. Added ", x[knots[end]], " ", xlabel, ", maxerrdB = ", maxerrdB)
        end
        iszero(max_err) && break
        max_err ≤ max_err_lim && (repeats += 1)
        max_err > max_err_lim && (repeats = 0)
    end

    showprogress && print("\u1b[1F") # Move to start of prev line
    return fvalues
end # function


end # Module
