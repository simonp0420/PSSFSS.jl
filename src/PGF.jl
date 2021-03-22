module PGF

export electric_modal_sum_funcs, magnetic_modal_sum_funcs, jksums

using LinearAlgebra: norm, ⋅, ×
using OffsetArrays
using ..Rings: Ring
using ..Sheets: MV2
using ..Layers: Layer
using FFTW: fft!
using ..Constants: tdigits

# Variables used by the spatial routines:
const jkringmax = 65 # Max. number of rings to sum over
  
# Variables used by the spectral routines:
const mmax_list = (32, 2048)


#=
"""
    ring(r::Integer)

Return vector of (m,n) pairs comprising the r'th summation ring.
"""
function ring(r::Integer)
    r == 0 && return [(0,0)]
    leftright = vec([(m,n) for m in (-r,r), n in -r:r])
    topbot = vec([(m,n) for m in 1-r:r-1, n in (-r,r)])
    return vcat(leftright,topbot)
end
=#




"""
    jksums(uρ₀₀, ψ₁, ψ₂, us₁, us₂, extract::Bool; convtest=1e-8) --> (jsum, ksum)
                                                                  
Compute the frequency-independent sums defined in Equations (7.32c) and (7.32d) of the 
theory documentation.

## Arguments:

- `uρ₀₀`: 2-vector containing the difference of the observation and source vectors, 
          multiplied by `u`, the smoothing parameter.
- `ψ₁`, `ψ₂`:  The unit cell incremental phase shifts (radians).
- `us₁`, `us₂`:  2-vectors containing unit cell direct lattice vectors multiplied 
                 by `u`, the smoothing parameter.
- `extract`: If true, directs this function to perform singularity extraction 
    in the jsum calculation.
- `convtest`: Relative convergence criterion.  If the latest loop contributions
              for both sums are smaller than `convtest` times the magnitude of the 
              latest respective sum values then the sums are assumed to have converged.

## Return values

- `jsum`: The complex sum appearing in the integral of Equations (7.22c) and (7.32c) of 
              the theory documentation.
- `ksum`: The complex sum appearing in the integral of Equation (7.22d) and (7.32d) of 
              the theory documentation.
"""
function jksums(uρ₀₀, ψ₁, ψ₂, us₁, us₂, extract, convtest=1e-8)
    rmin = 5 # min number rings to sum over
    small = 1e-5  # Test value for small-argument approximation

    arg = norm(uρ₀₀)
    rterm = exp(-arg)
    ksum = complex(rterm)
    if extract 
        if arg ≤ small
            jsum = complex((-arg/6 + 0.5)*arg - 1) # Small argument extraction formula
        else
            jsum = complex((rterm - 1) / arg)  # Large argument extraction formula
        end
    else
      jsum = complex(rterm / arg)  # No extraction formula
    end
    
    # Begin loop over summation lattice rings.
    rsave = 0; conv = 0.0; converged = false # Establish scope outside loop
    for r in 1:jkringmax  
        rsave = r
        jring = kring = zero(ComplexF64) # Initialize ring sums
        for (m,n) in Ring(r)
            arg = norm(uρ₀₀ - (m*us₁ + n*us₂))
            term = exp(-complex(arg, m*ψ₁ + n*ψ₂))
            jring += term / arg
            kring += term
        end
        jsum += jring
        ksum += kring
        # Test for convergence if we're far enough along:
        if r > rmin
            conv = max(abs(jring/jsum), abs(kring/ksum))
            converged = conv < convtest
            converged && break
        end
    end
    converged || @warn "Exceeded maximum number of loops." rsave conv maxlog=10
    return (jsum, ksum)
end



"""
    mysqrt(x)

Same as `sqrt` unless `sqrt(x)` is pure negative imaginary in which
case it returns `-sqrt(x)` (i.e., positive pure imaginary).
"""
mysqrt(x) = sqrt(x)
function mysqrt(z::Complex)
    ans = sqrt(z)
    return  real(ans) == 0 && imag(ans) < 0 ? -ans : ans
end





"""
    c3_calc(k0, u , μ₁, ϵ₁, μ₂, ϵ₂)

Compute the magnetic vector potential expansion coefficient `c3`, defined
in Equation (4.30) of the theory documentation. The units of `c3` are the 
square of the units of `u` (or `k0`).

## Arguments:

- `k0` Free-space wavenumber.
- `u` Smoothing factor.  `k0` and `u` can be of any 
"""
function c3_calc(k0, u , μ₁, ϵ₁, μ₂, ϵ₂)
    w1sq = k0*k0 * μ₁*ϵ₁ + u*u    # Eq. (4.25)
    w2sq = k0*k0 * μ₂*ϵ₂ + u*u    
    c3 = (μ₁*w2sq + μ₂*w1sq) / (2*(μ₁ + μ₂)) 
end
  

"""
    d3_calc(k0, u , μ₁, ϵ₁, μ₂, ϵ₂)

Compute the electric scalar potential expansion coefficient `d3`, defined
in Equation (4.33) of the theory documentation. The units of `d3` are the 
square of the units of `u` (or `k0`).

## Arguments:

- `k0` Free-space wavenumber.
- `u` Smoothing factor.  `k0` and `u` can be of any 
"""
function d3_calc(k0, u , μ₁, ϵ₁, μ₂, ϵ₂)
    w1sq = k0*k0 * μ₁*ϵ₁ + u*u  # Eq. (4.25)
    w2sq = k0*k0 * μ₂*ϵ₂ + u*u
    d3_num = μ₁*(w2sq*(2*ϵ₁ + ϵ₂) - w1sq*ϵ₁) + μ₂*(w1sq*(2*ϵ₂ + ϵ₁) - w2sq*ϵ₂)       
    d3 = d3_num / (2*(μ₁ + μ₂)*(ϵ₁ + ϵ₂))  
end

  
  
"""
    electric_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, convtest=5e-12) --> (Σm1_func, Σm2_func)
   
Return a pair of functions that efficiently compute the modal series for the magnetic vector 
potential and electric scalar potential as defined in Eqs. (5.19) of the theory documentation.

## Arguments:

- `k0`: Free-space wavenumber (1/meter).
- `u`:  Smoothing parameter (1/meter).
- `ψ₁`, `ψ₂`:  Unit cell incremental phase shifts (radians).
- `layers`  An `AbstractVector` of element type `Layer` containing the layer 
       parameters for the cascade structure.  Note that the first
       and last layer's thicknesses are not accounted for in this
       function.  They are assumed to be semi-infinite.
- `s`  Interface number (within layers) at which the FSS or PSS sheet is located.
- `β₁`, `β₂`:   2-vectors containing the reciprocal lattice basis 
              vectors.  Units are (1/meters).
- `β₀₀`    2-vector containing the principal (i.e. with index (0,0)) Floquet 
           vector transverse wavenumber which incorporates the intrinsic phase 
           shifts.  Units are (1/meters).
- `convtest` Relative convergence criterion. 

##  Return Values
    
A pair of functions that evaulate the two series defined in Eq. (5.19) of the
theory documentation.  Each function takes a single argument `ρdif`, a 2-vector
containing the difference of the observation and source point position vectors.

"""
function electric_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers::AbstractVector{Layer},
                                             s, β₁, β₂, β₀₀, convtest=5e-12)

    t1 = time_ns()
    nl = length(layers) # Number of layers.
    nl < 2 && error("Too few layers")
    k0sq = k0 * k0
    β²min = 1e-10*k0sq
    # Compute quantities defined in Eq. (5.20):
    μ̃ = 2 * layers[s].μᵣ * layers[s+1].μᵣ / (layers[s].μᵣ + layers[s+1].μᵣ) # Normalized to μ0.
    ϵ̄ = (layers[s].ϵᵣ + layers[s+1].ϵᵣ) / 2 # Normalized to ϵ0.
    area = 4π^2 / norm(β₁ × β₂) # unit cell area (m^2):
    # Obtain the appropriate Green's function expansion coefficients:
    c3 = c3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)
    d3 = d3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)

    m = mmax_list[2] ÷ 2
    table1 = OffsetArray(zeros(ComplexF64, 2m+1, 2m+1), -m:m, -m:m) 
    table2 = OffsetArray(zeros(ComplexF64, 2m+1, 2m+1), -m:m, -m:m)

    converged = false
    convrepeat = 40 # number of consecutive rings for which convergence must occur
    convlist = OffsetArray([false for i in 0:mmax_list[end]÷2], 0:mmax_list[end]÷2)
    mmax = mmax_list[1] ÷ 2
    mmax_old = -2
    test1 = test2 = 0.0 # Establish scope
    first = true
    while mmax < mmax_list[end] # Convergence loop
        if first
            mmax = 0
            first = false
        else
            mmax = nextprod([2,3,5], 1 + mmax)
            while 0 ≠ rem(mmax,2)
                mmax = nextprod([2,3,5], 1 + mmax)
            end
        end
        mmaxo2 = mmax÷2
        mmax_oldo2 = mmax_old÷2
        # Fill the tables:
        Threads.@threads for r in (mmax_oldo2+1):mmaxo2
            ringsum1 = zero(eltype(table1))
            ringsum2 = zero(eltype(table2))
            for (m,n) in Ring(r)
                βmn = β₀₀ + m*β₁ + n*β₂   # Modal transverse wave vector
                β² = βmn ⋅ βmn # magnitude squared
                β² = max(β², β²min)  # Avoid singularity
                κmn² = β² + u*u #  Eq. (4.24)
                κmn = sqrt(κmn²) 
                # Compute region 1 mode parameters:
                γ = mysqrt(β² - k0sq * layers[1].ϵᵣ * layers[1].μᵣ)
                # Calculate TE modal impedance of Region 1, divided by (jωμ₀):
                ZleftTE = layers[1].μᵣ / γ
                # Calculate TM modal impedance of Region 1, multiplied by (jωϵ₀):
                ZleftTM = γ / layers[1].ϵᵣ
                # Step left-looking impedances up to junction s using (5.14):
                for i in 2:s
                    γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                    Z0TE = layers[i].μᵣ / γ   # TE modal impedance (normalized)
                    Z0TM = γ / layers[i].ϵᵣ  # TM modal impedance (normalized)
                    tanhi = tanh(layers[i].width * γ)
                    # Eq. (5.14b):
                    ZleftTE = Z0TE * (ZleftTE + Z0TE*tanhi) / (Z0TE + ZleftTE*tanhi)
                    ZleftTM = Z0TM * (ZleftTM + Z0TM*tanhi) / (Z0TM + ZleftTM*tanhi)
                end
                # Compute region nl mode parameters:
                γ = mysqrt(β² - k0sq * layers[end].ϵᵣ * layers[end].μᵣ)
                ZrightTE = layers[end].μᵣ / γ # divided by (jωμ₀)
                ZrightTM = γ / layers[end].ϵᵣ # multiplied by (jωϵ₀)
                # Step right-looking impedances down to junction s using (5.14d):
                for i in nl-1:-1:s+1
                    γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                    Z0TE = layers[i].μᵣ / γ
                    Z0TM = γ / layers[i].ϵᵣ  
                    tanhi = tanh(layers[i].width * γ) 
                    ZrightTE = Z0TE * (ZrightTE + Z0TE*tanhi) / (Z0TE + ZrightTE*tanhi)
                    ZrightTM = Z0TM * (ZrightTM + Z0TM*tanhi) / (Z0TM + ZrightTM*tanhi)
                end
                # Compute (normalized) TLGF's using Eq. (5.13):
                ViTE = ZleftTE * ZrightTE / (ZleftTE + ZrightTE) # divided by jωμ₀
                ViTM = ZleftTM * ZrightTM / (ZleftTM + ZrightTM) # multiplied by jωϵ₀
                # Compute summands (apart from phase factor and 1/(2A) factor):
                ViTE = 2 * ViTE / μ̃ # 1st quantity in square brackets in (5.19a)
                table1[m,n] = ViTE - (1 + c3/κmn²) / κmn  # Eq. (5.19a)
                ringsum1 += table1[m,n]
                ViTM = 2 * ϵ̄ * ViTM # 1st quantity in square brackets in (5.19b)
                table2[m,n] = (ViTM + k0sq*ϵ̄*μ̃*ViTE) / β² - (1 + d3/κmn²) / κmn  # Eq. (5.19b)
                ringsum2 += table2[m,n]
            end
            # Check for convergence of this ring
            test1 = abs(ringsum1/table1[0,0])
            test2 = abs(ringsum2/table2[0,0])
            convlist[r] = test1 < convtest && test2 < convtest
        end
        
        #  Check for convergence
        if mmaxo2 ≥ convrepeat && all(@view convlist[(mmaxo2 - convrepeat + 1):mmaxo2])
            converged = true
            break
        else
            mmax_old = mmax
        end
    end

    !converged && @warn "Inadequate Convergence" test1 test2 convtest mmax maxlog=5

    # Create proper sized storage arrays for FFT routine:
    mmaxo2 = mmax÷2
    table1t = table1[-mmaxo2:mmaxo2-1,-mmaxo2:mmaxo2-1]
    table2t = table2[-mmaxo2:mmaxo2-1,-mmaxo2:mmaxo2-1]
    fft!(table1t)
    fft!(table2t) 
    # Adjust phase according to Equation (5.32).  Also, include factor of 1/(2*area)
    for q in 0:mmax-1
        qterm = q * (π - ψ₂/mmax)
        for p in 0:mmax-1
            pterm = p * (π - ψ₁/mmax)
            cfact = cis(pterm+qterm) / (2*area)
            table1t[p+1,q+1] *= cfact
            table2t[p+1,q+1] *= cfact
        end
    end
    # Create proper sized interpolation array---Note that we add an extra row
    # and extra column at both the beginning and end of each table to allow 
    # for extra points needed in the interpolation scheme.
    table1 = OffsetArray(zeros(ComplexF64, mmax+2, mmax+2), -1:mmax, -1:mmax)
    table1[0:mmax-1, 0:mmax-1] = table1t
    table2 = OffsetArray(zeros(ComplexF64, mmax+2, mmax+2), -1:mmax, -1:mmax)
    table2[0:mmax-1, 0:mmax-1] = table2t
    # Add extra row and column to cover all the way to ξ=1 and η=1:
    table1[mmax, 0:mmax-1] = cis(-ψ₁) * table1[0, 0:mmax-1]
    table1[0:mmax-1, mmax] = cis(-ψ₂) * table1[0:mmax-1, 0]
    table1[mmax, mmax] = cis(-(ψ₁+ψ₂)) * table1[0,0]
    table2[mmax, 0:mmax-1] = cis(-ψ₁) * table2[0, 0:mmax-1]
    table2[0:mmax-1, mmax] = cis(-ψ₂) * table2[0:mmax-1, 0]
    table2[mmax, mmax] = cis(-(ψ₁+ψ₂)) * table2[0,0]
    # Add extra row and column to cover all the way to ξ=-1/mmax and η=-1/mmax:
    table1[-1, 0:mmax] = cis(ψ₁) * table1[mmax-1, 0:mmax]
    table1[-1:mmax-1, -1] = cis(ψ₂) * table1[-1:mmax-1, mmax-1]
    table2[-1, 0:mmax] = cis(ψ₁) * table2[mmax-1, 0:mmax]
    table2[-1:mmax-1, -1] = cis(ψ₂) * table2[-1:mmax-1, mmax-1]
    # Use meaningful names:
    Σm1_func = make_Σm_func(table1, β₁, β₂, ψ₁, ψ₂)
    Σm2_func = make_Σm_func(table2, β₁, β₂, ψ₁, ψ₂)
    t2 = time_ns()
    tsec = round((t2-t1)/1e9; digits=tdigits)
    @info "      $tsec seconds to compute $mmax × $mmax electric modal tables"
    return (Σm1_func, Σm2_func)
end    



"""
    make_Σm_func(table::AbstractArray, β₁::MV2, β₂::MV2, ψ₁::Real, ψ₂::Real) -> Σm_func

Return a one-argument function `Σm` that evaluates one of the four modal series defined in
Equations (5.19) and (5.26) of the theory documentation. The single argument to `Σm` is a 2-vector
containing the difference of the observation and source positions. The evaluation is performed by
using a 6-point interpolation into a precomputed table.

## Arguments:

- `table`:  An `OffsetArray` generated by the function `electric_modal_sum_funcs` or `magnetic_modal_sum_funcs`, 
            with both axes consisting of `-1:mmax`.
- `β₁`, `β₂`:   2-vectors containing the reciprocal lattice basis 
              vectors.  Units are (1/meters).
- `ψ₁`, `ψ₂`: Incremental unit cell phase shifts (radians).

## Return value:

- `Σm_func`: A function that takes a single argument `ρdif`, a 2-vector containing the difference between
             observation and source points, and returns the value of the modal sum via interpolation in the 
             precomputed table.

"""
function make_Σm_func(table::AbstractArray, β₁::MV2, β₂::MV2, ψ₁::Real, ψ₂::Real)
    axes(table,1) == axes(table,2) || error("Non-square table")
    mmax = maximum(axes(table,1))
    twopi = 2π
    function Σm_func(ρdif)::eltype(table)
        let table=table, β₁=β₁, β₂=β₂, ψ₁=ψ₁, ψ₂=ψ₂, mmax=mmax, twopi=twopi
            # Obtain ξ₁ and ξ₂ using Equation (2.10) or (5.31) of the theory docs:
            ξ₁_orig = (β₁ ⋅ ρdif) / twopi
            ξ₁_orig = abs(ξ₁_orig) < 1e-8 ? 0.0 : ξ₁_orig
            ξ₂_orig = (β₂ ⋅ ρdif) / twopi
            ξ₂_orig = abs(ξ₂_orig) < 1e-8 ? 0.0 : ξ₂_orig
            # Adjust values to place in the interval [0,1):
            ξ₁ = mod(ξ₁_orig, 1)
            ξ₂ = mod(ξ₂_orig, 1)
            # Calculate number of unit cell shifts (see Eq. (5.30)):
            mshift = round(Int, ξ₁_orig - ξ₁)
            nshift = round(Int, ξ₂_orig - ξ₂)
            # Determine m,n,p,q so that the point (ξ₁,ξ₂) is in the square
            # bounded by (m/mmax,n/mmax) and ((m+1)/mmax,(n+1)/mmax2) with p and q
            # the fractional distances along the square:
            m = trunc(Int, mmax * ξ₁)
            p = mmax * ξ₁ - m 
            n = trunc(Int, mmax * ξ₂)
            q = mmax * ξ₂ - n 
            pp = 1 - p  # p' in interpolation formula
            qp = 1 - q  # q' in interpolation formula
            # Perform the interpolation using formula 25.2.67 of AMS-55:
            Σm = q*(q-1)/2 * table[m,n-1] +
                 p*(p-1)/2 * table[m-1,n] + 
                (1 + p*q - p*p - q*q) * table[m,n] +
                p*(p - 2*q + 1)/2 * table[m+1,n] +
                q*(q - 2*p + 1)/2 * table[m,n+1] +
                p*q * table[m+1,n+1] 
            # Add any phase shift due to range:
            if mshift ≠ 0 || nshift ≠ 0
                Σm *= cis(-(mshift*ψ₁ + nshift*ψ₂)) # Eq. (5.30)
            end
            return Σm
        end # let block
    end # closure
    return Σm_func
end




"""
    magnetic_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀,convtest=5e-12) --> (Σpm1_func, Σpm2_func)
   
Return a pair of functions that efficiently compute the modal series for the electric vector 
potential and magnetic scalar potential as defined in Eqs. (5.26) of the theory documentation.

## Arguments:

- `k0`: Free-space wavenumber (1/meter).
- `u`:  Smoothing parameter (1/meter).
- `ψ₁`, `ψ₂`:  Unit cell incremental phase shifts (radians).
- `layers`  An `AbstractVector` of element type `Layer` containing the layer 
       parameters for the cascade structure.  Note that the first
       and last layer's thicknesses are not accounted for in this
       function.  They are assumed to be semi-infinite.
- `s`  Interface number (within layers) at which the FSS or PSS sheet is located.
- `β₁`, `β₂`:   2-vectors containing the reciprocal lattice basis 
              vectors.  Units are (1/meters).
- `β₀₀`    2-vector containing the principal (i.e. with index (0,0)) Floquet 
           vector transverse wavenumber which incorporates the intrinsic phase 
           shifts.  Units are (1/meters).
- `convtest` Relative convergence criterion. 

##  Return Values
    
A pair of functions that evaulate the two series defined in Eq. (5.26) of the
theory documentation.  Each function takes a single argument `ρdif`, a 2-vector
containing the difference of the observation and source point position vectors.

"""
function magnetic_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers::AbstractVector{Layer},
                                            s, β₁, β₂, β₀₀, convtest=5e-12)
    t1 = time_ns()
    nl = length(layers) # Number of layers.
    nl < 2 && error("Too few layers")
    k0sq = k0 * k0
    β²min = 1e-10*k0sq
    # Compute quantities defined in (5.20):
    μ̃ = 2 * layers[s].μᵣ * layers[s+1].μᵣ / (layers[s].μᵣ + layers[s+1].μᵣ) # Normalized to μ0.
    ϵ̄ = (layers[s].ϵᵣ + layers[s+1].ϵᵣ) / 2 # Normalized to ϵ0.
    area = 4π^2 / norm(β₁ × β₂) # unit cell area (m^2):
    # Obtain the appropriate Green's function expansion coefficient:
    c3s = c3_calc(k0, u, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ)
    c3sp1 = c3_calc(k0, u, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ)
    d3s = d3_calc(k0, u, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ)
    d3sp1 = d3_calc(k0, u, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ)
    
    # Factors used in summands:
    f1 = 2 * ϵ̄
    f3 = layers[s].ϵᵣ * c3s + layers[s+1].ϵᵣ * c3sp1
    p1 = 2 / μ̃
    p3 = d3s / layers[s].μᵣ + d3sp1 / layers[s+1].μᵣ

    m = mmax_list[2] ÷ 2
    table1 = OffsetArray(zeros(ComplexF64, 2m+1, 2m+1), -m:m, -m:m) 
    table2 = OffsetArray(zeros(ComplexF64, 2m+1, 2m+1), -m:m, -m:m)

    converged = false
    convrepeat = 40 # number of consecutive rings for which convergence must occur
    convlist = OffsetArray([false for i in 0:mmax_list[end]÷2], 0:mmax_list[end]÷2)
    mmax = mmax_list[1] ÷ 2
    mmax_old = -2
    test1 = test2 = 0.0
    first = true
    while mmax < mmax_list[end] # Convergence loop
        if first
            mmax = 0
            first = false
        else
            mmax = nextprod([2,3,5], 1 + mmax)
            while 0 ≠ rem(mmax,2)
                mmax = nextprod([2,3,5], 1 + mmax)
            end
        end
        mmaxo2 = mmax÷2
        mmax_oldo2 = mmax_old÷2
        # Fill the tables:
        Threads.@threads for r in (mmax_oldo2+1):mmaxo2
            ringsum1 = zero(eltype(table1))
            ringsum2 = zero(eltype(table1))
            for (m,n) in Ring(r)
                βmn = β₀₀ + m*β₁ + n*β₂   # Modal transverse wave vector
                β² = βmn ⋅ βmn # magnitude squared
                β² = max(β², β²min)  # Avoid singularity
                κmn² = β² + u*u #  Eq. (4.24)
                κmn = sqrt(κmn²) 
                # Compute region 1 mode parameters:
                γ = mysqrt(β² - k0sq * layers[1].ϵᵣ * layers[1].μᵣ)
                # Calculate TE modal admittance of Region 1, multiplied by jωμ₀:
                YleftTE = γ / layers[1].μᵣ 
                # Calculate TM modal admittance of Region 1, divided by jωϵ₀:
                YleftTM = layers[1].ϵᵣ / γ
                # Step left-looking impedances up to junction s using Eq. (5.24):
                for i in 2:s
                    γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                    Y0TE = γ / layers[i].μᵣ  # TE modal admittance (normalized)
                    Y0TM = layers[i].ϵᵣ / γ  # TM modal admittance (normalized)
                    tanhi = tanh(layers[i].width * γ)
                    YleftTE = Y0TE * (YleftTE + Y0TE*tanhi) / (Y0TE + YleftTE * tanhi)  # Eq. (5.24d)
                    YleftTM = Y0TM * (YleftTM + Y0TM*tanhi) / (Y0TM + YleftTM * tanhi)  # Eq. (5.24d)
                end
                # Compute region nl mode parameters:
                γ = mysqrt(β² - k0sq * layers[nl].ϵᵣ * layers[nl].μᵣ)
                YrightTE = γ / layers[nl].μᵣ # TE modal admit. of Region nl multiplied by jωμ₀
                YrightTM = layers[nl].ϵᵣ / γ # TM modal admit. of Region nl divided by jωϵ₀
                # Step right-looking admittances down to junction s using Eq. (5.24b):
                for i in nl-1:-1:s+1
                    γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                    Y0TE = γ / layers[i].μᵣ # TE modal admittance (normalized)
                    Y0TM = layers[i].ϵᵣ / γ # TM modal admittance (normalized)
                    tanhi = tanh(layers[i].width * γ) 
                    YrightTE = Y0TE * (YrightTE + Y0TE*tanhi) / (Y0TE + YrightTE * tanhi)
                    YrightTM = Y0TM * (YrightTM + Y0TM*tanhi) / (Y0TM + YrightTM * tanhi)
                end
                # Compute summands (apart from phase factor and 1/A factor):
                table1[m,n] = YleftTM + YrightTM  - (f1 + f3/κmn²) / κmn  # Eq. (5.26a)
                ringsum1 += table1[m,n]
                table2[m,n] = ((YleftTM+YrightTM)*k0sq + YleftTE + YrightTE) / β² - 
                              (p1 + p3/κmn²) / κmn  # Eq. (5.26b)
                ringsum2 += table2[m,n]
            end
            # Check for convergence of this ring
            test1 = abs(ringsum1/table1[0,0])
            test2 = abs(ringsum2/table2[0,0])
            convlist[r] = test1 < convtest && test2 < convtest
        end

        #  Check for convergence
        if mmaxo2 ≥ convrepeat && all(@view convlist[(mmaxo2 - convrepeat + 1):mmaxo2])
            converged = true
            break
        else
            mmax_old = mmax
        end
    end

    !converged && @warn "Inadequate Convergence" test1 test2 convtest mmax maxlog=5

    # Create proper sized storage arrays for FFT routine:
    mmaxo2 = mmax÷2
    table1t = table1[-mmaxo2:mmaxo2-1,-mmaxo2:mmaxo2-1]
    table2t = table2[-mmaxo2:mmaxo2-1,-mmaxo2:mmaxo2-1]
    fft!(table1t)
    fft!(table2t) 
    # Adjust phase according to Equation (5.32).  Also, include factor of 1/area:
    for q in 0:mmax-1
        qterm = q * (π - ψ₂/mmax)
        for p in 0:mmax-1
            pterm = p * (π - ψ₁/mmax)
            cfact = cis(pterm+qterm) / area
            table1t[p+1,q+1] *= cfact
            table2t[p+1,q+1] *= cfact
        end
    end
    
    # Create proper sized interpolation array---Note that we add an extra row
    # and extra column at both the beginning and end of each table to allow 
    # for extra points needed in the interpolation scheme.
    table1 = OffsetArray(zeros(ComplexF64, mmax+2, mmax+2), -1:mmax, -1:mmax)
    table1[0:mmax-1, 0:mmax-1] = table1t
    table2 = OffsetArray(zeros(ComplexF64, mmax+2, mmax+2), -1:mmax, -1:mmax)
    table2[0:mmax-1, 0:mmax-1] = table2t
    # Add extra row and column to cover all the way to ξ=1 and η=1:
    table1[mmax, 0:mmax-1] = cis(-ψ₁) * table1[0, 0:mmax-1]
    table1[0:mmax-1, mmax] = cis(-ψ₂) * table1[0:mmax-1, 0]
    table1[mmax, mmax] = cis(-(ψ₁+ψ₂)) * table1[0,0]
    table2[mmax, 0:mmax-1] = cis(-ψ₁) * table2[0, 0:mmax-1]
    table2[0:mmax-1, mmax] = cis(-ψ₂) * table2[0:mmax-1, 0]
    table2[mmax, mmax] = cis(-(ψ₁+ψ₂)) * table2[0,0]
    # Add extra row and column to cover all the way to ξ=-1/mmax and η=-1/mmax:
    table1[-1, 0:mmax] = cis(ψ₁) * table1[mmax-1, 0:mmax]
    table1[-1:mmax-1, -1] = cis(ψ₂) * table1[-1:mmax-1, mmax-1]
    table2[-1, 0:mmax] = cis(ψ₁) * table2[mmax-1, 0:mmax]
    table2[-1:mmax-1, -1] = cis(ψ₂) * table2[-1:mmax-1, mmax-1]
    # Use meaningful names (p for "prime"):
    Σpm1_func = make_Σm_func(table1, β₁, β₂, ψ₁, ψ₂)
    Σpm2_func = make_Σm_func(table2, β₁, β₂, ψ₁, ψ₂)
    t2 = time_ns()
    tsec = round((t2-t1)/1e9; digits=tdigits)
    @info "      $tsec seconds to compute $mmax × $mmax magnetic modal tables"
    return (Σpm1_func, Σpm2_func)
end


"""
    direct_electric_modal_series(k0,u,ψ₁,ψ₂,layers,s β₁,β₂,β₀₀,ρdif) --> (Σm1,Σm2)
   
This function uses direct, brute-force summation to calculate the electric source
modal series needed for the potential Green's functions.  **THIS ROUTINE IS FOR
TEST PURPOSES ONLY!!!!**  This routine is **NOT** numerically efficient.
It is only used for comparison and testing purposes.

## Arguments:

- `k0`: Free-space wavenumber (1/meter).
- `u`:  Smoothing parameter (1/meter).
- `ψ₁`, `ψ₂`:  Unit cell incremental phase shifts (radians).
- `layers`:  An array of element type `Layer` containing the layer 
       parameters for the cascade structure.  Note that the first
       and last layer's thicknesses are not accounted for in this
       function.  They are assumed to be semi-infinite.
- `s`:  Interface number (within layers) at which the FSS or PSS sheet is located.
- `β₁`, `β₂`:   2-vectors containing the reciprocal lattice basis 
              vectors.  Units are (1/meters).
- `β₀₀`:    2-vector containing the principal (i.e. with index (0,0)) Floquet 
           vector transverse wavenumber which incorporates the intrinsic phase 
           shifts.  Units are (1/meters).
- `ρdif`:  2-vector containing the difference between observation and source 
           points.

##  Return Values
    
`Σm1`, `Σm2`:  Offset arrays of with indices 0:max_ring containing the direct 
               modal series defined in Eq. (5.19) of the theory documentation.  
               `Σm1[i]` contains the partial sum of ring `i`, and similarly 
               for `Σm2`.
"""
function direct_electric_modal_series(k0, u, ψ₁, ψ₂,
                                      layers::Vector{Layer}, s, β₁,β₂,β₀₀, ρdif)
    max_rings = 3100 # Max number of rings for direct_modal_series.
    Σm1 = OffsetArray(zeros(ComplexF64, max_rings+1), 0:max_rings) 
    Σm2 = OffsetArray(zeros(ComplexF64, max_rings+1), 0:max_rings) 

    mmax =  32 # Initial trial value
    mmax_old = 0 
    converged = false
    k0sq = k0 * k0
    nl = length(layers)
    # Compute quantities defined in Eq. (5.20):
    μ̃ = 2 * layers[s].μᵣ * layers[s+1].μᵣ / (layers[s].μᵣ + layers[s+1].μᵣ) # Normalized to μ0.
    ϵ̄ = (layers[s].ϵᵣ + layers[s+1].ϵᵣ) / 2 # Normalized to ϵ0.
    area = 4π^2 / norm(β₁ × β₂) # unit cell area (m^2):

    # Obtain the appropriate Green's function expansion coefficients:
    c3 = c3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)
    d3 = d3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)

    function summands(m,n)
        let β₀₀=β₀₀, u=u, layers=layers, s=s, nl=nl, k0sq=k0sq, c3=c3, d3=d3, ϵ̄=ϵ̄, μ̃=μ̃, β₁=β₁, β₂=β₂
            βmn = β₀₀ + m*β₁ + n*β₂ 
            β² = βmn ⋅ βmn
            β²min = 1e-10*k0sq
            β² = max(β², β²min)  # Avoid singularity
            κmn² = β² + u*u 
            κmn = sqrt(κmn²)
            # Compute region 1 mode parameters:
            γ = mysqrt(β² - k0sq * layers[1].ϵᵣ * layers[1].μᵣ)
            # Calculate TE modal impedance of Region 1, divided by (jωμ₀):
            ZleftTE = layers[1].μᵣ / γ
            # Calculate TM modal impedance of Region 1, multiplied by (jωϵ₀):
            ZleftTM = γ / layers[1].ϵᵣ
            # Step left-looking impedances up to junction s using (5.14):
            for i in 2:s
                γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                Z0TE = layers[i].μᵣ / γ   # TE modal impedance (normalized)
                Z0TM = γ / layers[i].ϵᵣ  # TM modal impedance (normalized)
                tanhi = tanh(layers[i].width * γ)
                # Eq. (5.14b):
                ZleftTE = Z0TE * (ZleftTE + Z0TE*tanhi) / (Z0TE + ZleftTE*tanhi)
                ZleftTM = Z0TM * (ZleftTM + Z0TM*tanhi) / (Z0TM + ZleftTM*tanhi)
            end
            # Compute region nl mode parameters:
            γ = mysqrt(β² - k0sq * layers[end].ϵᵣ * layers[end].μᵣ)
            ZrightTE = layers[end].μᵣ / γ # divided by (jωμ₀)
            ZrightTM = γ / layers[end].ϵᵣ # multiplied by (jωϵ₀)
            # Step right-looking impedances down to junction s using (5.14d):
            for i in nl-1:-1:s+1
                γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                Z0TE = layers[i].μᵣ / γ
                Z0TM = γ / layers[i].ϵᵣ  
                tanhi = tanh(layers[i].width * γ) 
                ZrightTE = Z0TE * (ZrightTE + Z0TE*tanhi) / (Z0TE + ZrightTE*tanhi)
                ZrightTM = Z0TM * (ZrightTM + Z0TM*tanhi) / (Z0TM + ZrightTM*tanhi)
            end
            # Compute (normalized) TLGF's using Eq. (5.13):
            ViTE = ZleftTE * ZrightTE / (ZleftTE + ZrightTE)
            ViTM = ZleftTM * ZrightTM / (ZleftTM + ZrightTM)
            # Compute summands (apart from phase factor and 1/(2A) factor):
            ViTE *= 2 / μ̃ 
            summand1 = ViTE - (1 + c3/κmn²) / κmn 
            ViTM *= 2 * ϵ̄ 
            summand2 = (ViTM + (k0sq * ϵ̄ * μ̃) * ViTE) / β² - (1 + d3/κmn²) / κmn
            return summand1, summand2
        end
    end


    Σ1, Σ2 = summands(0,0)
    Σm1[0] = Σ1 * cis(-(β₀₀ ⋅ ρdif))
    Σm2[0] = Σ2 * cis(-(β₀₀ ⋅ ρdif))
    
    # Begin loop over summation lattice rings.
    for r in 1:max_rings
        sumring1 = zero(ComplexF64)
        sumring2 = zero(ComplexF64)
        for (m,n) in Ring(r)
            cfact = cis(-((β₀₀ + m*β₁ + n*β₂) ⋅ ρdif))
            Σ1, Σ2 = summands(m,n)
            sumring1 += Σ1 * cfact
            sumring2 += Σ2 * cfact
        end
        # Done with ring r.  Add ring contributions to sums.
        #Σm1[r] = Σm1[r-1] + sumring1 
        #Σm2[r] = Σm2[r-1] + sumring2 
        Σm1[r] = sumring1 
        Σm2[r] = sumring2 
    end
#    cumsum!(Σm1, Σm1)
#    cumsum!(Σm2, Σm2)
    Σm1 /= (2 * area)
    Σm2 /= (2 * area)
    return Σm1, Σm2
end

"""
    direct_magnetic_modal_series(k0,u,ψ₁,ψ₂,layers,s β₁,β₂,β₀₀,ρdif) --> (Σpm1,Σpm2)
   
This function uses direct, brute-force summation to calculate the magnetic 
source modal series needed for the potential Green's functions.  
**THIS ROUTINE IS FOR TEST PURPOSES ONLY!!!!**  This routine is **NOT** 
numerically efficient. It is only used for comparison and testing purposes.

## Arguments:

- `k0`: Free-space wavenumber (1/meter).
- `u`:  Smoothing parameter (1/meter).
- `ψ₁`, `ψ₂`:  Unit cell incremental phase shifts (radians).
- `layers`:  An array of element type `Layer` containing the layer 
       parameters for the cascade structure.  Note that the first
       and last layer's thicknesses are not accounted for in this
       function.  They are assumed to be semi-infinite.
- `s`:  Interface number (within layers) at which the FSS or PSS sheet is located.
- `β₁`, `β₂`:   2-vectors containing the reciprocal lattice basis 
              vectors.  Units are (1/meters).
- `β₀₀`:    2-vector containing the principal (i.e. with index (0,0)) Floquet 
           vector transverse wavenumber which incorporates the intrinsic phase 
           shifts.  Units are (1/meters).
- `ρdif`:  2-vector containing the difference between observation and source 
           points.

##  Return Values
    
`Σm1`, `Σm2`:  Offset arrays of with indices 0:max_rings containing the direct 
                 modal series defined  in Eq. (5.26) of the theory documentation.  
                `Σm1[i]` contains the partial sum of ring `i`, and similarly 
                 for `Σm2`.
"""
function direct_magnetic_modal_series(k0, u, ψ₁, ψ₂,
                                      layers::Vector{Layer}, s, β₁,β₂,β₀₀, ρdif)
    max_rings = 3100 # Max number of rings for direct_modal_series.
    Σm1 = OffsetArray(zeros(ComplexF64, max_rings+1), 0:max_rings) 
    Σm2 = OffsetArray(zeros(ComplexF64, max_rings+1), 0:max_rings) 

    k0sq = k0 * k0
    nl = length(layers)
    # Compute quantities defined in Eq. (5.20):
    μ̃ = 2 * layers[s].μᵣ * layers[s+1].μᵣ / (layers[s].μᵣ + layers[s+1].μᵣ) # Normalized to μ0.
    ϵ̄ = (layers[s].ϵᵣ + layers[s+1].ϵᵣ) / 2 # Normalized to ϵ0.
    area = 4π^2 / norm(β₁ × β₂) # unit cell area (m^2):

    # Obtain the appropriate Green's function expansion coefficients:
    c3s = c3_calc(k0, u, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ)
    c3sp1 = c3_calc(k0, u, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ)
    d3s = d3_calc(k0, u, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ)
    d3sp1 = d3_calc(k0, u, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ)
    
    # Factors used in summands:
    f1 = 2 * ϵ̄
    f3 = layers[s].ϵᵣ * c3s + layers[s+1].ϵᵣ * c3sp1
    p1 = 2 / μ̃
    p3 = d3s / layers[s].μᵣ + d3sp1 / layers[s+1].μᵣ


    function summands(m,n)
        let β₀₀=β₀₀, u=u, layers=layers, s=s, nl=nl, k0sq=k0sq, c3s=
            c3s, c3sp1=c3sp1, d3s=d3s, d3sp1=d3sp1, f1=f1, f3=f3, p1=p1, p3=p3,
            β₁=β₁, β₂=β₂
            
            βmn = β₀₀ + m*β₁ + n*β₂ 
            β² = βmn ⋅ βmn
            β²min = 1e-10*k0sq
            β² = max(β², β²min)  # Avoid singularity
            κmn² = β² + u*u 
            κmn = sqrt(κmn²)
            # Compute region 1 mode parameters:
            γ = mysqrt(β² - k0sq * layers[1].ϵᵣ * layers[1].μᵣ)
            # Calculate TE modal admittance of Region 1, multiplied by jωμ₀:
            YleftTE = γ / layers[1].μᵣ 
            # Calculate TM modal admittance of Region 1, divided by jωϵ₀:
            YleftTM = layers[1].ϵᵣ / γ
            # Step left-looking impedances up to junction s using Eq. (5.24):
            for i in 2:s
                γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                Y0TE = γ / layers[i].μᵣ  # TE modal admittance (normalized)
                Y0TM = layers[i].ϵᵣ / γ  # TM modal admittance (normalized)
                tanhi = tanh(layers[i].width * γ)
                YleftTE = Y0TE * (YleftTE + Y0TE*tanhi) / (Y0TE + YleftTE * tanhi)  # Eq. (5.24d)
                YleftTM = Y0TM * (YleftTM + Y0TM*tanhi) / (Y0TM + YleftTM * tanhi)  # Eq. (5.24d)
            end
            # Compute region nl mode parameters:
            γ = mysqrt(β² - k0sq * layers[nl].ϵᵣ * layers[nl].μᵣ)
            YrightTE = γ / layers[nl].μᵣ # TE modal admit. of Region nl multiplied by jωμ₀
            YrightTM = layers[nl].ϵᵣ / γ # TM modal admit. of Region nl divided by jωϵ₀
            # Step right-looking impedances down to junction s using Eq. (5.24b):
            for i in nl-1:-1:s+1
                γ = mysqrt(β² - k0sq * layers[i].ϵᵣ * layers[i].μᵣ)
                Y0TE = γ / layers[i].μᵣ # TE modal admittance (normalized)
                Y0TM = layers[i].ϵᵣ / γ # TM modal admittance (normalized)
                tanhi = tanh(layers[i].width * γ) 
                YrightTE = Y0TE * (YrightTE + Y0TE*tanhi) / (Y0TE + YrightTE * tanhi)
                YrightTM = Y0TM * (YrightTM + Y0TM*tanhi) / (Y0TM + YrightTM * tanhi)
            end
            # Compute summands (apart from phase factor and 1/A factor):
            summand1 = YleftTM + YrightTM  - (f1 + f3/κmn²) / κmn  # Eq. (5.26a)
            summand2 = ((YleftTM+YrightTM)*k0sq + YleftTE + YrightTE) / β² - 
                                             (p1 + p3/κmn²) / κmn  # Eq. (5.26b)
            return summand1, summand2
        end
    end


    Σ1, Σ2 = summands(0,0)
    Σm1[0] = Σ1 * cis(-(β₀₀ ⋅ ρdif))
    Σm2[0] = Σ2 * cis(-(β₀₀ ⋅ ρdif))
    # Begin loop over summation lattice rings.
    for r in 1:max_rings
        sumring1 = zero(ComplexF64)
        sumring2 = zero(ComplexF64)
        for (m,n) in Ring(r)
            cfact = cis(-((β₀₀ + m*β₁ + n*β₂) ⋅ ρdif))
            Σ1, Σ2 = summands(m,n)
            sumring1 += Σ1 * cfact
            sumring2 += Σ2 * cfact
        end
        # Done with ring r.  Add ring contributions to sums.
        #Σm1[r] = Σm1[r-1] + sumring1 
        #Σm2[r] = Σm2[r-1] + sumring2 
        Σm1[r] = sumring1 
        Σm2[r] = sumring2 
    end
    Σm1 /= area
    Σm2 /= area
    return Σm1, Σm2
end


end # module
