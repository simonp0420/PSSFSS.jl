module FillZY
export fillz, filly


using Statistics: mean
using StaticArrays: SMatrix, MVector, @SVector, SVector
using OffsetArrays
using MetalSurfaceImpedance: Zsurface
using Unitful # for ustrip and u"m"
using LinearAlgebra: norm, ⋅
using PSSFSS.Constants: μ₀, ϵ₀, c₀, twopi, fourpi, tol, tdigits
using PSSFSS.Layers: Layer
using PSSFSS.Sheets: RWGSheet, facecount
using PSSFSS.RWG: RWGData
using PSSFSS.PGF: c3_calc, d3_calc
using PSSFSS.Zint: zint, filljk!, vtxcrd
using PSSFSS.PGF: electric_modal_sum_funcs, magnetic_modal_sum_funcs
using PSSFSS.Log: @logfile
using OhMyThreads: @tasks, @set, @init, DynamicScheduler, StaticScheduler

const next = (2, 3, 1)
const prev = (3, 1, 2)
const third = SMatrix{3,3}([0 3 2; 3 0 1; 2 1 0])
const symbols = (:I1, :I1_ξ, :I1_η, :I2, :J, :J_ξ, :J_η, :K, :K_ξ, :K_η, :rinv, :ρ_r)

"""
    vertexcoords_opposite_edge(edge::Int, face::Int, apert::RWGData)

    Return the `SV2` vector for the vertex opposite the given edge of the face.
"""
function vertexcoords_opposite_edge(edge::Int, face::Int, apert::RWGSheet)
    edges = @view apert.fe[:, face]
    ρs = vtxcrd(face, apert)
    if edge == edges[1]
        return ρs[1]
    elseif edge == edges[2]
        return ρs[2]
    elseif edge == edges[3]
        return ρs[3]
    else
        error("Bad edge $edge and face $face")
    end
end

"""
    fillz(k0,u,layers::AbstractVector{Layer},s,ψ₁,ψ₂,metal::RWGSheet,rwgdat::RWGData) -> zmat

Fill the generalized impedance matrix for an FSS of electric current type.

## Arguments:

- `k0`: Free-space wavenumber (rad/meter).
- `u`:  Green's function smoothing factor (1/meter).
- `layers`:  An array characterizing the dielectric layers surrounding the FSS sheet.
- `s`:  An integer indexing the interface within layers at which the FSS
          sheet is located. `s=1` implies the sheet is between `layers[1]` and 
          `layers[2]`, etc.
- `ψ₁`,`ψ₂`:  Variables containing the unit cell incremental phase shifts (radians).
- `metal`:  A variable which characterizes the metalization region of the FSS/PSS.
- `rwgdat`:  A variable which defines the basis functions.

## Return value

- `zmat`: Complex array of size `(Nbf,Nbf)`, where `Nbf` is the 
          number of basis functions. On exit, this array will have been 
          filled with the generalized impedance matrix of the moment 
          method formulation.

"""
function fillz(k0, u, layers::AbstractVector{Layer}, s, ψ₁, ψ₂, metal::RWGSheet, rwgdat::RWGData)


    closed = true              # Always use singularity extraction.
    nbf = size(rwgdat.bfe, 2)
    zmat = rwgdat.zorymat
    zmat .= zero(eltype(zmat))

    # Initialize the Green's functions expansion coefficients:
    c3 = c3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)
    d3 = d3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)

    # Calculate beta00, the fundamental transverse wave vector (1/meter)
    units_per_meter = ustrip(Float64, metal.units, 1u"m")
    β₁, β₂ = metal.β₁ * units_per_meter, metal.β₂ * units_per_meter
    β₀₀ = (ψ₁ * β₁ + ψ₂ * β₂) / twopi

    # Initialize functions for modal series:
    (Σm1_func, Σm2_func) = electric_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, 1e-7)

    # floquet_factor is indexed into using values in rwgdat.eci:
    floquet_factor = OffsetArray(SVector(1.0, 1.0, cis(-ψ₁), 1.0, cis(-ψ₂)), 0:4)

    # Set up aliases into data structures for more convenient reference:
    bff = rwgdat.bff
    bfe = rwgdat.bfe
    eci = rwgdat.eci

    nface = facecount(metal)
    i2s = CartesianIndices((nface, nface))
    nbf = size(rwgdat.bfe, 2)
    ϵᵣ₁ = layers[s].ϵᵣ
    μᵣ₁ = layers[s].μᵣ
    ϵᵣ₂ = layers[s+1].ϵᵣ
    μᵣ₂ = layers[s+1].μᵣ
    μ̃ = 2 * μᵣ₁ * μᵣ₂ / (μᵣ₁ + μᵣ₂) # Equation (4-11) (normalized to μ₀)
    ϵ̄ = (ϵᵣ₁ + ϵᵣ₂) / 2              # Equation (4-23) (normalized to ϵ₀)
    ω = k0 * c₀                      # Radian frequency (Radians/second)
    jω = im * ω
    A_factor = μ₀ / fourpi * μ̃
    Φ_factor = im / (ϵ̄ * twopi * ω * ϵ₀)

    # Update Zs if it depends on conductivity
    Zs = metal.Zs
    if metal.σ > 0
        fhz = k0 * c₀ / twopi
        Zs = Zsurface(fhz, metal.σ, metal.Rq, metal.disttype)
        metal.Zs = Zs
    end

    # Check whether or not the frequency-independent face/face integrals are up to date:
    if ψ₁ ≠ metal.ψ₁ || ψ₂ ≠ metal.ψ₂ || metal.u == 0 ||
       abs((u / units_per_meter - metal.u) / metal.u) > tol
        # Set up the values of the face/face integrals' params stored in metal:
        metal.ψ₁ = ψ₁
        metal.ψ₂ = ψ₂
        metal.u = u / units_per_meter
        # Fill the frequency-independent face/face integrals:
        t_spatial = time()
        filljk!(metal, rwgdat, closed)
        t_spatial = time() - t_spatial
        @logfile "      $(round(t_spatial,digits=tdigits)) seconds for spatial face integrals"
    end

    nufp = rwgdat.nufp
    # Compute frequency-dependent face-face integrals:
    nufp == length(metal.I1) || (metal.I1 = zeros(ComplexF64, nufp))
    nufp == length(metal.I1_ξ) || (metal.I1_ξ = zeros(ComplexF64, nufp))
    nufp == length(metal.I1_η) || (metal.I1_η = zeros(ComplexF64, nufp))
    nufp == length(metal.I2) || (metal.I2 = zeros(ComplexF64, nufp))

    t1 = time_ns()
    nthr = Threads.nthreads()
    nchunks = 2 * nthr
    @tasks for iufp in 1:rwgdat.nufp   # Loop over each unique face pair
        @set scheduler = DynamicScheduler(; nchunks)
        ifmifs = rwgdat.ufp2fp[iufp][1]  # Obtain index into face/face matrix.
        ifm, ifs = Tuple(i2s[ifmifs])  # indices of match and source triangles
        # Obtain the coordinates (in meters) of the source triangle's vertices:
        rs = vtxcrd(ifs, metal) ./ units_per_meter
        rm = vtxcrd(ifm, metal) ./ units_per_meter  # Coords (m) of the match tri. vertices.
        rmc = mean(rm)        # Match face centroid (meters).
        # Perform frequency-dependent integrals over source triangle 
        (metal.I1[iufp], metal.I1_ξ[iufp], metal.I1_η[iufp], metal.I2[iufp]) = zint(Σm1_func, Σm2_func, rs, rmc)
    end
    t2 = time_ns()
    t_spectral = t2 - t1
    @logfile "      $(round(t_spectral/1e9,digits=tdigits)) seconds for spectral face integrals"


    # Fill the interaction matrix
    pm = (1, -1)
    nthr = Threads.nthreads()
    nchunks = 2 * nthr
    t1 = time_ns()
    @tasks for bfci in CartesianIndices((nbf,nbf))   # Loop over basis function pairs
        @set scheduler = DynamicScheduler(; nchunks)

        mbf, sbf = Tuple(bfci) # match and source basis function indices

        sfp, sfm = @view bff[:, sbf] # plus and minus faces of source basis function
        sep, sem = @view bfe[:, sbf] # plus and minus edges of source basis function
        mfp, mfm = @view bff[:, mbf] # plus and minus faces of match basis function
        mep, mem = @view bfe[:, mbf] # plus and minus edges of match basis function
        for (ss, sf, se) in zip(pm, (sfp, sfm), (sep, sem)) # source sign, source face, source edge

            # Obtain the coordinates (in meters) of the source triangle's vertices:
            rs = vtxcrd(sf, metal) ./ units_per_meter
            rs_opp = vertexcoords_opposite_edge(se, sf, metal) ./ units_per_meter

            rs32 = rs[3] - rs[2]; rs12 = rs[1] - rs[2]
            area = 0.5 * (rs32[1] * rs12[2] - rs32[2] * rs12[1]) # source triangle signed area
            area48 = 48 * abs(area)  # Needed for surface loading

            source_flag = ss * floquet_factor[eci[se]]
            A_source_flag = A_factor * source_flag
            Φ_source_flag = Φ_factor * source_flag

            for (ms, mf, me) in zip(pm, (mfp,mfm), (mep,mem)) # match sign, match face, match edge
                match_flag = ms * conj(floquet_factor[eci[me]])
                rm = vtxcrd(mf, metal) ./ units_per_meter  # Coords (m) of the match tri. vertices.
                rmc = mean(rm)        # Match face centroid (meters).

                # Compute vector from active vertex to centroid of match triangle
                # (divided by 2) as in Eq. (7-15):
                rm_opp = vertexcoords_opposite_edge(me, mf, metal) ./ units_per_meter
                ρc2 = 0.5 * (rmc - rm_opp)

                iufp = rwgdat.ufpm[mf, sf] # Obtain unique face pair index
        
                # Recall the spatial face integrals:
                I1 = metal.I1[iufp]; I1_ξ = metal.I1_ξ[iufp]; I1_η = metal.I1_η[iufp]
                I2 = metal.I2[iufp]
                J = metal.J[iufp]; J_ξ = metal.J_ξ[iufp]; J_η = metal.J_η[iufp]
                K = metal.K[iufp]; K_ξ = metal.K_ξ[iufp]; K_η = metal.K_η[iufp]
                rinv = metal.rinv[iufp]; ρ_r = metal.ρ_r[iufp]
                I1_ζ = I1 - I1_ξ - I1_η
                J_ζ = J - J_ξ - J_η
                K_ζ = K - K_ξ - K_η

                # Compute singular contribution for this edge (the middle term in 
                # square brackets in Equation (7-21) using (B.2):
                Asing = ρ_r + u * rinv * (rmc - rs_opp)

                # Compute Eq. (7-26) (but correct sign is carried in source flags):
                I1_i = rs[1] * I1_ξ + rs[2] * I1_η + rs[3] * I1_ζ - rs_opp * I1
                J_i = rs[1] * J_ξ + rs[2] * J_η + rs[3] * J_ζ - rs_opp * J
                K_i = rs[1] * K_ξ + rs[2] * K_η + rs[3] * K_ζ - rs_opp * K

                # Compute Equation (7-21):
                A_i = A_source_flag * (fourpi * I1_i + Asing + u * J_i + c3 / u * K_i)

                # Compute Equation (7-31):
                Φ_i = Φ_source_flag * (fourpi * I2 + u * (rinv + J) + d3 / u * K)

                # Compute the dot product in Eq (7-15) apart from sign:
                dotprod = ρc2 ⋅ A_i

                # Add contribution to the impedance matrix
                zmat[mbf, sbf] += match_flag * (jω * dotprod - Φ_i)

                # Add surface loading, if applicable:
                if sf == mf && !iszero(Zs)  
                    ρ2s = metal.ρ[metal.e2[se]] / units_per_meter
                    ρ1s = metal.ρ[metal.e1[se]] / units_per_meter
                    ls = norm(ρ2s - ρ1s)
                    if me == se # Self edge
                        lother1 = norm(rs_opp - ρ1s)
                        lother2 = norm(rs_opp - ρ2s)
                        Zload = Zs / area48 *
                                (3 * (lother1^2 + lother2^2) - ls^2)   # Eq. (7-34)
                    else
                        ρ2m = metal.ρ[metal.e2[me]] / units_per_meter
                        ρ1m = metal.ρ[metal.e1[me]] / units_per_meter
                        lm = norm(ρ2m - ρ1m)
                        ρ1other = vertexcoords_opposite_edge(se,sf,metal) / units_per_meter
                        ρ2other = vertexcoords_opposite_edge(me,mf,metal) / units_per_meter
                        lother = norm(ρ2other - ρ1other)
                        Zload = Zs / area48 * source_flag * match_flag *
                                (ls^2 + lm^2 - 3 * lother^2) # Eq. (7-35)
                    end
                    zmat[mbf,sbf] += Zload
                end # test for surface loading
            end # match sign, match face, match edge
        end # source sign, source face, source edge
    end # # Loop over each unique basis function pairs
    t2 = time_ns()
    tsec = round((t2 - t1) / 1e9; digits=tdigits)
    @logfile "      $tsec seconds to fill $(size(zmat,1)) × $(size(zmat,2)) matrix entries"
    return zmat
end


"""
    filly(k0,u,layers::AbstractVector{Layer},s,ψ₁,ψ₂,apert::RWGSheet,rwgdat::RWGData) -> ymat

Fill the generalized impedance matrix for an FSS of electric current type.

## Arguments:

- `k0`: Free-space wavenumber (rad/meter).
- `u`:  Green's function smoothing factor (1/meter).
- `layers`:  An array characterizing the dielectric layers surrounding the FSS sheet.
- `s`:  An integer indexing the interface within layers at which the FSS
          sheet is located. `s=1` implies the sheet is between `layers[1]` and 
          `layers[2]`, etc.
- `ψ₁`,`ψ₂`:  Variables containing the unit cell incremental phase shifts (radians).
- `apert`:  A variable which characterizes the aperture region of the FSS/PSS.
- `rwgdat`:  A variable which defines the basis functions.

## Return value

- `ymat`: Complex array of size `(Nbf,Nbf)`, where `Nbf` is the 
          number of basis functions. On exit, this array will have been 
          filled with the generalized admittance matrix of the moment 
          method formulation.

"""
function filly(k0, u, layers::AbstractVector{Layer}, s, ψ₁, ψ₂, apert, rwgdat)

    closed = true              # Always use singularity extraction.
    nbf = size(rwgdat.bfe, 2)
    ymat = rwgdat.zorymat
    ymat .= zero(eltype(ymat))

    # Initialize the Green's functions expansion coefficients:
    c3s = c3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ)
    c3sp1 = c3_calc(k0, u, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)
    d3s = d3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ)
    d3sp1 = d3_calc(k0, u, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)

    # Calculate beta00, the fundamental transverse wave vector (1/meter)
    units_per_meter = ustrip(Float64, apert.units, 1u"m")
    β₁, β₂ = apert.β₁ * units_per_meter, apert.β₂ * units_per_meter
    β₀₀ = (ψ₁ * β₁ + ψ₂ * β₂) / twopi

    # Initialize functions for modal series:
    (Σm1_func, Σm2_func) = magnetic_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, 1e-7)

    # floquet_factor is indexed into using values in rwgdat.eci:
    floquet_factor = OffsetArray(SVector(1.0, 1.0, cis(-ψ₁), 1.0, cis(-ψ₂)), 0:4)

    # Set up aliases into data structures for more convenient reference:
    bfe = rwgdat.bfe
    bff = rwgdat.bff
    eci = rwgdat.eci

    nface = facecount(apert)
    i2s = CartesianIndices((nface, nface))
    nbf = size(bfe, 2)
    ϵᵣ₁ = layers[s].ϵᵣ
    μᵣ₁ = layers[s].μᵣ
    ϵᵣ₂ = layers[s+1].ϵᵣ
    μᵣ₂ = layers[s+1].μᵣ
    μ̃ = 2.0 * μᵣ₁ * μᵣ₂ / (μᵣ₁ + μᵣ₂) # Eq (4-11) (normalized to μ₀)
    ϵ̄ = 0.5 * (ϵᵣ₁ + ϵᵣ₂) # Equation (1-23) (normalized to ϵ₀)
    ω = k0 * c₀              # Radian frequency (Radians/second)
    jω = im * ω
    I1fact = π / ϵ̄
    KFfact = (c3s * ϵᵣ₁ + c3sp1 * ϵᵣ₂) / (2 * ϵ̄ * u)
    I2fact = π * μ̃
    KPfact = μ̃ / (2 * u) * (d3s / μᵣ₁ + d3sp1 / μᵣ₂)
    F_factor = -ϵ₀ / π * ϵ̄ 
    Ψ_factor = 2im / (π * ω * μ₀ * μ̃)

    # Check whether or not the apert face/face integrals are up to date:
    if ψ₁ ≠ apert.ψ₁ || ψ₂ ≠ apert.ψ₂ || apert.u == 0 ||
       abs((u / units_per_meter - apert.u) / apert.u) > tol
        # Set up the values of the face/face integrals' params stored in apert:
        apert.ψ₁ = ψ₁
        apert.ψ₂ = ψ₂
        apert.u = u / units_per_meter # Units are 1/(local length units)
        # Fill the frequency-independent face/face integrals:
        t_spatial = time()
        filljk!(apert, rwgdat, closed)
        t_spatial = time() - t_spatial
        @logfile "      $(round(t_spatial,digits=tdigits)) seconds for spatial face integrals"
    end

    nufp = rwgdat.nufp
    # Compute frequency-dependent face-face integrals:
    nufp == length(apert.I1) || (apert.I1 = zeros(ComplexF64, nufp))
    nufp == length(apert.I1_ξ) || (apert.I1_ξ = zeros(ComplexF64, nufp))
    nufp == length(apert.I1_η) || (apert.I1_η = zeros(ComplexF64, nufp))
    nufp == length(apert.I2) || (apert.I2 = zeros(ComplexF64, nufp))

    t1 = time_ns()
    nthr = Threads.nthreads()
    nchunks = 2 * nthr
    @tasks for iufp in 1:rwgdat.nufp   # Loop over each unique face pair
        @set scheduler = DynamicScheduler(; nchunks)
        ifmifs = rwgdat.ufp2fp[iufp][1]  # Obtain index into face/face matrix.
        ifm, ifs = Tuple(i2s[ifmifs])  # indices of match and source triangles
        # Obtain the coordinates (in meters) of the source triangle's vertices:
        rs = vtxcrd(ifs, apert) ./ units_per_meter
        rm = vtxcrd(ifm, apert) ./ units_per_meter  # Coords (m) of the match tri. vertices.
        rmc = mean(rm)        # Match face centroid (meters).
        # Perform frequency-dependent integrals over source triangle 
        (apert.I1[iufp], apert.I1_ξ[iufp], apert.I1_η[iufp], apert.I2[iufp]) = zint(Σm1_func, Σm2_func, rs, rmc)
    end
    t2 = time_ns()
    t_spectral = t2 - t1
    @logfile "      $(round(t_spectral/1e9,digits=tdigits)) seconds for spectral face integrals"

    # Fill the interaction matrix
    pm = (1, -1)
    nthr = Threads.nthreads()
    nchunks = 2 * nthr
    t1 = time_ns()

    @tasks for bfci in CartesianIndices((nbf,nbf))
        @set scheduler = DynamicScheduler(; nchunks)
        
        mbf, sbf = Tuple(bfci) # match and source basis function indices

        sfp, sfm = @view bff[:, sbf] # plus and minus faces of source basis function
        sep, sem = @view bfe[:, sbf] # plus and minus edges of source basis function
        mfp, mfm = @view bff[:, mbf] # plus and minus faces of match basis function
        mep, mem = @view bfe[:, mbf] # plus and minus edges of match basis function
        for (ss, sf, se) in zip(pm, (sfp, sfm), (sep, sem)) # source sign, source face, source edge
            # Obtain the coordinates (in meters) of the source triangle's vertices:
            rs = vtxcrd(sf, apert) ./ units_per_meter
            rs_opp = vertexcoords_opposite_edge(se, sf, apert) ./ units_per_meter
            source_flag = ss * floquet_factor[eci[se]]
            F_source_flag = F_factor * source_flag
            Ψ_source_flag = Ψ_factor * source_flag

            for (ms, mf, me) in zip(pm, (mfp,mfm), (mep,mem)) # match sign, match face, match edge
                match_flag = ms * conj(floquet_factor[eci[me]])
                rm = vtxcrd(mf, apert) ./ units_per_meter  # Coords (m) of the match tri. vertices.
                rmc = mean(rm)        # Match face centroid (meters).
                iufp = rwgdat.ufpm[mf, sf] # Obtain unique face pair index

                # Recall the spatial face integrals:
                I1 = apert.I1[iufp]; I1_ξ = apert.I1_ξ[iufp]; I1_η = apert.I1_η[iufp]
                I2 = apert.I2[iufp]
                J = apert.J[iufp]; J_ξ = apert.J_ξ[iufp]; J_η = apert.J_η[iufp]
                K = apert.K[iufp]; K_ξ = apert.K_ξ[iufp]; K_η = apert.K_η[iufp]
                rinv = apert.rinv[iufp]; ρ_r = apert.ρ_r[iufp]
                I1_ζ = I1 - I1_ξ - I1_η
                J_ζ = J - J_ξ - J_η
                K_ζ = K - K_ξ - K_η

                # Compute vector from active vertex to centroid of match triangle
                # (divided by 2) as in Eq. (7-53):
                rm_opp = vertexcoords_opposite_edge(me, mf, apert) ./ units_per_meter
                ρc2 = 0.5 * (rmc - rm_opp)
                    
                # Compute singular contribution for this edge (the middle term in 
                # square brackets in Equation (7-57)) using (B-2):
                Fsing = ρ_r + u * rinv * (rmc - rs_opp)

                # Compute Eq. (7-26) (but correct sign is carried in source flags):
                I1_i = rs[1] * I1_ξ + rs[2] * I1_η + rs[3] * I1_ζ - rs_opp * I1
                J_i = rs[1] * J_ξ + rs[2] * J_η + rs[3] * J_ζ - rs_opp * J
                K_i = rs[1] * K_ξ + rs[2] * K_η + rs[3] * K_ζ - rs_opp * K

                # Compute Equation (7-57):
                F_i = F_source_flag * (I1fact * I1_i + Fsing + u * J_i + KFfact * K_i)

                # Compute Equation (7-60):
                Ψ_i = Ψ_source_flag * (I2fact * I2 + u * (rinv + J) + KPfact * K)

                # Compute one of the dot products in Eq (7-53) apart from sign:
                dotprod = ρc2 ⋅ F_i

                # Add contribution to the admittance matrix
                ymat[mbf,sbf] += match_flag * (-im*ω*dotprod - Ψ_i)

            end  # match sign, match face, match edge
        end  # source sign, source face, source edge
    end  # threaded loop over basis function match/source pairs

    t2 = time_ns()
    tsec = round((t2 - t1) / 1e9; digits=tdigits)
    @logfile "      $tsec seconds to fill $(size(ymat,1)) × $(size(ymat,2)) matrix entries"

    return ymat
end # function


end # module
