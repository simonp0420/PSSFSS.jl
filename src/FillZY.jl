module FillZY
export fillz, filly


using Statistics: mean
using StaticArrays: SMatrix
using OffsetArrays
using Unitful # for ustrip and u"m"
using LinearAlgebra: norm, ⋅
using ..Constants: μ₀, ϵ₀, c₀, twopi, fourpi, tol
using ..Layers: Layer
using ..Sheets: RWGSheet #, MV2, SV2
using ..RWG: RWGData
using ..PGF: c3_calc, d3_calc
using ..Zint: zint, filljk!, vtxcrd
using ..PGF: electric_modal_sum_funcs, magnetic_modal_sum_funcs

const next = (2,3,1)
const prev = (3,1,2)
const third = SMatrix{3,3}([0 3 2; 3 0 1; 2 1 0])



"""
    fillz(k0,u,layers::Vector{Layer},s,ψ₁,ψ₂,metal::RWGSheet,rwgdat::RWGData) -> zmat

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
function fillz(k0,u,layers::Vector{Layer},s,ψ₁,ψ₂,metal::RWGSheet,rwgdat::RWGData)
    

    closed = true              # Always use singularity extraction.
    nbf = size(rwgdat.bfe, 2)
    zmat = zeros(ComplexF64, nbf, nbf)

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
    floquet_factor = OffsetArray(zeros(ComplexF64, 5), 0:4) 
    floquet_factor[0] = 1.0       # Edges not on unit cell boundary
    floquet_factor[1] = 1.0       # Edges at ξ=0 boundary.
    floquet_factor[2] = cis(-ψ₁)  # Edges at ξ=1 boundary.
    floquet_factor[3] = 1.0       # Edges at η=0 boundary.
    floquet_factor[4] = cis(-ψ₂)  # Edges at η=1 boundary.
    
    # Set up aliases into data structures for more convenient reference:
    bfe = rwgdat.bfe
    bff = rwgdat.bff
    ebf = rwgdat.ebf
    eci = rwgdat.eci
    
    nface = size(metal.fv,2)
    i2s = CartesianIndices((nface,nface))
    nbf = size(rwgdat.bfe,2)
    ϵᵣ₁ = layers[s].ϵᵣ
    μᵣ₁ = layers[s].μᵣ
    ϵᵣ₂ = layers[s+1].ϵᵣ
    μᵣ₂ = layers[s+1].μᵣ
    μ̃ = 2 * μᵣ₁ * μᵣ₂ / (μᵣ₁ + μᵣ₂) # Equation (4-11) (normalized to μ₀)
    ϵ̄ = (ϵᵣ₁ + ϵᵣ₂)/2              # Equation (4-23) (normalized to ϵ₀)
    ω = k0 * c₀                  # Radian frequency (Radians/second)
    
    # Check whether or not the frequency-independent face/face integrals are up to date:
    if ψ₁ ≠ metal.ψ₁ ||  ψ₂ ≠ metal.ψ₂ || metal.u == 0 ||
        abs((u / units_per_meter - metal.u) / metal.u) > tol 
        # Set up the values of the face/face integrals' params stored in metal:
        metal.ψ₁ = ψ₁
        metal.ψ₂ = ψ₂
        metal.u = u / units_per_meter 
        # Fill the frequency-independent face/face integrals:
        t_spatial = time()
        filljk!(metal, rwgdat, closed)
        t_spatial = time() - t_spatial
        @info "Spatial face integrals used $(round(t_spatial,digits=5)) sec"
    end

    t1 = time_ns()
    for iufp in 1:rwgdat.nufp  # Loop over each unique face pair
        ifmifs = rwgdat.ufp2fp[iufp][1]  # Obtain index into face/face matrix
        rowcol = i2s[ifmifs]
        ifm, ifs = rowcol[1], rowcol[2] # indices of match and source triangles
        is = @view metal.fe[:,ifs]     # Obtain the three edges of the source triangle.
        # Obtain the coordinates (in meters) of the source triangle's vertices:
        rs = vtxcrd(ifs, metal) ./ units_per_meter
        # Calculate the signed area of the source triangle:
        rs32 = rs[3] - rs[2]
        rs12 = rs[1] - rs[2]
        area = 0.5 * (rs32[1] * rs12[2] - rs32[2] * rs12[1])
        area48 = 48 * abs(area)  # Needed for surface loading
      
        self_tri = (ifm == ifs) # Source and match tri are the same?
      
        if self_tri
            # Obtain the length of each side of the source=match
            # triangle. To be used later in surface loading calculation.
            ls = (norm(rs32), norm(rs[1] - rs[3]), norm(rs12))
        end
      
        clsflg = self_tri || closed  # Extract if self or if always.
      
        i_m = @view metal.fe[:,ifm]   # Obtain the three edges of the match triangle.
        rm = vtxcrd(ifm, metal) ./ units_per_meter # Coordinates (m) of the match tri. vertices.
        rmc = mean(rm)        # Match face centroid (meters).

        # Perform frequency-dependent integrals over source triangle 
        (I1, I1_ξ, I1_η, I2) = zint(Σm1_func, Σm2_func, rs, rmc)
        # Recall the frequency-independent spatial integrals:
        J = metal.J[iufp]
        J_ξ = metal.J_ξ[iufp]
        J_η = metal.J_η[iufp]
        K = metal.K[iufp]
        K_ξ = metal.K_ξ[iufp]
        K_η = metal.K_η[iufp]
        rinv = metal.rinv[iufp]
        ρ_r = metal.ρ_r[iufp]
        I1_ζ = I1 - I1_ξ - I1_η
        J_ζ = J - J_ξ - J_η
        K_ζ = K - K_ξ - K_η
      
        # Compute vector from each vertex to centroid of match triangle
        # (divided by 2) as in Eq. (7-15):
        ρc2 =  0.5 * [rmc - rm[i] for i in 1:3]
      
        # Loop over the face pairs in this equivalence class:
        for (i0, ifmifs) in enumerate(rwgdat.ufp2fp[iufp])
            rowcol = i2s[ifmifs]
            ifm, ifs = rowcol[1], rowcol[2] # indices of match and source triangles
            is = @view metal.fe[:,ifs]      # Obtain the three edges of source tri.
            i_m = @view metal.fe[:,ifm]      # Obtain the three edges of the match tri.
            # Loop over each edge of the source triangle.
            for isl in 1:3   # isl is local edge index, is[isl] is global edge index.
                sbf = ebf[is[isl]]  # Source basis function for this source edge.
                sbf == 0 && continue  # no basis func. is defined
                if bff[1,sbf] == ifs # "Plus" source triangle:
                    source_flag = floquet_factor[eci[is[isl]]]
                elseif bff[2,sbf] == ifs  # "Minus" source triangle:
                    source_flag = -floquet_factor[eci[is[isl]]]
                else
                    error("Impossible situation in source edge loop")
                end
                A_source_flag = μ₀ / fourpi * source_flag * μ̃ 
                Φ_source_flag = im / ϵ̄ * (source_flag / (twopi*ω*ϵ₀))

                # Compute singular contribution for this edge (the middle term in 
                # square brackets in Equation (7-21) using (B.2):
                Asing = ρ_r  +  u * rinv * (rmc - rs[isl])
        
                # Compute Eq. (7-26) (but correct sign is carried in source flags):
                I1_i = rs[1] * I1_ξ + rs[2] * I1_η + rs[3] * I1_ζ - rs[isl] * I1
                J_i = rs[1] * J_ξ + rs[2] * J_η + rs[3] * J_ζ - rs[isl] * J
                K_i = rs[1] * K_ξ + rs[2] * K_η + rs[3] * K_ζ - rs[isl] * K
          
                # Compute Equation (7-21):
                A_i = A_source_flag * (fourpi * I1_i + Asing + u * J_i + c3 / u * K_i)
                # Compute Equation (7-31):
                Φ_i = Φ_source_flag * (fourpi * I2 + u * (rinv + J) + d3 / u * K)
          
                # Now loop over each edge of the match triangle:
                for iml in 1:3 # iml is local edge index, i_m[iml] is global index.
                    self_edge = i_m[iml] == is[isl] # Match edge = Source edge?
                    mbf = ebf[i_m[iml]]  # Match basis function for this match edge.
                    mbf == 0 && continue  # no basis function defined for this match tri edge
                    if bff[1,mbf] == ifm   # "Plus" match triangle.
                        match_flag = conj(floquet_factor[eci[i_m[iml]]])
                    elseif bff[2,mbf] == ifm  # "Minus" match tri.
                        match_flag = -conj(floquet_factor[eci[i_m[iml]]])
                    else
                        error("Impossible situation in match edge loop")
                    end
            
                    # Compute one of the dot products in Eq (7-15) apart from sign:
                    dotprod = ρc2[iml] ⋅ A_i
                    # Add contribution to the impedance matrix
                    zmat[mbf,sbf] += match_flag * (im*ω*dotprod - Φ_i)
                    # Add surface loading, if applicable:
                    if self_tri && metal.fr[ifm] ≠ 0
                        if self_edge
                            Zload = metal.fr[ifs] / area48 * 
                                (3*(ls[next[isl]]^2 + ls[prev[isl]]^2) - ls[isl]^2) * one(ComplexF64)  # Eq. (7-34)
                        else
                            Zload = metal.fr[ifs] / area48 * source_flag * match_flag * 
                                (ls[isl]^2  +  ls[iml]^2  - 3 * ls[third[isl,iml]]^2) # Eq. (7-35)
                        end
                        zmat[mbf,sbf] += Zload
                    end
                end 
            end # loop over isl, source triangle edges
        end # face pairs in this equivalence class
    end # facepair loop
    t2 = time_ns()
    tsec = round((t2-t1)/1e9; digits=5)
    @info "$tsec seconds to fill $(size(zmat,1)) × $(size(zmat,2)) matrix entries"
    return zmat
end


"""
    filly(k0,u,layers::Vector{Layer},s,ψ₁,ψ₂,apert::RWGSheet,rwgdat::RWGData) -> ymat

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
function filly(k0, u, layers::Vector{Layer}, s, ψ₁, ψ₂, apert, rwgdat)

    closed = true              # Always use singularity extraction.

    nbf = size(rwgdat.bfe, 2)
    ymat = zeros(ComplexF64, nbf, nbf)

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
    (Σm1_func, Σm2_func) = magnetic_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀,1e-7)

    # floquet_factor is indexed into using values in rwgdat.eci:
    floquet_factor = OffsetArray(zeros(ComplexF64, 5), 0:4) 
    floquet_factor[0] = 1.0       # Edges not on unit cell boundary
    floquet_factor[1] = 1.0       # Edges at ξ=0 boundary.
    floquet_factor[2] = cis(-ψ₁)  # Edges at ξ=1 boundary.
    floquet_factor[3] = 1.0       # Edges at η=0 boundary.
    floquet_factor[4] = cis(-ψ₂)  # Edges at η=1 boundary.
    
    # Set up aliases into data structures for more convenient reference:
    bfe = rwgdat.bfe
    bff = rwgdat.bff
    ebf = rwgdat.ebf
    eci = rwgdat.eci
    
    nface = size(apert.fv,2)
    i2s = CartesianIndices((nface,nface))
    nbf = size(rwgdat.bfe,2)
    ϵᵣ₁ = layers[s].ϵᵣ
    μᵣ₁ = layers[s].μᵣ
    ϵᵣ₂ = layers[s+1].ϵᵣ
    μᵣ₂ = layers[s+1].μᵣ
    μ̃ = 2.0 * μᵣ₁ * μᵣ₂ / (μᵣ₁ + μᵣ₂) # Eq (4-11) (normalized to μ₀)
    ϵ̄ = 0.5 * (ϵᵣ₁ + ϵᵣ₂) # Equation (1-23) (normalized to ϵ₀)
    ω = k0 * c₀              # Radian frequency (Radians/second)
    I1fact = π / ϵ̄
    KFfact = (c3s*ϵᵣ₁ + c3sp1*ϵᵣ₂) / (2 * ϵ̄ * u)
    I2fact = π * μ̃
    KPfact = μ̃ / (2*u) * (d3s / μᵣ₁ + d3sp1 / μᵣ₂)
    
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
        @info "Spatial face integrals used $(round(t_spatial,digits=5)) sec"
    end

    t1 = time_ns()
    for iufp in 1:rwgdat.nufp  # Loop over each unique face pair
        ifmifs = rwgdat.ufp2fp[iufp][1]  # Obtain index into face/face matrix.
        rowcol = i2s[ifmifs]
        ifm, ifs = rowcol[1], rowcol[2] # indices of match and source triangles
        is = @view apert.fe[:,ifs]     # Obtain the three edges of the source triangle.
        # Obtain the coordinates (in meters) of the source triangle's vertices:
        rs = vtxcrd(ifs, apert) ./ units_per_meter
        # Calculate the signed area of the source triangle:
        rs32 = rs[3] - rs[2]
        rs12 = rs[1] - rs[2]
        area = 0.5 * (rs32[1] * rs12[2] - rs32[2] * rs12[1])
        area48 = 48 * abs(area)  # Needed for surface loading
      
        self_tri = (ifm == ifs) # Source and match tri are the same?
        clsflg = self_tri || closed  # Extract if self or if always.
      
        i_m = @view apert.fe[:,ifm]   # Obtain the three edges of the match triangle.
        rm = vtxcrd(ifm, apert) ./ units_per_meter  # Coords (m) of the match tri. vertices.
        rmc = mean(rm)        # Match face centroid (meters).
      
        # Perform frequency-dependent integrals over source triangle 
        (I1, I1_ξ, I1_η, I2) = zint(Σm1_func, Σm2_func, rs, rmc)
        # Recall the frequency-independent spatial integrals:
        J = apert.J[iufp]
        J_ξ = apert.J_ξ[iufp]
        J_η = apert.J_η[iufp]
        K = apert.K[iufp]
        K_ξ = apert.K_ξ[iufp]
        K_η = apert.K_η[iufp]
        rinv = apert.rinv[iufp]
        ρ_r = apert.ρ_r[iufp]
        I1_ζ = I1 - I1_ξ - I1_η
        J_ζ = J - J_ξ - J_η
        K_ζ = K - K_ξ - K_η
      
        # Compute vector from each vertex to centroid of match triangle
        # (divided by 2) as in Eq. (7.53):
        ρc2 =  0.5 * [rmc - rm[i] for i in 1:3]
      
        # Loop over the face pairs this equivalence class:
        for (i0, ifmifs) in enumerate(rwgdat.ufp2fp[iufp])
            rowcol = i2s[ifmifs]
            ifm, ifs = rowcol[1], rowcol[2] # indices of match and source triangles
            is = @view apert.fe[:,ifs]      # Obtain the three edges of source tri.
            i_m = @view apert.fe[:,ifm]      # Obtain the three edges of the match tri.
            # Loop over each edge of the source triangle.
            for isl in 1:3 # isl is local edge index, is[isl] is global edge index.
                sbf = ebf[is[isl]]  # Source basis function for this edge.
                sbf == 0 && continue  # no basis function is defined
                if bff[1,sbf] == ifs # "Plus" source triangle:
                    source_flag = floquet_factor[eci[is[isl]]]
                elseif bff[2,sbf] == ifs # "Minus" source triangle:
                    source_flag = -floquet_factor[eci[is[isl]]]
                else
                    error("Impossible situation in source edge loop")
                end
                F_source_flag = -ϵ₀ / π * ϵ̄ * source_flag
                Ψ_source_flag = 2im / (π*ω*μ₀*μ̃) * source_flag 
          
                # Compute singular contribution for this edge (the middle term in 
                # square brackets in Equation (7-57)) using (B-2):
                Fsing = ρ_r  +  u * rinv * (rmc - rs[isl]) 
        
                # Compute Eq. (7-26) (but correct sign is carried in source flags):
                I1_i = rs[1] * I1_ξ + rs[2] * I1_η + rs[3] * I1_ζ - rs[isl] * I1
                J_i = rs[1] * J_ξ + rs[2] * J_η + rs[3] * J_ζ - rs[isl] * J
                K_i = rs[1] * K_ξ + rs[2] * K_η + rs[3] * K_ζ - rs[isl] * K
          
                # Compute Equation (7-57):
                F_i = F_source_flag * (I1fact * I1_i + Fsing + u * J_i + KFfact * K_i)
          
                # Compute Equation (7-60):
                Ψ_i = Ψ_source_flag * (I2fact * I2 + u * (rinv  + J) + KPfact * K)
          
                # Now loop over each edge of the match triangle:
                for iml in 1:3 # iml is local edge index, im[iml] is global index.
                    self_edge = i_m[iml] == is[isl] # Match edge = Source edge?
                    mbf = ebf[i_m[iml]]  # Match basis function for this edge.
                    mbf == 0 && continue # no basis function defined for this edge
                    if bff[1,mbf] == ifm  # "Plus" match triangle.
                        match_flag = conj(floquet_factor[eci[i_m[iml]]])
                    elseif bff[2,mbf] == ifm # "Minus" match tri.
                        match_flag = -conj(floquet_factor[eci[i_m[iml]]])
                    else
                        error("Impossible situation in match edge loop")
                    end
            
                    # Compute one of the dot products in Eq (7-53) apart from sign:
                    dotprod = ρc2[iml] ⋅ F_i
            
                    # Add contribution to the admittance matrix#
                    ymat[mbf,sbf] += match_flag * (-im*ω*dotprod - Ψ_i)
                end
            end # loop over source edges
        end
    end # loop over face pairs
    t2 = time_ns()
    tsec = round((t2-t1)/1e9; digits=5)
    @info "$tsec seconds to fill $(size(ymat,1)) × $(size(ymat,2)) matrix entries"

    return ymat
end      




end # module




