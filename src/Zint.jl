module Zint

using StaticArrays: SVector, MVector
using ..Sheets: SV2, MV2, RWGSheet
using ..Layers: Layer
using ..RWG: RWGData
using ..PGF: jksums
using LinearAlgebra: ⋅, norm
using Statistics: mean

export filljk!, zint



zhatcross(t) = [-t[2], t[1]]
zhatcross(t::MV2) = MV2(-t[2], t[1])
zhatcross(t::SV2) = SV2(-t[2], t[1])




"""
    filljk!(metal::RWGSheet, rwgdat::RWGData, closed::Bool)
    
Compute matrices of frequency-independent integrals needed in filling 
the generalized impedance matrix.

## Arguments:

- `metal`: A variable of `RWGSheet` containing the face 
  information, unit cell incremental phase shifts, and 
  Green's function smoothing parameter.
- `rwgdat`:  A variable of type `RWGData` from which we make
  use of the arrays `ufpm` and `ufp2fp`.
- `closed`:    A logical flag which, if true, instructs this function
  to perform singularity extraction and analytic integration
  for all source triangles. If false then analytic singularity
  integration is performed ONLY for the case when source and 
  observation triangles are the same.

##  Outputs:

- `metal`:  The following fields of `metal` are modified:
  `J`, `J_ξ`, `J_η`, `K`, `K_ξ`, `K_η`.  These are matrices
  of face-pair scalar integrals defined in Equations 
  (7.22) through (7.27) of the theory documentation. Note these
  integrals are unitless.  Also modified are the fields `ρ_r` (the 
  vector of singular vector integrals defined in Equation (B.3b), 
  divided by twice the source triangle area to make them unitless, 
  and `rinv`, a vector of singular scalar integrals defined in Eq. 
  (B.3a), divided by u times twice the source triangle 
  area to make them unitless.

"""
function filljk!(metal::RWGSheet, rwgdat::RWGData, closed::Bool)
    ξ = SVector(0.33333333333333330, 0.10128650732345633, 0.79742698535308730,
                0.10128650732345633, 0.47014206410511505, 0.05971587178976989,
                0.47014206410511505)
    η = SVector(0.33333333333333333, 0.79742698535308730, 0.10128650732345633,
                0.10128650732345633, 0.05971587178976989, 0.47014206410511505,
                0.47014206410511505)
    wght = SVector(0.1125, 0.06296959027241358, 0.06296959027241358, 0.06296959027241358,
                   0.06619707639425308, 0.06619707639425308, 0.06619707639425308)

    next = (2, 3, 1) 

    nface = size(metal.fe, 2) # Number of faces in triangulated sheet.
    nedge = length(metal.e1) # Number of edges in triangulated sheet.
    nufp = rwgdat.nufp # Number of unique face pairs
    J = zeros(ComplexF64, nufp)
    J_ξ = zeros(ComplexF64, nufp)
    J_η = zeros(ComplexF64, nufp)
    K = zeros(ComplexF64, nufp)
    K_ξ = zeros(ComplexF64, nufp)
    K_η = zeros(ComplexF64, nufp)
    ρ_r = zeros(typeof(MV2(0.,0.)), nufp)
    rinv = zeros(Float64, nufp) 

    i2s = CartesianIndices((nface,nface))

    ulocal = metal.u  # Obtain smoothing parameter in units that
    #                 # are consistent with metal's length units.
    us1 = ulocal * metal.s₁
    us2 = ulocal * metal.s₂

    Threads.@threads for iufp ∈ 1:rwgdat.nufp  # Loop over each unique face pair
        ifmifs = rwgdat.ufp2fp[iufp][1]  # Obtain index into face/face matrix
        rowcol = i2s[ifmifs]
        ifm, ifs = rowcol[1], rowcol[2] # indices of match and source triangles
        is = @view metal.fe[:,ifs]     # Obtain the three edges of the source triangle.
        
        #vtxcrd!(rs, ifs, metal) # Obtain coordinates of source tri's vertices
        rs = vtxcrd(ifs, metal) # Obtain coordinates of source tri's vertices
        # Calculate twice the signed area of the source triangle
        rs32 = rs[3] - rs[2]
        rs12 = rs[1] - rs[2]
        area2 = zhatcross(rs32) ⋅ rs12
        unz = 1.0   # Sign of unit surface normal (in z direction)
        if area2 < 0
            area2 = abs(area2)
            unz = -1.0
        end

        rs31 = rs[3] - rs[1]
        rs21 = -rs12
      
        self_tri = ifm == ifs
        clsflg = self_tri || closed  # Extract if self or if always
        im = @view metal.fe[:,ifm]   # Obtain the three edges of the observation triangle
        rm = vtxcrd(ifm, metal) # observation triangle vertex coordinates
        rmc = mean(rm) # observation face centroid (meters)
        
        for i ∈ 1:length(ξ)   # Loop for numerical integration
            rt = rs[1] + rs21 * ξ[i] + rs31 * η[i] # Source point
            uρ00 = ulocal * (rmc - rt)
            Jsum, Ksum = jksums(uρ00, metal.ψ₁, metal.ψ₂, us1, us2, clsflg) # spatial sums
            # Accumulate results of numerical integration:
            Jsumwght = Jsum * wght[i]
            J_ξ[iufp] += Jsumwght * ξ[i]
            J_η[iufp] += Jsumwght * η[i]
            J[iufp] += Jsumwght
            Ksumwght = Ksum * wght[i]
            K_ξ[iufp] += Ksumwght * ξ[i]
            K_η[iufp] += Ksumwght * η[i]
            K[iufp] += Ksumwght
        end
        if clsflg
            # compute closed-form of singular integrals from Eqs (B.3)
            for i ∈ 1:3  # Loop over sides of source triangle
                ip1 = next[i]  # Next edge in cyclic list.
                # Find components of l, unit tangent vector to edge #i
                lvec = rs[ip1] - rs[i]
                lvec /= norm(lvec)  # Make into unit vector.
                
                # Find unit outward normal to edge i, lying in plane z = 0
                uvec = -unz * zhatcross(lvec) # so (ux,uy) = (ly * unz, -lx * unz)
                
                # Compute signed perp. distance from observation point to edge #i
                p0 = uvec ⋅ (rs[i] - rmc)
                p0sq = p0 * p0
                
                # Compute lplus and lminus as defined in the reference
                lp = lvec ⋅ (rs[ip1] - rmc) 
                lm = lvec ⋅ (rs[i] - rmc) 
                ledge = lp - lm   # Length of the edge
                
                pplus = sqrt(p0sq + lp * lp)
                pminus = sqrt(p0sq + lm * lm)
                
                # Check for special case of observation point on extension of edge #i
                if abs(p0/ledge) < 1e-4
                    p0 = 0.0
                    factor = lp * pplus  -  lm * pminus
                else
                    factor = p0 * log((pplus + lp) / (pminus + lm))
                    rinv[iufp] += factor
                    factor = factor * p0 + pplus * lp - pminus * lm
                end
                ρ_r[iufp] += uvec * factor
            end
            ρ_r[iufp] *= 0.5
            # rinv is now the quantity defined in (B.3a)
            # ρ_r is now the quantity defined in (B.3b)
            ρ_r[iufp] /= area2  # Make it unitless
            rinv[iufp] /= (area2 * ulocal) # Make it unitless
        end
      
    end

    metal.J = J 
    metal.J_ξ = J_ξ 
    metal.J_η = J_η 
    metal.K = K 
    metal.K_ξ = K_ξ 
    metal.K_η = K_η 
    metal.ρ_r = ρ_r 
    metal.rinv = rinv 

    return nothing
end  

#=
"""
 Return the coordinates (in local units) of the triangle vertices for face iface.
"""
@inline function vtxcrd!(ρvec, iface, metal)
    vi = @view metal.fv[:,iface] # Vertex indices
    ρvec[:] = metal.ρ[vi]
end
=#


"""
 Return the coordinates (in local units) of the triangle vertices for face iface.
"""
@inline function vtxcrd(iface, metal)
    vi = @view metal.fv[:,iface] # Vertex indices
    @view metal.ρ[vi]
end


"""
    zint(Σm1_func, Σm2_func, rs, rmc)  --> (I1, I1_ξ, I1_η, I2)

Compute frequency-dependent integrals needed to fill the generalized impedance matrix.

## Arguments

- `Σm1_func`, `Σm2_func`:  Functions that evaluate the modal series.
- `rs`:  Vector of length three containing the source triangle vertices in meters.
         Each vertex is a 2-vector such as a StaticArrays SV2.
- `rmc`  A 2-vector containing the field observation (match) point in meters.

## Return Values:

- `I1`, `I1_ξ`, `I1_η`, `I2`: Complex, frequency-dependent, spectral integrals defined
                              by Eqs (7.22a), (7.27), and (7.32a).  Their units are (1/m).
"""
function zint(Σm1_func::Function, Σm2_func::Function, rs, rmc)
    ξ = SVector(0.33333333333333330, 0.10128650732345633, 0.79742698535308730,
                0.10128650732345633, 0.47014206410511505, 0.05971587178976989,
                0.47014206410511505)
    η = SVector(0.33333333333333333, 0.79742698535308730, 0.10128650732345633,
                0.10128650732345633, 0.05971587178976989, 0.47014206410511505,
                0.47014206410511505)
    wght = SVector(0.1125, 0.06296959027241358, 0.06296959027241358, 0.06296959027241358,
                   0.06619707639425308, 0.06619707639425308, 0.06619707639425308)

    I1_ξ = zero(ComplexF64)
    I1_η = zero(ComplexF64)
    I1 = zero(ComplexF64)
    I2 = zero(ComplexF64)

    for i ∈ 1:length(ξ)
        rt = rs[1] + (rs[2] - rs[1]) * ξ[i]  +  (rs[3] - rs[1]) * η[i] # Source point
        ρdif = rmc - rt 
        sig1 = Σm1_func(rhodif)
        sig2 = Σm2_func(rhodif)
        sig1wght = sig1 * wght[i]
        I1_ξ += sig1wght * ξ[i]
        I1_η += sig1wght * η[i]
        I1 += sig1wght
        I2 += sig2 * wght[i]
    end
    return   (I1, I1_ξ, I1_η, I2)
end  # function

end # module






