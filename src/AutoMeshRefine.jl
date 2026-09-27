module AutoMeshRefine

using PSSFSS.Sheets: facecount, RWGSheet
using PSSFSS.RWG: RWGData
using PSSFSS.Zint: vtxcrd
using Unitful: @u_str, ustrip
using OffsetArrays: OffsetArray
using StaticArrays: SVector

"""
    tricharge!(charges, currents, rwg, sheet, ω)

Compute the charge density on each triangle face of the sheet.

## Input Arguments

* `charges`: A pre-allocated, complex vector of length `ntri`, where `ntri == facecount(sheet)`.
  The computed charge densities will be saved into this vector after the function completes execution,
  with units of Coulomb per square meter (i.e. [C/m²]).
* `currents`: A vector of complex current coefficients of length `size(rwg.bfe, 2)`, i.e. the number of
  basis functions defined in `rwg`.
* `rwg`: An instance of `PSSFSS.RWGData` defining the Rao-Wilton-Glisson basis function structure.
* `sheet`: An instance of `PSSFSS.RWGSheet` defining the FSS/PSS sheet triangulation. It is assumed herein
  that `sheet.ψ₁` and `sheet.ψ₂` are initialized properly.
* `ω`: The radian frequency (rad/sec).
"""
function tricharge!(
    charges::AbstractVector{<:Complex},
    currents::AbstractVector{<:Complex},
    rwg::RWGData,
    sheet::RWGSheet,
    ω::AbstractFloat)

    ntri = facecount(sheet)
    ntri == length(charges) || error("charges length not compatible with sheet data")
    nbf = size(rwg.bfe, 2)
    nbf == length(currents) || error("currents vector length not compatible with RWG data")

    ψ₁, ψ₂ = sheet.ψ₁, sheet.ψ₂
    units_per_meter = ustrip(Float64, sheet.units, 1u"m")
    # floquet_factor is indexed into using values in rwg.eci:
    floquet_factor = OffsetArray((im / ω) * SVector(1.0, 1.0, cis(-ψ₁), 1.0, cis(-ψ₂)), 0:4)

    charges .= zero(eltype(charges))
    for face in 1:ntri
        rs = vtxcrd(face, sheet) ./ units_per_meter
        rs32 = rs[3] - rs[2]; rs12 = rs[1] - rs[2]
        areainv = inv(0.5 * abs(rs32[1] * rs12[2] - rs32[2] * rs12[1]))
        edges = @view sheet.fe[:, face]
        for edge in edges
            bf = rwg.ebf[edge] # basis function index
            iszero(bf) && continue
            ffactor = floquet_factor[rwg.eci[edge]]
            poscharge = areainv * ffactor * currents[bf]
            if face == rwg.bff[1, bf]
                charges[face] += poscharge
            elseif face == rwg.bff[2, bf]
                charges[face] -= poscharge
            else
                error("Impossible situation for edge $edge, face $face, basis function $bf")
            end
        end
    end
end



"""
    tricharge2!(charges, currents, rwg, sheet, ω)

Compute the charge density on each triangle face of the sheet (alternate method).

## Input Arguments

* `charges`: A pre-allocated, complex vector of length `ntri`, where `ntri == facecount(sheet)`.
  The computed charge densities will be saved into this vector after the function completes execution,
  with units of Coulomb per square meter (i.e. [C/m²]).
* `currents`: A vector of complex current coefficients of length `size(rwg.bfe, 2)`, i.e. the number of
  basis functions defined in `rwg`.
* `rwg`: An instance of `PSSFSS.RWGData` defining the Rao-Wilton-Glisson basis function structure.
* `sheet`: An instance of `PSSFSS.RWGSheet` defining the FSS/PSS sheet triangulation. It is assumed herein
  that `sheet.ψ₁` and `sheet.ψ₂` are initialized properly.
* `ω`: The radian frequency (rad/sec).
"""
function tricharge2!(
    charges::AbstractVector{<:Complex},
    currents::AbstractVector{<:Complex},
    rwg::RWGData,
    sheet::RWGSheet,
    ω::AbstractFloat)

    ntri = facecount(sheet)
    ntri == length(charges) || error("charges length not compatible with sheet data")
    nbf = size(rwg.bfe, 2)
    nbf == length(currents) || error("currents vector length not compatible with RWG data")

    ψ₁, ψ₂ = sheet.ψ₁, sheet.ψ₂
    units_per_meter = ustrip(Float64, sheet.units, 1u"m")
    # floquet_factor is indexed into using values in rwg.eci:
    floquet_factor = OffsetArray((im / ω) * SVector(1.0, 1.0, cis(-ψ₁), 1.0, cis(-ψ₂)), 0:4)

    charges .= zero(eltype(charges))
    for bf in 1:nbf
        for (fsign, edge, face) in zip((1, -1), (@view rwg.bfe[:, bf]), (@view rwg.bff[:, bf]))
            rs = vtxcrd(face, sheet) ./ units_per_meter
            rs32 = rs[3] - rs[2]; rs12 = rs[1] - rs[2]
            areainv = inv(0.5 * abs(rs32[1] * rs12[2] - rs32[2] * rs12[1]))
            ffactor = floquet_factor[rwg.eci[edge]]
            charge = fsign * areainv * ffactor * currents[bf]
            charges[face] += charge
        end
    end
end

end # module