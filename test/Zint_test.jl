using PSSFSS
using PSSFSS.Elements: s₁s₂2β₁β₂ 
using PSSFSS.Zint: filljk!, vtxcrd
using PSSFSS.RWG: setup_rwg
#using StaticArrays
using Cubature
using LinearAlgebra: ⋅, norm
using Test


"""
    duffy(f,reltol=1e-10)

Perform the integral  ``\\int_0^1 \\int_0^{1-y} f(x,y) \\, dx \\, dy`` via a call to `pcubature` of the 
`Cubature` package using a Duffy transform.
"""
function duffy(f,reltol=1e-10)
    a = [0.,0.]
    b = [1., 1.]
    hcubature(x -> (1-x[2]) * f([(1-x[2])*x[1],x[2]]), a, b; reltol)
end

"""
    duffy2(f,reltol=1e-10)

Perform the integral  ``\\int_0^1 \\int_0^{1-y} f(x,y) \\, dx \\, dy`` via a call to `hcubature` of the 
`Cubature` package using a Duffy transform.  Here `f(ξ,η)` returns a 2-vector.
"""
function duffy2(f,reltol=1e-10)
    a = [0.,0.]
    b = [1., 1.]
    hcubature(2, (x,v) -> v[:] = (1-x[2]) * f([(1-x[2])*x[1],x[2]]), a, b; reltol)
end


function fρ_r(ξ,η,rm,rs)
    r = rs[1] + ξ*(rs[2]-rs[1]) + η*(rs[3]-rs[1]) - rm
    return r/norm(r)
end

function frinv(ξ,η,rm,rs)
    r = rs[1] + ξ*(rs[2]-rs[1]) + η*(rs[3]-rs[1]) - rm
    return 1/norm(r)
end

function twice_area(r1,r2,r3)
    r21 = r2-r1
    r31 = r3-r1
    return abs(r21[1]*r31[2] - r21[2]*r31[1])
end
twice_area(rvec) = twice_area(rvec[1],rvec[2],rvec[3])

sheet = polyring(ntri=100, sides=32, a=[0.4957], b=[0.5557],
                 s1=[1.1414, 0.0], s2=[0.5707, 0.9885], units=inch)
rwgdata = setup_rwg(sheet)
nface = size(sheet.fv,2)
nufp = rwgdata.nufp
i2s = CartesianIndices((nface,nface))
ψ₁, ψ₂ = 0., 0.  # Incremental phase shifts (rad)
s₁, s₂ = sheet.s₁, sheet.s₂ # Lattice vectors (mil)
β₁, β₂ = s₁s₂2β₁β₂(s₁,s₂) # Reciprocal lattice vectors (1/mil)
u = 0.4 * max(norm(β₁), norm(β₂)) #  smoothing parameter (1/mil)
sheet.u = u
sheet.ψ₁ = ψ₁
sheet.ψ₂ = ψ₂
extract = true
filljk!(sheet, rwgdata, extract)
conv = 1e-7


@testset "filljk! Test" begin
    for iufp = [1,2,3, nufp-3, nufp-2, nufp-1, nufp]
        ifmifs = rwgdata.ufp2fp[iufp][1]  # Obtain index into face/face matrix
        obs_src = i2s[ifmifs]
        ifm, ifs = obs_src[1], obs_src[2] # face indices: match and source
        rs = vtxcrd(ifs, sheet) # Obtain coordinates of source tri's vertices
        area2 = twice_area(rs)
        rm = sum(vtxcrd(ifm, sheet))/3 # Obtain coordinates of match tri's centroid
        if ifm == ifs
            # Singularity must occur at integration limit for hcubature
            rsnew = [rm, rs[1], rs[2]]
            a1 = twice_area(rsnew)
            rinv1, err_rinv1 = duffy(ξη -> frinv(ξη[1],ξη[2],rm,rsnew))
            ρ_r1, err_ρ1 = duffy2(ξη -> fρ_r(ξη[1],ξη[2],rm,rsnew))
            rsnew = [rm, rs[1], rs[3]]
            a2 = twice_area(rsnew)
            rinv2, err_rinv2 = duffy(ξη -> frinv(ξη[1],ξη[2],rm,rsnew))
            ρ_r2, err_ρ2 = duffy2(ξη -> fρ_r(ξη[1],ξη[2],rm,[rm, rs[1], rs[3]]))
            rsnew = [rm, rs[2], rs[3]]
            a3 = twice_area(rsnew)
            rinv3, err_rinv3 = duffy(ξη -> frinv(ξη[1],ξη[2],rm,rsnew))
            ρ_r3, err_ρ3 = duffy2(ξη -> fρ_r(ξη[1],ξη[2],rm,rsnew))
            rinv_duffy = (a1*rinv1 + a2*rinv2 + a3*rinv3) / area2
            ρ_r_duffy = (a1*ρ_r1 + a2*ρ_r2 + a3*ρ_r3) / area2
        else            
            rinv_duffy, err_rinv = duffy(ξη -> frinv(ξη[1],ξη[2],rm,rs))
            ρ_r_duffy, err_ρ_r = duffy2(ξη -> fρ_r(ξη[1],ξη[2],rm,rs))
        end
        @test abs((u*sheet.rinv[iufp] - rinv_duffy) / rinv_duffy) < conv
        @test norm((sheet.ρ_r[iufp] - ρ_r_duffy) / ρ_r_duffy) < conv
    end
end

