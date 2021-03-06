"""
# Module RWG 
This module defines the modified Rao, Wilton, and Glisson triangle 
subdomain basis function derived data type.  The definition of a basis 
function has been generalized slightly to allow the "plus" and "minus" 
triangles to be nonadjacent, and to identify those half-basis functions 
whose defining edge lies at the ξ=1 or η=1 boundaries of the unit cell.
These modifications allow one to handle the analysis of a structure in a 
periodic unit cell, as required in phased array and frequency selective 
surface formulations.  Reference: P. S. Simon, "Modified RWG basis functions 
for analysis of periodic structures," *2002 IEEE MTT-S International Microwave 
Symposium Digest* (Cat. No.02CH37278), Seattle, WA, USA, 2002, pp. 2029-2032 
vol. 3, doi: 10.1109/MWSYM.2002.1012266.
"""
module RWG
export RWGData, setup_rwg, edge_current_unit_vector, rwgbfft!


using ..Sheets: RWGSheet, SV2
using StaticArrays: SVector, MVector, SArray
using LinearAlgebra: ⋅, norm
using NearestNeighbors: KDTree, inrange
using OffsetArrays
using Statistics: mean
using ..PSSFSSLen

mutable struct RWGData
    #  Basis function edge indices.  The value in bfe[1,i] is index of the
    #  edge of the "plus" triangle associated with basis function i.
    #  The value in bfe[2,i] is the index of the edge of the "minus" triangle
    #  associated with basis function i. In most cases, these
    #  two values (edges) will be identical.
    bfe::Array{Int,2}
    
    #  Basis function face indices.  The value in bff[1,i] is index of the
    #  "plus" triangle face associated with basis function i.
    #  The value in bff[2,i] is the index of the "minus" triangle face
    #  associated with basis function i. 
    bff::Array{Int,2}
    
    #  Edge basis function map.
    #  ebf[i], if nonzero, is the index of the basis function associated with 
    #  edge i.
    ebf::Array{Int,1}
    
    #  Edge cell index.
    #  eci[i] takes on the values 0, 1, 2, 3, or 4.   The values have the 
    #  following meanings:
    #  0  The edge does not lie along a unit cell boundary.
    #  1  The edge lies along the ξ=0 unit cell boundary.
    #  2  The edge lies along the ξ=1 unit cell boundary.
    #  3  The edge lies along the η=0 unit cell boundary.
    #  4  The edge lies along the η=1 unit cell boundary.
    eci::Array{Int,1}
    
    #  Unique face pair matrix.
    #  ufpm[i,j] contains the unique face pair index for observation face i 
    #  with respect to source face j.  Two sets of face pairs are considered 
    #  to be equivalent if they can be made to superimpose via a rigid 
    #  translation, and the nodes in the corresponding triangles of each pair
    #  are numbered in the same order. Since the periodic Green's functions are 
    #  translationally invariant, the integrals involving equivalent face pairs
    #  will have identical values and thus need be computed only a single time.
    ufpm::Array{Int,2}
    
    #  ufp2fp[i] contains the vector of face pair indices for equivalence class i. The 
    #  face pair index uses column major ordering to enumerate the elements 
    #  of a matrix of order Nface × Nface.
    ufp2fp::Array{Array{Int,1},1}

    zorymat::Array{ComplexF64,2} # MoM matrix
    rhs::Vector{ComplexF64}  # MoM right-hand side 
    bfftstore::Array{SArray{Tuple{2},ComplexF64,1,2}, 3} # Basis function Fourier Transforms

    nufp::Int  # Number of unique face pairs.
end # mutable struct
  
  
"""
    setup_rwg(sheet::RWGSheet)::RWGdata

This function accepts the sheet geometry data structure as created
by the function `get_sheet_data` and creates a instance of `RWGdata` as the
function return value.  When `sheet.fufp` is `true`, it directs this 
function to search for unique face pairs.  This search can be time consuming 
and thus the default action, when `sheet.fufp` is .false., 
is to skip the search.  The tradeoff is the greater time needed to fill the 
interaction matrix when all face pairs are considered unique.
"""
function setup_rwg(sheet::RWGSheet, leafsize::Int=9)::RWGData
    tol = 1.e-5 # Comparison tolerance
    ieη0 = Int[] # List edges at η=0
    ieη1 = Int[] # List edges at η=1
    ieξ0 = Int[] # List edges at ξ=0
    ieξ1 = Int[] # List edges at ξ=1

    nface = size(sheet.fe, 2) # Number of faces in triangulated sheet.
    nedge = length(sheet.e1) # Number of edges in triangulated sheet.

    eci = zeros(Int, nedge)


    # Count the number of basis functions:
    nbf = 0
    # To begin with, we locate any edge that is adjacent to a pair of faces...
    for fp in 1:nface, fm in fp+1:nface  # Loop over "plus" and "minus" triangles
        for gep in @view sheet.fe[:,fp]  # Global index of "plus" edge.
            for gem in @view sheet.fe[:,fm] # Global index of "minus" edge.
                gep == gem && (nbf += 1)    # Faces share common edge
            end
        end
    end
    # We must now add to this count the number of edges that are located at 
    # the ξ=0  or the η=0 boundaries of the unit cell.  In the process of 
    # locating such edges, we will also set the correct values in the eci, 
    # ieξ0, ieξ1, ieη0, and  ieη1 arrays.
    if sheet.ξη_check
        for gep in 1:nedge  # Loop over each edge of the structure.
            # Obtain the ξ and η values for each vertex of edge ie:
            (ξinit, ξterm, ηinit, ηterm) = get_ie_ξη(gep, sheet)
            if abs(ξinit) < tol  &&  abs(ξterm) < tol
                # Edge is located on the ξ=0 boundary of unit cell.
                nbf += 1         # Bump count of # basis functions.
                eci[gep] = 1     # Store edge code.
                push!(ieξ0, gep) # Update list of edges at ξ = 0.
            elseif abs(ξinit-1) < tol  && abs(ξterm-1) < tol
                # Edge is located on the ξ=1 boundary of unit cell.
                eci[gep] = 2
                push!(ieξ1, gep)  # Update list of edges at ξ = 1.
            elseif abs(ηinit) < tol  && abs(ηterm) < tol
                # Edge is located on the η=0 boundary of unit cell.
                nbf += 1
                eci[gep] = 3 
                push!(ieη0, gep)  # Update list of edges at η = 0.
            elseif abs(ηinit-1) < tol && abs(ηterm-1) < tol
                # Edge is located on the η=1 boundary of unit cell.
                eci[gep] = 4          # Store edge code.
                push!(ieη1, gep)  # Update list of edges at η = 1.
            end
        end # for
    end # if
    length(ieξ0) == length(ieξ1) || error("Inconsistent # edges at ξ=0 and ξ=1")
    length(ieη0) == length(ieη1) || error("Inconsistent # edges at η=0 and η=1")

    bfe = zeros(Int, 2,nbf)
    bff = zeros(Int, 2,nbf)
    ebf = zeros(Int, nedge)

    # Loop over pairs of triangles.  Note that the order of the loops dictates 
    # that the face with lower index will be the "plus" face.
    i = 0  # Basis function index
    for fp in 1:nface, fm = fp+1:nface            # Loop over "plus" and "minus" triangles.
        for gep in view(sheet.fe, :, fp), gem in view(sheet.fe, :, fm) # Global edge index for "+" and "-" triangles.
            if gep == gem   # Faces share common edge:
                i += 1        # Bump basis function index.
                bfe[1,i] = gep # "Plus" triangle edge index for i'th basis funct
                bfe[2,i] = gem # "Minus" triangle edge index for i'th basis funct
                bff[1,i] = fp  # "Plus" triangle face index for i'th basis funct
                bff[2,i] = fm  # "Minus" triangle face index for i'th basis funct
                ebf[gep] = i
            end
        end
    end

    #  Now treat the basis functions located at unit cell boundaries.
    for ge0 in ieη0  # Loop over edges at η=0 unit cell boundary.
        # Search for the triangle face containing edge ge0:
        face0 = 0  # Initialize index of face containing edge at η=0.
        for j0 in 1:nface
            if ge0 in @view sheet.fe[:,j0]
                face0 = j0  # Save face index.
                break
            end
        end
        face0 == 0 && error("Unable to locate face containing η=0 edge $ge0")

        # Obtain the initial and terminal ξ values of edge at η = 0:
        (ξinit, ξterm, _, _) = get_ie_ξη(ge0, sheet) 
        ξmin0, ξmax0 = extrema((ξinit, ξterm))

        # Search for the corresponding edge at η = 1:
        ge1 = 0   # Initialize global index of corresponding edge at η = 1.
        for ge in ieη1  # Loop over each edge at the η=1 unit cell boundary.
            # Obtain the initial and terminal ξ values:
            (ξinit, ξterm, _, _) = get_ie_ξη(ge, sheet) 
            ξmin1, ξmax1 = extrema((ξinit, ξterm))
            # Search for edges occupying the same ξ interval:
            if abs(ξmin1-ξmin0) < tol && abs(ξmax1-ξmax0) < tol
                ge1 = ge     # Save global edge index for edge at η = 1.
                break
            end
        end
        ge1 == 0 && error("Unable to locate edge at η=1 corresponding to edge $ge0")
        # Search for the triangle face containing edge ge1:
        face1 = 0  # Initialize index of face containing edge at η=1.
        for j1 in 1:nface
            if ge1 in @view sheet.fe[:,j1]
                face1 = j1  # Save face index.
                break        # Jump out of search loop.
            end
        end
        face1 == 0 && error("Unable to locate face containing η=1 edge $ge1")
        i += 1      # Bump basis function index.
        bfe[1,i] = ge0 # "Plus" triangle edge index for i'th basis function.
        bfe[2,i] = ge1 # "Minus" triangle edge index for i'th basis function.
        bff[1,i] = face0  # "Plus" triangle face index for i'th basis function.
        bff[2,i] = face1  # "Minus" triangle face index for i'th basis function.
        ebf[ge0] = i
        ebf[ge1] = i
    end # for edges at η=0
    #
    for ge0 in ieξ0  # Loop over edges at ξ=0 unit cell boundary.
        # Search for the triangle face containing edge ge0:
        face0 = 0  # Initialize index of face containing edge at η=0.
        for j0 = 1:nface
            if ge0 in @view sheet.fe[:,j0]
                face0 = j0  # Save face index.
                break
            end
        end
        face0 == 0 && error("Unable to locate face containing ξ=0 edge $ge0")
        # Obtain the initial and terminal η values of edge at ξ = 0:
        (_, _, ηinit, ηterm) = get_ie_ξη(ge0, sheet)
        ηmin0, ηmax0 = extrema((ηinit, ηterm))
        # Search for the corresponding edge at ξ = 1:
        ge1 = 0   # Initialize global index of corresponding edge at ξ = 1.
        for ge in ieξ1  # Loop over each edge at the ξ=1 unit cell boundary.
            # Obtain the initial and terminal η values:
            (_, _, ηinit, ηterm) = get_ie_ξη(ge, sheet) 
            ηmin1, ηmax1 = extrema((ηinit, ηterm))
            # Search for edges occupying the same η interval:
            if abs(ηmin1-ηmin0) < tol && abs(ηmax1-ηmax0) < tol
                ge1 = ge     # Save global edge index for edge at ξ = 1.
                break
            end
        end
        ge1 == 0 && error("Unable to locate edge at ξ=1 corresponding to edge $ge0")
        # Search for the triangle face containing edge ge1:
        face1 = 0  # Initialize index of face containing edge at η=1.
        for j1 = 1:nface
            if ge1 in @view sheet.fe[:,j1]
                face1 = j1  # Save face index.
                break
            end
        end
        face1 == 0 && error("Unable to locate face containing ξ=1 edge $ge1")
        i += 1      # Bump basis function index.
        bfe[1,i] = ge0 # "Plus" triangle edge index for i'th basis function.
        bfe[2,i] = ge1 # "Minus" triangle edge index for i'th basis function.
        bff[1,i] = face0  # "Plus" triangle face index for i'th basis function.
        bff[2,i] = face1  # "Minus" triangle face index for i'th basis function.
        ebf[ge0] = i
        ebf[ge1] = i
    end # ξ0_loop

    i == nbf || error("Inconsistent number of basis functions")
    zorymat = zeros(ComplexF64, nbf, nbf)
    rhs = zeros(ComplexF64, nbf)
    bfftstore = zeros(SArray{Tuple{2},ComplexF64,1,2}, 0,0,0)
    nufp =  nface*nface 
    ufpm = reshape(collect(1:nufp), (nface,nface))
    if !sheet.fufp  # Don't search for unique face pairs. Assume all are unique.
        ufp2fp = [ [i] for i in 1:nufp]
        return RWGData(bfe, bff, ebf, eci, ufpm, ufp2fp, zorymat, rhs, bfftstore, nufp)
    end
        
    
    #  The remaining code in this function sets up the
    #  two arrays ufpm and ufp2fp.  These are used to identify
    #  redundant face/pairs.  ufpm(i) contains the unique face/pair
    #  index (or the equivalence class (E.C.) index) of face/pair i.  
    #  The face pairs are numbered sequentially in the same manner as 
    #  the elements of a square matrix of dimension nface (column major 
    #  ordering). The row index is used for the observation face number and 
    #  the column index is the source face number.  The array ufp2fp is 
    #  of length nufp (the number of unique face pairs, or equivalence 
    #  classes). Each element of ufp2fp contains a pointer to an allocated 
    #  array of face/pair indices.  Row i contains a list of the face/pairs 
    #  that are members of equivalence class i.  Two face pairs belong to 
    #  the same equivalence class if the source triangles can be overlaid
    #  using a rigid translation, the source triangle nodes are listed
    #  in the same order, and if the observation point (centroid 
    #  of the match triangle) is in the same relative position wrt the 
    #  source triangle.  This version uses the NearestNeighbors package:

      #  Allocate the unique face pairs matrix and some scratch arrays.
    ufpm = zeros(Int, (nface,nface))
    
    r1 = @view sheet.ρ[sheet.fv[1,:]]
    r2 = @view sheet.ρ[sheet.fv[2,:]]
    r3 = @view sheet.ρ[sheet.fv[3,:]]
    centroid = [(r1[n] + r2[n] + r3[n]) / 3 for n in eachindex(r1)]
    data = Array{SVector{6,Float64}}(undef, nface^2)
    #  Characterize each face pair:
    mn = 0
    for n in 1:nface, m in 1:nface
        mn += 1 # Bump single face/face index.
        # Calculate test vectors for face pair (m,n):
        rmn1 = centroid[m] - r1[n]
        rmn2 = centroid[m] - r2[n]
        rmn3 = centroid[m] - r3[n]
        data[mn] = SVector{6,Float64}(rmn1[1], rmn1[2], rmn2[1], rmn2[2], rmn3[1], rmn3[2])
    end

    kdtree = KDTree(data, leafsize=leafsize)
    ufp2fp = Vector{Int}[]
    mn = 0  # Initialize face/pair index.
    nufp = 0
    found = fill(false, nface^2)
    r = 1e-6 * norm(centroid[1] - centroid[2])
    for n in 1:nface, m in 1:nface
        mn += 1  # Bump global face/pair counter.
        found[mn] && continue
        nufp += 1
        idxs = inrange(kdtree, data[mn], r, true)
        found[idxs] .= true
        push!(ufp2fp, idxs)
        ufpm[idxs] .= nufp
    end
    nufp == length(ufp2fp) || error("Miscount of nufp")
    
    return RWGData(bfe, bff, ebf, eci, ufpm, ufp2fp, zorymat, rhs, bfftstore, nufp)    

end # function setup_rwg



zhatcross(t) = SV2([-t[2], t[1]])
zhatcross(t::SV2) = SV2(-t[2], t[1])



"""
    edge_current_unit_vector(ie::Integer, rwgdat::RWGData, metal::RWGSheet)::SV2

Evaluate a unit vector u = SV2([ux,uy]) in the positive reference direction for 
the basis function associated with edge ie of the triangulated sheet.
"""
function edge_current_unit_vector(ie::Integer, rwgdat::RWGData, metal::RWGSheet)::SV2
    bf = rwgdat.ebf[ie]    #  Basis function index associated with edge ie
    bf == 0 && error("No basis function for edge $ie")
    # Begin by assuming the normal is parallel to \zhat \cross (\r_2 - \r_1):
    n1 = metal.e1[ie]   # Initial node of edge ie.
    n2 = metal.e2[ie]   # Terminal node of edge ie.
    ρ21 = metal.ρ[n2] - metal.ρ[n1]
    d = norm(ρ21)
    u = zhatcross(ρ21/d)
    # Now check that dot product of assumed unit vector with \vecrho^+ evaluated
    # at one of the edge nodes is positive:
    f = rwgdat.bff[1,bf]# Obtain the "plus" face adjacent to edge ie.
    # Find the free vertex of this face:
    nfree = 0 # Establish scope
    for i in 1:3
      nfree = metal.fv[i,f]
      nfree != n1  &&  nfree != n2 && break
    end
    ρ = metal.ρ[n1] - metal.ρ[nfree]
    if ρ ⋅ u < 0 
      return -u
    end
    return u
end 


        
"""
    get_ie_ξη(ie::Int, sheet::RWGSheet) -> (ξinit, ξterm, ηinit, ηterm)

Evaluate the ξ and η coordinates of the initial and terminal vertices 
of a given edge. The position vector of a point is represented as 
ξ*s1 + η*s2, where s1 and s2 are the direct lattice vectors.
 
## Arguments:

- `ie`:  Global edge index.
- `sheet`: The RWGSheet object containing the triangulation info.

## Return value:

A 4-tuple containing

- `ξinit`: The ξ value of the initial vertex of edge `ie`.
- `ξterm`: The ξ value of the terminal vertex of edge `ie`.
- `η_init`: The η value of the initial vertex of edge ie.
- `η_term`: The η value of the terminal vertex of edge ie.

"""
@inline function get_ie_ξη(ie::Int, sheet::RWGSheet)
    node_init = sheet.e1[ie]  
    node_term = sheet.e2[ie]  
    ρinit = sheet.ρ[node_init]
    ρterm = sheet.ρ[node_term]
    ξinit = (sheet.β₁ ⋅ ρinit) / (2π)
    ηinit = (sheet.β₂ ⋅ ρinit) / (2π)
    ξterm = (sheet.β₁ ⋅ ρterm) / (2π)
    ηterm = (sheet.β₂ ⋅ ρterm) / (2π)
    return (ξinit, ξterm, ηinit, ηterm)
end


"""
    rwgbfft!(ft, rwgdat::RWGData, sheet::RWGSheet, k::AbstractVector, ψ₁::Real, ψ₂::Real) -> ft::Vector
    
Compute the 2D fourier transform (FT) of the set of modified Rao-Wilton-Glisson basis functions 
defined in `rwgdat` and `sheet`, evaluated at the transform variable `k`. `ψ₁` and `ψ₂` are the 
unit cell incremental phase shifts in radians.

## Arguments:

- `ft`: A vector of length `nbf` whose elements are 2-vectors, where `nbf` is the number of 
    basis functions.  The contents of `ft` will be overwritten with the basis function FTs.
- `rwgdat`: Contains the basis function definitions.
- `sheet`: The contains the triangulation info.
- `k`: A 2-vector containing the transform variable value at which the fourier transforms are to be 
    evaluated. Units are (1/meter).
- `ψ₁`, `ψ₂`:  Unit cell incremental phase shifts in units of radians.

## Return value:

A length `nbf` complex vector containing the Fourier transforms of the basis functions.  Units are meters^2. 
`nbf` is the number of basis functions.

## Reference
Kim McInturff and Peter S. Simon, "The Fourier transform of linearly varying functions with polygonal
support," IEEE Trans. Antennas Propagat., Vol. 39, no. 9, Sept. 1991, pp. 1441-1443.

"""
function rwgbfft!(ft, rwgdat::RWGData, sheet::RWGSheet, k::AbstractVector, ψ₁::Real, ψ₂::Real)
    next = SVector{3, Int}(2,3,1)
    nbf = size(rwgdat.bfe,2) # Number of basis functions
    length(ft) == nbf || error("length(ft) is $(length(ft)) instead of $nbf")
    ft .*= 0
    floquet_factor = OffsetArray(SVector{5, ComplexF64}(1, 1, cis(-ψ₁), 1, cis(-ψ₂)),   0:4)
    # Meanings:
    # floquet_factor[0] Edges not on unit cell boundary
    # floquet_factor[1] Edges at ξ=0 boundary.
    # floquet_factor[2] Edges at ξ=1 boundary.
    # floquet_factor[3] Edges at η=0 boundary.
    # floquet_factor[4] Edges at η=1 boundary.

    nface = size(sheet.fv,2) # Number of triangular faces
    kmagsq = k ⋅ k
    kmag = sqrt(kmagsq)
    one_meter = ustrip(Float64, sheet.units, 1u"m")

    lvec = zeros(SV2, 3)
    rc = zeros(SV2, 3)
    j0kl2 = MVector{3,Float64}(0.0, 0.0, 0.0)
    rtrm2 = zeros(SV2, 3)
    cphasv = MVector{3,ComplexF64}(0.0im, 0.0im, 0.0im)
    csum = MVector{2,ComplexF64}(0.0+0.0im, 0.0+0.0im)
    centroid = SVector{2,Float64}(0.0, 0.0)
    ft0 = SVector{2,ComplexF64}(0.0+0.0im, 0.0+0.0im)
    for iface in 1:nface
        r =  vtxcrd_m(iface, sheet, one_meter)
        lvec .= (r[next[i]] - r[i] for i in 1:3) # Edge vectors
        rc .= (r[i] + 0.5 * lvec[i] for i in 1:3)  # Edge centers
        if kmag * maximum(norm.(lvec)) < 1e-4 # small k
            centroid = mean(r)   #  Compute centroid coordinates.
            cphase = cis(k ⋅ centroid) #  Phase factor at centroid.
            for i in 1:3 #  Loop over three edges of the triangle
                ie = sheet.fe[i,iface] # Global index for edge opposite r(i).
                ib = rwgdat.ebf[ie]    # Global basis function index.
                ib == 0 && continue
                ft0 = 0.5 * (centroid - r[i])
                if rwgdat.bff[1,ib] == iface  # Plus triangle
                    ft[ib] += ft0 * (cphase * floquet_factor[rwgdat.eci[ie]])
                elseif rwgdat.bff[2,ib] == iface # Minus triangle
                    ft[ib] -= ft0 * (cphase * floquet_factor[rwgdat.eci[ie]])
                else
                    error("Impossible situation!")
                end
            end
        else # Case where k is not small
            kfact = 2k / kmagsq
            darea = zdotaxb(lvec[1], lvec[2]) # Twice the directed area
            denom = darea * kmagsq
            zdotlxk = [zdotaxb(l, k) for l in lvec]
            ctrm1 = [im * rci - kfact for rci in rc]
            for i in 1:3  # Loop over three edges
                dotkl2 = 0.5 * (k ⋅ lvec[i])
                j0kl2[i] = j₀(dotkl2)
                rtrm2[i] = 0.5 * zdotlxk[i] * j₁(dotkl2) * lvec[i]
                cphasv[i] = cis(k ⋅ rc[i])
            end
            #
            #  Loop over edges of triangle associated with basis functions.
            #
            for i in 1:3 
                ie = sheet.fe[i,iface] # Global index for edge opposite r[i].
                ib = rwgdat.ebf[ie]    # Global basis function index.
                ib == 0 && continue # Skip if no basis func for this edge.
                cjr = im * r[i]
                csum .= complex(0.0,0.0)
                for n in 1:3  # Perform sum over n as shown in Equation (2-11):
                    ctrm3 = (zhatcross(lvec[n]) + zdotlxk[n]*(ctrm1[n]-cjr)) * j0kl2[n]
                    csum += cphasv[n] * (ctrm3 - rtrm2[n])
                end
                #csum = csum * norm(lvec(next(i))) # Needed for orig. defn. of RWG 
                #                                          # basis funct.
                #  Add to sum total fourier transform with proper sign:
                if rwgdat.bff[1,ib] == iface  # Plus triangle.
                    ft[ib] += csum * (floquet_factor[rwgdat.eci[ie]] / denom)
                else
                    ft[ib] -= csum * (floquet_factor[rwgdat.eci[ie]] / denom)
                end
            end
        end
    end
    
    return ft
end # function
        
        
"""
    j₀(x)
    
Spherical Bessel function of the first kind of order 0 and argument x.
"""
j₀(x::Real) =     abs(x) < 1e-3 ? begin x² = x*x; 1 - x²/6*(1-x²/20) end : sin(x)/x
j₀(x::Complex) = abs2(x) < 1e-6 ? begin x² = x*x; 1 - x²/6*(1-x²/20) end : sin(x)/x

"""
    j₁(x)
    
Spherical Bessel function of the first kind of order 1 and argument x.
"""
j₁(x::Real) =     abs(x) < 1e-3 ? x * (1/3 - x*x/30) : begin s,c = sincos(x); (s/x - c) / x end
j₁(x::Complex) = abs2(x) < 1e-6 ? x * (1/3 - x*x/30) : begin s,c = sincos(x); (s/x - c) / x end
    
"""
    vtxcrd_m(iface::Int, sheet::Sheet, one_meter::Real)

Return the coordinates (in meters) of the triangle vertices for face iface.
"""
@inline function vtxcrd_m(iface, sheet, one_meter)
    vi = @view sheet.fv[:,iface] # Vertex indices
    sheet.ρ[vi] ./ one_meter
end


zdotaxb(a,b) = a[1] * b[2] - a[2] * b[1]

            
end # module
