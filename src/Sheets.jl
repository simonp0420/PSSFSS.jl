module Sheets

export RWGSheet, read_sheet_data, write_sheet_data, find_unique_periods
export rotate!, translate!, combine, recttri, SV2, MV2, nodecount, facecount, edgecount
export export_sheet, STL_ASCII, STL_BINARY

using StaticArrays: SVector, MVector, SMatrix
using ..PSSFSSLen
using JLD2
using LinearAlgebra: norm
using RecipesBase
using Printf: @printf


const MV2 = MVector{2,Float64}
const SV2 = SVector{2,Float64}

abstract type Sheet end

"""
    RWGSheet
A type that represents a zero-thickness sheet of periodically patterned metalization.  Particular instances
are created by calling a constructor function for a specific type of sheet geometry.  These include:
`diagstrip`, `jerusalemcross`, `loadedcross`, `manji`, `meander`, `pecsheet`, `pmcsheet`, `polyring`,
`rectstrip`, `sinuous`, and `splitring`.
"""
mutable struct RWGSheet <: Sheet
    style::String
    units::PSSFSSLength # Length unit
    s₁::SV2 # Direct lattice vector (specified units)
    s₂::SV2 # Direct lattice vector (specified units)
    β₁::SV2 # Reciprocal lattice vector (1/(specified units))
    β₂::SV2 # Reciprocal lattice vector (1/(specified units))
    dx::Float64 # Unit cell displacment in x (in specified units)
    dy::Float64 # Unit cell displacment in y (in specified units)
    rot::Float64 # Rotation angle for unit cell (deg)
    ρ::Vector{SV2} # Node coordinates
    e1::Vector{Int} # Edge connect. list. e1[i] is the initial node of edge i
    e2::Vector{Int} # Edge connect. list. e2[i] is the terminal node of edge i
    fv::Array{Int,2} #  Face/vertex list. fv[:,i] lists vertices of face i
    fe::Array{Int,2} #  Face/edge list. fe[:,i] lists edges of face i
    # The following fields are storage for face/face integrals:
    I1::Vector{ComplexF64}
    I1_ξ::Vector{ComplexF64}
    I1_η::Vector{ComplexF64}
    I2::Vector{ComplexF64}
    J::Vector{ComplexF64}
    J_ξ::Vector{ComplexF64}
    J_η::Vector{ComplexF64}
    K::Vector{ComplexF64}
    K_ξ::Vector{ComplexF64}
    K_η::Vector{ComplexF64}
    ρ_r::Vector{SV2}
    rinv::Vector{Float64}
    # Parameters that the face/face integrals depend on:
    ψ₁::Float64  # Incremental phase shift (radians)
    ψ₂::Float64  # Incremental phase shift (radians)
    u::Float64   # Smoothing parameter (1/(specified units))

    class::Char # Sheet classifier. 'J' for electric current, 'M' for magnetic current, 'E' for PEC/E-wall, 'H' for PMC/H-wall
    info::String  # Informational comment
    # The following flag tells rwg_setup whether (.true.) or not (.false.)
    # to check for consistent edges at xi or eta = 0 and 1.  The default
    # value (.TRUE.) means that the check should be performed.
    ξη_check::Bool
    # The following flag tells rwg_setup whether (.true.) or not (.false.)
    # to Find Unique Face Pairs.
    fufp::Bool
    σ::Float64 # Bulk, DC conductivity for J-class surfaces [S/m].  If < 0, ignore it.
    Rq::Float64 # RMS surface roughness for J-class [m]
    disttype::Symbol # :normal or :rayleigh
    Zs::ComplexF64 # Surface impedance [Ω] for J-class. If σ₀ > 0, recompute for each frequency.
end # struct
import Base.==
==(sh1::RWGSheet, sh2::RWGSheet) = all((getfield(sh1, f) == getfield(sh2, f) for
                                        f in fieldnames(RWGSheet)))

# Add a zero-argument constructor:
RWGSheet() = RWGSheet("", u"mm",            # style, units
    SV2([0.0, 0.0]),       # s₁
    SV2([0.0, 0.0]),       # s₂
    SV2([0.0, 0.0]),       # β₁
    SV2([0.0, 0.0]),       # β₂
    0.0, 0.0, 0.0,        # dx, dy, rot
    SV2[],                # ρ
    Int[], Int[],         # e1, e2
    Array{Int}(undef, 0, 0), # fv
    Array{Int}(undef, 0, 0), # fe
    ComplexF64[],         # I1
    ComplexF64[],         # I1_ξ
    ComplexF64[],         # I1_η
    ComplexF64[],         # I2
    ComplexF64[],         # J
    ComplexF64[],         # J_ξ
    ComplexF64[],         # J_η
    ComplexF64[],         # K
    ComplexF64[],         # K_ξ
    ComplexF64[],         # K_η
    Array{SV2}(undef, 0),  # ρ_r
    Float64[],            # rinv
    0.0, 0.0, 0.0,        # ψ₁, ψ₂, u
    ' ', "",              # class, info
    true, false,          # ξη_check, fufp
    -Inf, 0.0, :normal, 0.0im)     # σ, Rq, disttype, Zs

"""
    nodecount(s::RWGSheet)

Return the number of unique triangle vertices in the `RWGSheet` triangulation.
"""
nodecount(s::RWGSheet) = length(s.ρ)


"""
    facecount(s::RWGSheet)

Return the number of triangle faces in the `RWGSheet` triangulation.
"""
facecount(s::RWGSheet) = size(s.fv, 2)


"""
    edgecount(s::RWGSheet)

Return the number of triangle edges in the `RWGSheet` triangulation.
"""
edgecount(s::RWGSheet) = length(s.e1)

function Base.show(io::IO, ::MIME"text/plain", s::RWGSheet)
    if s.class == 'E'
        print(io, "RWGSheet: perfect electric conducting wall")
    elseif s.class == 'H'
        print(io, "RWGSheet: perfect magnetic conducting wall")
    else
        if s.σ > 0
            print(io, "RWGSheet: style=", s.style, ", class=", s.class, ", ", nodecount(s), " nodes, ", edgecount(s),
            " edges, ", facecount(s), " faces, σ=", s.σ, " S/m, Rq=", s.Rq, " m (:", s.disttype, ")")
        else
            print(io, "RWGSheet: style=", s.style, ", class=", s.class, ", ", nodecount(s), " nodes, ", edgecount(s),
                " edges, ", facecount(s), " faces, Zs=", s.Zs, " Ω")
        end
    end
end

"""
    read_sheet_data(filename::AbstractString)::RWGSheet

Read the sheet triangulation and unit cell data from a `JLD2` file named in `filename`.
"""
function read_sheet_data(filename::AbstractString)::RWGSheet
    jldopen(filename, "r") do file
        try
            return file["sheet"]
        catch
            @error "$(filename) does not contain sheet data"
        end
    end
end # function

"""
    write_sheet_data(filename::AbstractString, sheet::RWGSheet)

Write the sheet triangulation and unit cell data to a `JLD2` file named in `filename`.
"""
function write_sheet_data(filename::AbstractString, sheet::RWGSheet)
    jldopen(filename, "w") do file
        file["sheet"] = sheet
    end
end

abstract type CAD_Export end
struct STL_ASCII <: CAD_Export end
struct STL_BINARY <: CAD_Export end

# The code for export_sheet was adapted from that in package MeshIO
export_sheet(fname::AbstractString, sheet::RWGSheet, t) = error("Unknown CAD export type $t")

"""
    export_sheet(fname::AbstractString, sheet::RWGSheet, export_type)

Export an `RWGSheet` triangulation to an STL CAD file.  `export_type` may be either
`STL_ASCII` or `STL_BINARY`.
"""
function export_sheet(fname::AbstractString, sheet::RWGSheet, export_type::Type{STL_ASCII})
    open(fname, "w") do io
        write(io, "solid vcg\n") # header
        for i in eachcol(sheet.fv)
            @printf(io, "  facet normal %e %e %e\n", 0.0, 0.0, 1.0)
            write(io, "    outer loop\n")
            for v in @view sheet.ρ[i]
                @printf(io, "      vertex  %e %e %e\n", v[1], v[2], 0.0)
            end
            write(io, "    endloop\n")
            write(io, "  endfacet\n")
        end
        write(io,"endsolid vcg\n")
    end
    return nothing
end

function export_sheet(fname::AbstractString, sheet::RWGSheet, export_type::Type{STL_BINARY})
    open(fname, "w") do io
          # Implementation made according to https://en.wikipedia.org/wiki/STL_%28file_format%29#Binary_STL
        for i in 1:80 # write empty header
            write(io, 0x00)
        end

        write(io, UInt32(facecount(sheet))) # write triangle count
        n = (0.0f0, 0.0f0, 1.0f0)
        for i in eachcol(sheet.fv)
            foreach(j -> write(io, n[j]), 1:3)
            for point in @view sheet.ρ[i]
                point3 = (Float32(point[1]), Float32(point[2]), 0.0f0)
                foreach(p -> write(io, p), point3)
            end
            write(io, 0x0000) # write 16bit empty bit
        end
    end
    return nothing
end



"""
    find_unique_periods(junction::Vector{Int}, sheets)

Find the unique unit cells for the sheets used in the FSS analysis.

# Arguments
- `junction`:  An integer array of length `(Nlayer-1)` containing
               in location `i` the index of the FSS sheet located
               at the interface of dielectric layers `i` and `i+1`. If
               no sheet is present there, the value is 0.
- `sheets`:     An iterable that contains the FSS sheets.
# Return Value
- `upa`            (Unique Periodicity Array) An integer array of the
                   same length of junction, containing zeros in the same
                   locations.  The nonzero entries correspond to sheet
                   locations, and are numbered consecutively according
                   to the equivalence class of the sheet at that location.
                   Two sheets are equivalent if they have the same unit cell.
"""
function find_unique_periods(junction::Vector{Int}, sheets)
    all(t isa Sheet for t in sheets) || error("Elements of sheets must be of type Sheet")
    one_meter = map(x -> ustrip(Float64, x.units, 1.0u"m"), sheets)
    s1s2 = mapreduce(x -> [x.s₁[1] x.s₁[2] x.s₂[1] x.s₂[2]], vcat, sheets) # Each row is s1x s1y s2x s2y
    s1s2 = s1s2 ./ one_meter # All rows now are comparable (in meters)
    s1s2 = round.(s1s2, sigdigits=8)

    upa = zeros(Int, length(junction))
    Nup = 0  # Initialize Number of Unique Periodicities.
    for i in 1:length(upa)  # Step through each junction.
        isht = junction[i] # Sheet index.
        ((isht == 0) || (sheets[isht].style == "NULL")) && continue
        # Compare s1 and s2 of current (isht) sheet with previous sheets.
        for n in 1:Nup  # compare to one member of each equivalence class.
            # Find a sheet that is in equivalence class n:
            nsht = junction[findfirst(isequal(n), upa)]  # Index of sheet to be compared.
            # Compare unit cell of sheet nsht with that of sheet isht:
            if view(s1s2, isht, :) == view(s1s2, nsht, :)
                # Sheets are in the same equivalence class
                upa[i] = n  # Store equiv. class number.
                @goto NextOuterFor
            end
        end
        # If execution fell through to here, we found a sheet that is not
        # in an existing equivalence class.
        Nup += 1 # Bump count of equiv. classes.
        upa[i] = Nup  # Store equiv. class number.
        @label NextOuterFor
    end
    return upa
end


"""
    rotate!(sh::RWGSheet, rot::Real)

Rotate a sheet by rot degrees (counter-clockwise).  The entire unit cell and its
contents are rotated.
"""
function rotate!(sh::RWGSheet, rot::Real)
    rot == 0 && return
    s, c = sincosd(rot)
    rotmat = SMatrix{2,2}([c -s; s c])
    sh.s₁ = rotmat * sh.s₁
    sh.s₂ = rotmat * sh.s₂
    sh.β₁ = rotmat * sh.β₁
    sh.β₂ = rotmat * sh.β₂
    for n in eachindex(sh.ρ)
        sh.ρ[n] = rotmat * sh.ρ[n]
    end
    sh.rot = rot
    return sh
end

"""
    orient!(sh::RWGSheet, rot::Real, center::AbstractVector)

Rotate a sheet by rot degrees (counter-clockwise) about the point specified by `center`.
Only the content of the unit cell (the triangulation) is rotated.  The unit cell is left unchanged.
"""
function orient!(sh::RWGSheet, rot::Real, center::AbstractVector)
    rot == 0 && return
    s, c = sincosd(rot)
    rotmat = SMatrix{2,2}([c -s; s c])
    for n in eachindex(sh.ρ)
        sh.ρ[n] = center + rotmat * (sh.ρ[n] - center)
    end
    return sh
end


"""
    translate!(sh::RWGSheet, dx, dy)

Translate a sheet by dx in x and dy in y.
"""
function translate!(sh::RWGSheet, dx::Real, dy::Real)
    dx == dy == 0 && (return sh)
    tvec = [dx, dy]
    for n in eachindex(sh.ρ)
        sh.ρ[n] = tvec + sh.ρ[n]
    end
    return sh
end


"""
    combine(sh1::RWGSheet, sh2::RWGSheet, dup_coor::Char, dup_coor_value::Real)

Combine the triangulations stored in sheets `sh1` and `sh2`.

# Arguments:
- `sh1`, `sh2`:  sheets having
              initialized values for fields `units`, `ρ`,
              `e1`, `e2` `fv`, `fe`, and possibly `fr`.  It is assumed
              that the two triangulations do not overlap except possibly
              along a line defined by `dup_coor` and `dup_coor_value`,
              as discussed below.  If they do coincide along such
              a line, then they must share the same set of vertices
              and edges along this line.  These duplicate vertices
              and edges will be removed by this routine.
- `dup_coor`    Either 'x' or 'y' to indicate at which coordinate
              constant line the two triangulations may overlap,
              requiring redundant edges and nodes to be removed,
              or ' ' indicating that no search for duplicate nodes
              is required.
- `dup_coor_value`  The value of the coordinate at which the two
              input triangulations overlap.
# Return value
- `sh3`         A `RWGSheet` instance with the following member arrays
                initialized: `units`, `ρ`, `ec`, `fv`, `fe`, and `fr`.
"""
function combine(sh1::RWGSheet, sh2::RWGSheet, dup_coor::Char, dup_coor_value::Real)
    sh1.units == sh2.units || error("Inconsistent units for sh1 and sh2")
    # Count number of vertices located at the duplicate coordinate.
    # Save vertex indices of matching points in vcen1 and vcen2
    #
    Nv1 = length(sh1.ρ)
    Ne1 = length(sh1.e1)
    Nf1 = size(sh1.fe, 2)
    Nv2 = length(sh2.ρ)
    Ne2 = length(sh2.e1)
    Nf2 = size(sh2.fe, 2)

    # Find coincident vertices
    Nvcen = 0
    vcen1 = Int[] # coincident vertices in sh1
    vcen2 = Int[] # coincident vertices in sh2
    if dup_coor == ' '
    elseif dup_coor == 'x'
        ρind = 1
    elseif dup_coor == 'y'
        ρind = 2
    else
        throw(ArgumentError("Illegal value for dup_coor: '$dup_coor'"))
    end
    if dup_coor ≠ ' '
        tol = 0.5e-4 * norm(sh1.ρ[sh1.e1[1]] - sh1.ρ[sh1.e2[1]])
        for i in 1:Nv2
            test = sh2.ρ[i][ρind]
            if abs(test - dup_coor_value) < tol
                # Test to see if there is a sh1 vertex at same coordinate:
                n1match = 0
                for n1 in 1:Nv1
                    if norm(sh1.ρ[n1] - sh2.ρ[i]) < tol
                        n1match = n1
                        break
                    end
                end
                if n1match ≠ 0
                    Nvcen += 1
                    push!(vcen2, i)
                    push!(vcen1, n1match)
                end
            end
        end

        # Locate and save indices in sh1 and sh2 of edges along the center line:
        Necen = 0 # Initialize count of shared center edges
        ecen1 = Int[] # Coincident edges in sh1
        ecen2 = Int[] # Coincident edges in sh1
        for i in 1:Ne2 # Loop over sh2 edges
            if (sh2.e1[i] in vcen2) && (sh2.e2[i] in vcen2)
                # Edge i of sh2 is a shared edge.
                Necen += 1
                push!(ecen2, i)
                i1 = findfirst(==(sh2.e1[i]), vcen2)
                i2 = findfirst(==(sh2.e2[i]), vcen2)
                # vcen1[k] and vcen2[k] refer to nodes in sh1 and sh2, resp, but at the same coordinates
                vcen1nodes = extrema((vcen1[i1], vcen1[i2])) # Same node locations but wrt sh1
                # Now find matching sh1 edge:
                for j1 in 1:Ne1
                    if extrema((sh1.e1[j1], sh1.e2[j1])) == vcen1nodes
                        push!(ecen1, j1)
                        # ecen1[k] and ecen2[k]  refer to edges in sh1 and sh2, resp, but at the same locations
                        @goto sh2edges
                    end
                end
            else
                continue
            end
            error("Unable to find matching duplicate edge in combine_sheet")
            @label sh2edges
        end
    end

    # Allocate triangulation arrays in new sheet:
    sh3 = RWGSheet()
    sh3.e1 = zeros(Int, Ne1 + Ne2) # Will remove Necen of these later
    sh3.e2 = zeros(Int, Ne1 + Ne2) # Will remove Necen of these later
    sh3.ρ = Vector{SV2}(undef, Nv1 + Nv2 - Nvcen)
    sh3.fe = zeros(Int, 3, Nf1 + Nf2)
    sh3.fv = zeros(Int, 3, Nf1 + Nf2)

    # Copy vertex locations:
    sh3.ρ[1:Nv1] = sh1.ρ
    if dup_coor == ' '
        sh3.ρ[(1+Nv1):end] = sh2.ρ
    else
        i2 = 1 + Nv1  # Node counter
        for i in 1:Nv2
            if !(i in vcen2)
                sh3.ρ[i2] = sh2.ρ[i]
                i2 += 1
            end
        end
    end

    # Copy edge-node matrices
    eoffset = Ne1
    voffset = Nv1
    sh3.e1[1:Ne1] = sh1.e1
    sh3.e2[1:Ne1] = sh1.e2
    sh3.e1[Ne1+1:end] .= sh2.e1 .+ voffset # Offset will be corrected later if necessary
    sh3.e2[Ne1+1:end] .= sh2.e2 .+ voffset # Offset will be corrected later if necessary

    # Correct if there are duplicate edges
    sh3e1new = @view sh3.e1[Ne1+1:end]
    sh3e2new = @view sh3.e2[Ne1+1:end]
    if dup_coor ≠ ' '
        for i in eachindex(sh2.e1, sh2.e2, sh3e1new, sh3e2new)
            if i ∉ ecen2 # Edges in ecen2 will be removed later
                i2 = findfirst(==(sh2.e1[i]), vcen2)
                if isnothing(i2)
                    # ordinary point not on shared border:
                    sh3e1new[i] -= count(<(sh2.e1[i]), vcen2)
                else
                    # initial point of edge i of sh2 is a duplicate vertex:
                    sh3e1new[i] = vcen1[i2]
                end
                #
                i2 = findfirst(==(sh2.e2[i]), vcen2)
                if isnothing(i2)
                    # ordinary point not on shared border:
                    sh3e2new[i] -= count(<(sh2.e2[i]), vcen2)
                else
                    # initial point of edge i of sh2 is a duplicate vertex:
                    sh3e2new[i] = vcen1[i2]
                end
            end
        end
        # Remove duplicate edges
        delindices = Ne1 .+ ecen2
        deleteat!(sh3.e1, delindices)
        deleteat!(sh3.e2, delindices)
    end

    # Copy face/vertex matrix
    sh3.fv[:, 1:Nf1] = sh1.fv
    sh3.fv[:, 1+Nf1:end] = sh2.fv .+ voffset # offset will be corrected later
    if dup_coor ≠ ' '
        #  Correct duplicate vertices from sh2:
        sh3fvnew = @view sh3.fv[:, 1+Nf1:end]
        for n in eachindex(sh3fvnew, sh2.fv)
            v2 = sh2.fv[n]
            ncen = findfirst(==(v2), vcen2)
            if isnothing(ncen)
                # ordinary, non-duplicate point
                sh3fvnew[n] -= count(<(v2), vcen2) # correct offset
            else
                sh3fvnew[n] = vcen1[ncen] # Replace duplicate node
            end
        end
    end

    # Copy face/edge matrix:
    eoffset = length(sh1.e1)
    sh3.fe[:, 1:Nf1] = sh1.fe
    sh3.fe[:, 1+Nf1:end] .= sh2.fe .+ eoffset # offset will be corrected later
    #  Correct duplicate edges from sh2:
    if dup_coor ≠ ' '
        sh3fenew = @view sh3.fe[:, 1+Nf1:end]
        for n in eachindex(sh2.fe, sh3fenew)
            e2 = sh2.fe[n]
            ncen = findfirst(==(e2), ecen2)
            if isnothing(ncen)
                # ordinary (non-duplicate) edge
                sh3fenew[n] -= count(<(e2), ecen2) # correct the offset
            else
                sh3fenew[n] = ecen1[ncen]
            end
        end
    end

    test_fefv(sh3)

    return sh3
end


"""
    recttri(rhobl::SVector{2,Float64}, rhotr::SVector{2,Float64}, nx::Int, ny::Int)

Create a variable of type `RWGSheet` that contains the triangulation for
a rectangular strip.  The fields `ρ`, `e1`, `e2`, `fv`, and `fe` properly initialized.
"""
function recttri(rhobl::SV2, rhotr::SV2, nx::Int, ny::Int)
    nodecount = (nx + 1) * (ny + 1)  # Number of nodes.
    edgecount = 3 * nx * ny + nx + ny  # Number of edges.
    facecount = 2 * nx * ny  # Number of faces.

    sh = RWGSheet()
    sh.ρ = Vector{SV2}(undef, nodecount)

    # Set the node coordinates:
    drho = (rhotr - rhobl) ./ [nx, ny]
    n = 0  # Initialize node index.
    for j in 0:ny
        yj = j * drho[2]
        for i in 0:nx
            n += 1
            sh.ρ[n] = rhobl + SV2([i * drho[1], yj])
        end
    end

    sh.e1 = zeros(Int, edgecount)
    sh.e2 = zeros(Int, edgecount)
    e = 0  # Initialize edge index.
    # Do the horizontal edges:
    for j in 0:ny
        kadd = j * (nx + 1)
        for i in 1:nx
            e += 1
            sh.e1[e] = i + kadd
            sh.e2[e] = sh.e1[e] + 1
        end
    end
    # Do the vertical edges:
    for j in 1:ny
        kadd = (j - 1) * (nx + 1) + 1
        for i in 0:nx
            e += 1
            sh.e1[e] = i + kadd
            sh.e2[e] = sh.e1[e] + (nx + 1)
        end
    end
    # Do the diagonal edges:
    for j in 1:ny
        kadd1 = (j - 1) * (nx + 1)
        kadd2 = 1 + j * (nx + 1)
        for i in 1:nx
            e += 1
            sh.e1[e] = i + kadd1
            sh.e2[e] = i + kadd2
        end
    end

    # Done with edges.  Begin setting up faces.
    # Allocate arrays whose length depends only on the number of faces:
    sh.fv = zeros(Int, 3, facecount)
    sh.fe = zeros(Int, 3, facecount)
    # Set up the face/vertex and face/edge matrices:
    nhe = nx * ny + nx  # Number of horizontal edges.
    nve = nx * ny + ny  # Number of vertical edges.
    nde = nx * ny       # Number of diagonal edges
    f = 0  # Initialize face index.
    for j in 1:ny
        nadd1 = (j - 1) * (nx + 1)
        nadd2 = 1 + j * (nx + 1)
        for i in 1:nx
            f += 1  # Bump face index (upper left face).
            sh.fv[1, f] = i + nadd1  # Lower Left vertex.
            sh.fv[2, f] = i + nadd2  # Upper right vertex.
            sh.fv[3, f] = i + nadd2 - 1 # Upper left vertex.
            sh.fe[1, f] = i + j * nx # Upper edge.
            sh.fe[2, f] = i + nhe + nadd1 # Left edge
            sh.fe[3, f] = i + nhe + nve + (j - 1) * nx # Diagonal edge
            f += 1  # Bump face index (lower right face).
            sh.fv[1, f] = sh.fv[1, f-1] # Lower Left vertex.
            sh.fv[2, f] = 1 + sh.fv[1, f-1] # Lower right vertex.
            sh.fv[3, f] = sh.fv[2, f-1] # Upper right vertex.
            sh.fe[1, f] = 1 + sh.fe[2, f-1] # Right edge.
            sh.fe[2, f] = sh.fe[3, f-1] # Diagonal edge.
            sh.fe[3, f] = sh.fe[1, f-1] - nx # Bottom edge
        end
    end

    return sh
end


"Plot recipe for RWGSheet"
@recipe function f(sh::RWGSheet; edges=true, faces=false, nodes=false,
    edgenumbers=false, facenumbers=false, nodenumbers=false,
    unitcell=false, rep=(1, 1), fontsize=9)
    # set a default value for an attribute with `-->`.  Force it with `:=`.
    xguide --> "x ($(sh.units))"
    yguide --> "y ($(sh.units))"
    aspect_ratio := :equal

    if isa(rep[1], Int)
        mrange = 1:rep[1]
    elseif isa(rep[1], UnitRange)
        mrange = rep[1]
    else
        error("Illegal type for rep[1]")
    end
    if isa(rep[2], Int)
        nrange = 1:rep[2]
    elseif isa(rep[2], UnitRange)
        nrange = rep[2]
    else
        error("Illegal type for rep[2]")
    end


    for m in mrange, n in nrange
        x0, y0 = (m - 1) * sh.s₁ + (n - 1) * sh.s₂

        # Add series for faces
        if faces
            for i in 1:size(sh.fv, 2)
                points = sh.ρ[sh.fv[:, i]]
                x = x0 .+ [point[1] for point in points]
                y = y0 .+ [point[2] for point in points]
                @series begin
                    seriestype := :shape
                    # ignore series in legend and color cycling
                    primary := false
                    linecolor := nothing
                    fillcolor --> :blue
                    fillalpha --> 0.8
                    markershape := :none
                    x, y
                end
            end
        end

        # Add series for edges
        if edges
            x = Float64[]
            y = Float64[]
            for i in 1:length(sh.e1)
                points = sh.ρ[[sh.e1[i], sh.e2[i]]]
                push!(x, NaN)
                push!(y, NaN)
                append!(x, [p[1] for p in points])
                append!(y, [p[2] for p in points])
            end
            x .+= x0
            y .+= y0
            @series begin
                seriestype := :path
                # ignore series in legend and color cycling
                primary := false
                linecolor --> :black
                linestyle := :solid
                fillcolor := nothing
                fillalpha := 0
                markershape := :none
                x, y
            end
        end

        # Add series for unit cell
        if unitcell
            points = [0 * sh.s₁, sh.s₁, sh.s₁ + sh.s₂, sh.s₂]
            x = x0 .+ [point[1] for point in points]
            push!(x, x[1])
            y = y0 .+ [point[2] for point in points]
            push!(y, y[1])
            @series begin
                seriestype := :path
                # ignore series in legend and color cycling
                primary := false
                linecolor := :blue
                linestyle := :dot
                fillcolor := nothing
                fillalpha := 0
                markershape := :none
                x, y
            end
        end

        # Add series for nodes
        if nodes
            x = x0 .+ [p[1] for p in sh.ρ]
            y = y0 .+ [p[2] for p in sh.ρ]
            @series begin
                seriestype := :scatter
                # ignore series in legend and color cycling
                primary := false
                linecolor := nothing
                markercolor --> :black
                markershape --> :circle
                markersize --> 1
                x, y
            end
        end

        # Add series for node numbers
        if nodenumbers
            x = x0 .+ [p[1] for p in sh.ρ]
            y = y0 .+ [p[2] for p in sh.ρ]
            @series begin
                seriestype := :scatter
                # ignore series in legend and color cycling
                primary := false
                linecolor := nothing
                markersize := 0
                markeralpha := 0
                markercolor := nothing
                markershape := :none
                annotations := [(x[i], y[i], string(i), fontsize) for i in 1:length(x)]
                x, y
            end
        end

        # Add series for edge numbers
        if edgenumbers
            x = zeros(Float64, length(sh.e1))
            y = zeros(Float64, length(sh.e1))
            for i in 1:length(sh.e1)
                x[i], y[i] = 0.5 * sum(sh.ρ[[sh.e1[i], sh.e2[i]]])
            end
            x .+= x0
            y .+= y0
            @series begin
                seriestype := :scatter
                # ignore series in legend and color cycling
                primary := false
                linecolor := nothing
                markersize := 0
                markeralpha := 0
                markercolor := nothing
                markershape := :none
                annotations := [(x[i], y[i], string(i), fontsize) for i in 1:length(x)]
                x, y
            end
        end

        # Add series for face numbers
        if facenumbers
            x = zeros(Float64, size(sh.fv, 2))
            y = zeros(Float64, size(sh.fv, 2))
            for i in 1:size(sh.fv, 2)
                x[i], y[i] = (1 / 3) * sum(sh.ρ[sh.fv[:, i]])
            end
            x .+= x0
            y .+= y0
            @series begin
                seriestype := :scatter
                # ignore series in legend and color cycling
                primary := false
                linecolor := nothing
                markersize := 0
                markeralpha := 0
                markercolor := nothing
                markershape := :none
                annotations := [(x[i], y[i], string(i), fontsize) for i in 1:length(x)]
                x, y
            end
        end
    end
end

function test_fefv(sheet::RWGSheet)
    ip1 = (2, 3, 1)
    ip2 = (3, 1, 2)
    for iface in axes(sheet.fv, 2)
        for ivlocal in axes(sheet.fv, 1) # local node index within current face
            ivglobal = sheet.fv[ivlocal, iface] # Global node index
            oppedge= sheet.fe[ivlocal, iface] # Should be global index of edge opposite this vertex
            edgen1n2 = Set((sheet.e1[oppedge], sheet.e2[oppedge])) # global node indices of edge vertices
            facen1n2 = Set((sheet.fv[ip1[ivlocal], iface], sheet.fv[ip2[ivlocal], iface])) # global node indices of other face vertices
            if edgen1n2 ≠ facen1n2
                @show edgen1n2
                @show facen1n2
                println("Bad triangle fe/fv matching for face #", iface, " local vertex ", ivlocal)
                println("face edges: " )
                for e in @view sheet.fe[:,iface]
                    println("  edge ", e, " from node ", sheet.e1[e], " to node ", sheet.e2[e])
                end
                print("face nodes: ")
                foreach(x -> print(x," "), @view sheet.fv[:,iface])
                println()
                error("bad")
            end
        end
    end # iface loop
end


end # module
