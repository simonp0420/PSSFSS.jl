module Sheets

export RWGSheet, read_sheet_data, write_sheet_data, find_unique_periods
export rotate!, combine, recttri, SV2, MV2

using StaticArrays: SVector, MVector, SMatrix
using ..PSSFSSLen
using JLD2
using LinearAlgebra: norm
using RecipesBase

const MV2 = MVector{2,Float64}
const SV2 = SVector{2,Float64}

abstract type Sheet end


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
  fr::Vector{Float64} # Face resistance list. fr[i] is the sheet resistance of face i (Ω/□)
  # The following fields are storage for face/face integrals:
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
  
  class::Char # Type of sheet. 'J' for electric current, 'M' for magnetic current
  info::String  # Informational comment
  # The following flag tells rwg_setup whether (.true.) or not (.false.)
  # to check for consistent edges at xi or eta = 0 and 1.  The default
  # value (.TRUE.) means that the check should be performed.
  ξη_check::Bool 
  # The following flag tells rwg_setup whether (.true.) or not (.false.)
  # to Find Unique Face Pairs.
  fufp::Bool
end # struct
import Base.==
==(sh1::RWGSheet, sh2::RWGSheet) = all((getfield(sh1,f)==getfield(sh2,f) for 
                                                    f in fieldnames(RWGSheet)))

# Add a zero-argument constructor:
RWGSheet() = RWGSheet("", u"mm",            # style, units
                      SV2([0.0,0.0]),       # s₁
                      SV2([0.0,0.0]),       # s₂
                      SV2([0.0,0.0]),       # β₁
                      SV2([0.0,0.0]),       # β₂
                      0.0, 0.0, 0.0,        # dx, dy, rot
                      SV2[],                # ρ
                      Int[], Int[],         # e1, e2
                      Array{Int}(undef,0,0), # fv
                      Array{Int}(undef,0,0), # fe
                      Float64[],            # fr
                      ComplexF64[],         # J
                      ComplexF64[],         # J_ξ
                      ComplexF64[],         # J_η
                      ComplexF64[],         # K
                      ComplexF64[],         # K_ξ
                      ComplexF64[],         # K_η
                      Array{SV2}(undef,0),  # ρ_r
                      Float64[],            # rinv
                      0.0, 0.0, 0.0,        # ψ₁, ψ₂, u
                      ' ', "",              # class, info
                      true, false)          # ξη_check, fufp
                    

Base.show(io::IO, ::MIME"text/plain", s::RWGSheet) =
    print(io, "RWGSheet: style=", s.style, ", class=", s.class, ", ", length(s.ρ), " nodes, ", length(s.e1), 
            " edges, ", size(s.fv,2), " faces")
             
                    
"""
    read_sheet_data(filename::AbstractString)::RWGSheet
    
Read the sheet geometry data from a `JLD2` file named in `filename`.
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

Write the sheet geometry data to a `JLD2` file named in `filename`.
"""   
function write_sheet_data(filename::AbstractString, sheet::RWGSheet)
    jldopen(filename, "w") do file
        file["sheet"] = sheet
    end
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
    s1s2 = vcat(map(x -> hcat(x.s₁..., x.s₂...), sheets)...) # Each row is s1x s1y s2x s2y
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

Rotate a sheet by rot degrees (counter-clockwise).
"""
function rotate!(sh::RWGSheet, rot::Real)
    rot == 0 && return
    s,c = sincosd(rot)
    rotmat = SMatrix{2,2}([c -s;s c])
    sh.s₁ = rotmat * sh.s₁
    sh.s₂ = rotmat * sh.s₂
    sh.β₁ = rotmat * sh.β₁
    sh.β₂ = rotmat * sh.β₂
    for n in eachindex(sh.ρ)
        sh.ρ[n] = rotmat * sh.ρ[n]
    end
    sh.rot = rot
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
    Nvcen = 0
    Necen = 0
    vcen1 = Int[]
    vcen2 = Int[]
    if dup_coor ≠ ' '
        tol = 0.5e-4 * norm(sh1.ρ[sh1.e1[1]] - sh1.ρ[sh1.e2[1]])
        for i in 1:length(sh2.ρ)
            if dup_coor == 'x'
                test = sh2.ρ[i][1]
            elseif dup_coor == 'y'
                test = sh2.ρ[i][2]
            else
                test = typemax(typeof(test))
            end
            if abs(test - dup_coor_value) < tol
                # Test to see if there is a sh1 vertex at same coordinate:
                n1match = 0
                for n1 in 1:length(sh1.ρ)
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
        Necen = Nvcen - 1 # Number of shared edges.
        ecen1 = zeros(Int, Necen); ecen2 = zeros(Int, Necen)
        # Locate and save indices in sh1 and sh2 of edges along the center line:
        Necen = 0
        for i in 1:length(sh2.e1) # Loop over sh2 edges
            if (sh2.e1[i] in vcen2) && (sh2.e2[i] in vcen2)
                # Edge i of sh2 is a shared edge.
                Necen += 1
                ecen2[Necen] = i
                i1 = findfirst(x->x==sh2.e1[i], vcen2)
                i2 = findfirst(x->x==sh2.e2[i], vcen2)
                # Now find matching sh1 edge:
                for j1 in 1:length(sh1.e1)
                    if  ((sh1.e1[j1] == vcen1[i1]) && (sh1.e2[j1] == vcen1[i2])) ||
                        ((sh1.e2[j1] == vcen1[i1]) && (sh1.e1[j1] == vcen1[i2]))
                        ecen1[Necen] = j1
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
    sh3.e1 = zeros(Int, length(sh1.e1) + length(sh2.e1) - Necen)
    sh3.e2 = zeros(Int, length(sh1.e2) + length(sh2.e2) - Necen)
    sh3.ρ = Vector{SV2}(undef, length(sh1.ρ) + length(sh2.ρ) - Nvcen)
    sh3.fe = zeros(Int, 3, size(sh1.fe,2) + size(sh2.fe,2))
    sh3.fv = zeros(Int, 3, size(sh1.fv,2) + size(sh2.fv,2))

    # Copy vertex locations:
    sh3.ρ[1:length(sh1.ρ)] = sh1.ρ
    if dup_coor == ' ' 
        sh3.ρ[(1+length(sh1.ρ)):end] = sh2.ρ
    else
        i2 = 1 + length(sh1.ρ)  # Node counter 
        for i in 1:length(sh2.ρ)
            if !(i in vcen2[1:Nvcen])
                sh3.ρ[i2] = sh2.ρ[i]
                i2 += 1
            end
        end
    end
    # Copy edge-node matrices
    eoffset = length(sh1.e1)
    voffset = length(sh1.ρ)
    sh3.e1[1:length(sh1.e1)] = sh1.e1
    sh3.e2[1:length(sh1.e2)] = sh1.e2
    if dup_coor == ' ' 
        sh3.e1[(eoffset+1):end] = sh2.e1 .+ voffset
        sh3.e2[(eoffset+1):end] = sh2.e2 .+ voffset
    else
        for i in 1:length(sh2.e1)
            if i in ecen2
                # Edge i of sh2 is located on duplication line.
                # Do not include this duplicate edge in sh3, but
                # decrement the edge index offset:
                eoffset -= 1
            else
                if sh2.e1[i] in vcen2
                    # The initial point of edge i of sh2 is a duplicate vertex.
                    i2 = findfirst(x->x==sh2.e1[i], vcen2)
                    sh3.e1[i+eoffset] = vcen1[i2]
                else
                    # Ordinary point
                    sh3.e1[i + eoffset] = sh2.e1[i] + voffset
                end
                #
                if sh2.e2[i] in vcen2
                    # The terminal point of edge i of sh2 is on the duplication edge.
                    i2 = findfirst(x->x==sh2.e2[i], vcen2)
                    sh3.e2[i+eoffset] = vcen1[i2]
                else
                    # Ordinary point:
                    sh3.e2[i+eoffset] = sh2.e2[i] + voffset
                end
            end
        end
        # Correct vertex indices:
        for i in Nvcen:-1:1
            sh3.e1[sh3.e1 .> (length(sh1.ρ) + vcen2[i])] .-= 1
            sh3.e2[sh3.e2 .> (length(sh1.ρ) + vcen2[i])] .-= 1
        end
    end
    # Copy face/vertex matrix
    sh3.fv[:, 1:size(sh1.fv,2)] = sh1.fv
    sh3.fv[:, 1+size(sh1.fv,2):end] = sh2.fv .+ voffset # offset will be corrected later
    #  Correct duplicate vertices from sh2:
    for n2 in 1:Nvcen # Examine each duplicate vertex in sh2
        n1 = vcen1[n2] # Initialize index of matching vertex in sh1
        # Replace all references to vertex vcen2[n2] with ref to vcen1[n1]
        sh3.fv[sh3.fv .== vcen2[n2] + voffset] .= n1
    end
    # Correct vertex indices:
    for i in Nvcen:-1:1
        sh3.fv[sh3.fv .> voffset + vcen2[i]] .-= 1
    end
    # Copy face/edge matrix:
    foffset = size(sh1.fv,2)
    eoffset = length(sh1.e1)
    sh3.fe[:,1:size(sh1.fe,2)] = sh1.fe
    sh3.fe[:,1+size(sh1.fe,2):end] = sh2.fe .+ eoffset # offset will be corrected later
    #  Correct duplicate edges from sh2:
    if dup_coor ≠ ' '
        for n2 in 1:length(ecen2) # Examine each duplicate edge in sh2
            e2 = ecen2[n2]  # Index of duplicate edge in sh2.
            e1 = ecen1[n2]  # Index of duplicate edge in sh1.
            sh3.fe[sh3.fe .== (e2 + eoffset)] .= e1
        end
        # Correct edge indices:
        for i in length(ecen2):-1:1
            sh3.fe[sh3.fe .> eoffset + ecen2[i]] .-= 1
        end
    end
    return sh3
end


"""
    recttri(rhobl::SVector{2,Float64}, rhotr::SVector{2,Float64}, nx::Int, ny::Int)

Create a variable of type `RWGSheet` that contains the triangulation for 
a rectangular strip.  The fields `ρ`, `e1`, `e2`, `fv`, and `fe` properly initialized.
"""
function recttri(rhobl::SV2, rhotr::SV2, nx::Int, ny::Int)   
    nodecount = (nx+1) * (ny+1)  # Number of nodes.
    edgecount = 3*nx*ny + nx + ny  # Number of edges.
    facecount = 2*nx*ny  # Number of faces.
    
    sh = RWGSheet()
    sh.ρ = Vector{SV2}(undef, nodecount)
    
    # Set the node coordinates:
    drho = (rhotr - rhobl) ./ [nx, ny]
    n = 0  # Initialize node index.
    for j in 0:ny
        yj = j * drho[2]
        for i in 0:nx
            n += 1 
            sh.ρ[n] = rhobl + SV2([i*drho[1], yj])
        end
    end

    sh.e1 = zeros(Int, edgecount)
    sh.e2 = zeros(Int, edgecount)
    e = 0  # Initialize edge index.
    # Do the horizontal edges:
    for j in 0:ny
        kadd = j * (nx+1)
        for i in 1:nx
            e += 1
            sh.e1[e] = i + kadd
            sh.e2[e] = sh.e1[e] + 1
        end
    end
    # Do the vertical edges:
    for j in 1:ny
        kadd = (j-1) * (nx+1) + 1
        for i in 0:nx
            e += 1
            sh.e1[e] = i + kadd
            sh.e2[e] = sh.e1[e] + (nx+1)
        end
    end
    # Do the diagonal edges:
    for j in 1:ny
        kadd1 = (j-1) * (nx+1)
        kadd2 = 1 + j * (nx+1)
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
    nhe = nx*ny + nx  # Number of horizontal edges.
    nve = nx*ny + ny  # Number of vertical edges.
    nde = nx*ny       # Number of diagonal edges
    f = 0  # Initialize face index.
    for j in 1:ny
        nadd1 = (j-1) * (nx+1)
        nadd2 = 1 + j * (nx+1)
        for i in 1:nx
            f += 1  # Bump face index (upper left face).
            sh.fv[1,f] = i + nadd1  # Lower Left vertex.
            sh.fv[2,f] = i + nadd2  # Upper right vertex.
            sh.fv[3,f] = i + nadd2 - 1 # Upper left vertex.
            sh.fe[1,f] = i + j*nx # Upper edge.
            sh.fe[2,f] = i + nhe + nadd1 # Left edge
            sh.fe[3,f] = i + nhe + nve + (j-1)*nx # Diagonal edge
            f += 1  # Bump face index (lower right face).
            sh.fv[1,f] = sh.fv[1,f-1] # Lower Left vertex.
            sh.fv[2,f] = 1 + sh.fv[1,f-1] # Lower right vertex.
            sh.fv[3,f] = sh.fv[2,f-1] # Upper right vertex.
            sh.fe[1,f] = 1 + sh.fe[2,f-1] # Right edge.
            sh.fe[2,f] = sh.fe[3,f-1] # Diagonal edge.
            sh.fe[3,f] = sh.fe[1,f-1] - nx # Bottom edge
        end
    end

    return sh
end

@recipe function f(sh::RWGSheet; edges=true, faces=false, nodes=false, 
                    edgenumbers=false, facenumbers=false, nodenumbers=false,
                    unitcell=false, rep=(1,1), fontsize=9)
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
        x0, y0 = (m-1)*sh.s₁ + (n-1)*sh.s₂

        # Add series for faces
        if faces
            for i in 1:size(sh.fv,2)
                points = sh.ρ[sh.fv[:,i]]
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
                push!(x,NaN)
                push!(y,NaN)
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
            points = [0*sh.s₁, sh.s₁, sh.s₁+sh.s₂, sh.s₂]
            x = x0 .+ [point[1] for point in points]; push!(x, x[1])
            y = y0 .+ [point[2] for point in points]; push!(y, y[1])
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
                markersize -->  1
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
                x,y
            end
        end

        # Add series for edge numbers
        if edgenumbers
            x = zeros(Float64, length(sh.e1))
            y = zeros(Float64, length(sh.e1))
            for i in 1:length(sh.e1)
                x[i], y[i] = 0.5*sum(sh.ρ[[sh.e1[i], sh.e2[i]]])
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
                x,y
            end
        end

        # Add series for face numbers
        if facenumbers
            x = zeros(Float64, size(sh.fv,2))
            y = zeros(Float64, size(sh.fv,2))
            for i in 1:size(sh.fv,2)
                x[i], y[i] = (1/3)*sum(sh.ρ[sh.fv[:,i]])
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
                x,y
            end
        end
    end
end

end # module
