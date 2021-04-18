module Elements

export rectstrip, polyring, meander, loadedcross, jerusalemcross, pecsheet, pmcsheet

using ..PSSFSSLen: mm, cm, inch, mil, PSSFSSLength
using ..Sheets: RWGSheet, rotate!, combine, recttri, SV2
using ..Meshsub: meshsub
using StaticArrays: SA
using LinearAlgebra: ×, norm, ⋅

macro testpos(var)
  return:($(esc(var)) > 0 || error($(esc(string(var))) * " must be positive!"))
end


zhatcross(t) = [-t[2], t[1]]
zhatcross(t::SV2) = SV2(-t[2], t[1])


"""
    s₁s₂2β₁β₂(s₁,s₂) -> (β₁, β₂)

Compute the reciprocal lattice vectors from the direct lattice vectors.
Inputs and outputs are mutable 2-vectors from StaticArrays.
"""
function s₁s₂2β₁β₂(s₁,s₂)
    s1 = [s₁...,0.0] # 3-vector
    s2 = [s₂...,0.0] # 3-vector
    fact = 2π / norm(s1 × s2)
    β1 = fact .* s2 × [0,0,1]
    β2 = fact .* [0,0,1] × s1
    β₁ = SV2(β1[1:2])
    β₂ = SV2(β2[1:2])
    return β₁, β₂
end


"""
    check_optional_kw_arguments!(kwargs:: AbstractDict{Symbol,T} where T)

Check the validity of the optional keyword arguments passed to one of the 
user-callable, specific sheet constructor functions.  If any of the arguments
were not passed, assign appropriate default values. 
"""
function check_optional_kw_arguments!(kwargs:: AbstractDict{Symbol,T} where T)
    defaults = Dict(:class=>'J', :dx=>0.0, :dy=>0.0, :rot=>0.0, :Rsheet=>0.0, :save=>"", :fufp=>false)
    validkws = keys(defaults)

    badkws = setdiff(keys(kwargs), validkws)
    if !isempty(badkws)
        error("Illegal keywords: ", join(badkws, ", "))
    end
    for (key,val) in defaults
        haskey(kwargs, key) || (kwargs[key] = val)
    end

    key = :class; val = kwargs[key]
    val ≠ 'J' && val ≠ 'M' && error("Illegal value $val for $key")

    for key in [:dx, :dy, :rot, :Rsheet]
        val = kwargs[key]
        val isa Real || error("$key must be a Real")
        key == :Rsheet && val < 0 && error("$key must be nonnegative")
    end

    key = :save; val = kwargs[key]
    val isa AbstractString || error("$key must be an AbstractString")

    return
end

const optional_kwargs = 
"""
## Optional arguments:
- `class::Char='J'`  Specify the class, either `'J'` or `'M'`.. If `'J'`,  the unknowns are electric surface 
           currents, as used to model a wire or metallic patch-type FSS.  If `'M'`,  the unknowns are
           magnetic surface currents, as used to model a slot or aperture in a perfectly conducting plane.
- `dx::Real=0.0`, `dy::Real=0.0`:  These specify the offsets in the x and y directions applied to the entire 
           unit cell and its contents.  Length units are as specified in the `units` keyword. 
- `rot::Real=0.0`:  Counterclockwise rotation angle in degrees applied to the entire unit cell and its contents. 
           This rotation is applied prior to any offsets specified in `dx` and `dy`.
- `Rsheet::Real=0.0`:  The surface resistance of the FSS conductor in units of Ohm per square.  
            This is only meaningful for a sheet of class `'J'`.
- `fufp::Bool`:  This keyword is not usually required. 
                `fufp` is mnemonic for "Find Unique Face Pairs".  If true, the code will search the 
                triangulation for classes of triangle
                pairs that are the equivalent in the toeplitz sense.  I.e., if triangle pairs (A,B) and (C,D) belong
                to the same equivalence class,  the six vertices in the pair (A,B) can be made to coincide 
                with those of pair (C,D) by a simple translation. If there are many such equivalent pairs, 
                a significant decrease in matrix fill time ensues by exploiting the equivalence.  The tradeoff
                is the time needed to identify them.  The default value is `true` for the `strip` and 
                `meander` styles (those employing structured meshes) and `false` for the remaining styles 
                 (those employing unstructured meshes).
- `save::String=""` Specifies a file name to which the sheet triangulation and unit cell data is to be written,
                   typically to be plotted later.
        
"""
        


"""
    pecsheet()

Return a variable of type `RWGSheet` that contains a perfect electric conducting sheet (i.e. an "E-wall").

"""
function pecsheet()::RWGSheet
    sheet = RWGSheet()
    sheet.style = "NULL"
    sheet.class = 'E'
    return sheet
end # function

"""
    pmcsheet()

Return a variable of type `RWGSheet` that contains a perfect magnetic conducting sheet (i.e. an "H-wall").

"""
function pmcsheet()::RWGSheet
    sheet = RWGSheet()
    sheet.style = "NULL"
    sheet.class = 'H'
    return sheet
end # function




"""
    rectstrip(;Lx::Real, Ly::Real, Nx::Int, Ny::Int, Px::Real, Py::Real, units::PSSFSSLength, kwargs...)

Return a variable of type `RWGSheet` that contains the triangulation for a rectangular strip.

# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `Lx` and `Ly`:  Lengths of the strip in the x and y directions.
- `Px` and `Py`:  Lengths (periods) of the rectangular unit cell in the x and y directions.
- `Nx` and `Ny`:  Number of line segments in the x and y directions, for dividing up the strip into
                  rectangles, which are  triangulated by adding a diagonal to each rectangle.
    
$(optional_kwargs)
"""
function rectstrip(; Lx::Real, Ly::Real, Nx::Int, Ny::Int, Px::Real, Py::Real, units::PSSFSSLength, 
                kwarg...)::RWGSheet
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(Lx); @testpos(Ly); @testpos(Nx); @testpos(Ny); @testpos(Px); @testpos(Py)

    sheet = RWGSheet()
    sheet.style = "rectstrip"
    sheet.units = units

    sheet.s₁ = SV2([Px,0.0])
    sheet.s₂ = SV2([0.0,Py])
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)
        
    nodecount = (Nx+1) * (Ny+1)  # Number of nodes.
    edgecount = 3*Nx*Ny + Nx + Ny  # Number of edges.
    facecount = 2*Nx*Ny  # Number of faces.

    # Setup nodes
    sheet.ρ = SV2.([zeros(2) for i in 1:nodecount])
    x0 = 0.5 * (Px - Lx)  # Center strip in unit cell
    y0 = 0.5 * (Py - Ly)  # Center strip in unit cell
    n = 0  # Initialize node index.
    for j in 0:Ny
        yj = j * (Ly / Ny)
        for i in 0:Nx
            n += 1  # Bump node index.
            sheet.ρ[n] = SV2([x0 + i * (Lx / Nx), y0 + yj])
        end
    end
    
    # Set up the edge matrices:
    sheet.e1 = zeros(Int, edgecount)
    sheet.e2 = zeros(Int, edgecount)
    e = 0  # Initialize edge index.
    # Do the horizontal edges:
    for j in 0:Ny
        kadd = j * (Nx+1)
        for i in 1:Nx
            e += 1
            sheet.e1[e] = i + kadd
            sheet.e2[e] = i + kadd + 1
        end
    end
    # Do the vertical edges:
    for j in 1:Ny
        kadd = (j-1) * (Nx+1) + 1
        for i in 0:Nx
            e += 1 
            sheet.e1[e] = i + kadd
            sheet.e2[e] = i + kadd + (Nx+1)
        end
    end
    # Do the diagonal edges:
    for j in 1:Ny
        kadd1 = (j-1) * (Nx+1)
        kadd2 = 1 + j * (Nx+1)
        for i in 1:Nx
            e += 1
            sheet.e1[e] = i + kadd1
            sheet.e2[e] = i + kadd2
        end
    end

    # Done with edges.  Begin setting up faces
    sheet.fv = zeros(Int, (3,facecount))
    sheet.fe = zeros(Int, (3,facecount))
    sheet.fr = zeros(Float64, facecount)
    Nhe = Nx*Ny + Nx  # Number of horizontal edges.
    Nve = Nx*Ny + Ny  # Number of vertical edges.
    Nde = Nx*Ny       # Number of diagonal edges
    f = 0  # Initialize face index.
    for j in 1:Ny
        nadd1 = (j-1) * (Nx+1)
        nadd2 = 1 + j * (Nx+1)
        for i in 1:Nx
            f += 1  # Bump face index (upper left face).
            sheet.fv[1,f] = i + nadd1  # Lower Left vertex.
            sheet.fv[2,f] = i + nadd2  # Upper right vertex.
            sheet.fv[3,f] = i + nadd2 - 1 # Upper left vertex.
            sheet.fe[1,f] = i + j*Nx  # Upper edge.
            sheet.fe[2,f] = i + (Nhe + nadd1)  # Left edge
            sheet.fe[3,f] = i + (Nhe + Nve + (j-1)*Nx) # Diagonal edge
            f += 1  # Bump face index (lower right face).
            sheet.fv[1,f] = sheet.fv[1,f-1]  # Lower Left vertex.
            sheet.fv[2,f] = 1 + sheet.fv[1,f]  # Lower right vertex.
            sheet.fv[3,f] = sheet.fv[2,f-1] # Upper right vertex.
            sheet.fe[1,f] = 1 + sheet.fe[2,f-1]  # Right edge.
            sheet.fe[2,f] = sheet.fe[3,f-1]  # Diagonal edge.
            sheet.fe[3,f] = sheet.fe[1,f-1] - Nx # Bottom edge
        end
    end
    
    # Done with faces.  Set the face sheet resistance values.
    Rsheet = kwargs[:Rsheet]
    sheet.fr .= Rsheet           # Broadcast value to entire array.

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.,0.]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.ξη_check = (Lx == Px || Ly == Py)

    return sheet

end # function



"""
    polyring(;s1::Vector, s2::Vector, a::Vector, b::Vector, sides::Int ,ntri::Int ,orient::Real, units::PSSFSSLength, kwargs...) --> RWGSheet

Return a variable of type `RWGSheet` that contains the triangulation for one or more concentric annular regions bounded by polygons.

# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `s1` and `s2`:  2-vectors containing the unit cell lattice vectors.
- `a` and `b`:  n-vectors (n>=1) of the same length providing the inner and outer radii, respectively of the polygonal rings.
               Entries in `a` and `b` must be strictly increasing, except for possibly `b[end]` as discussed 
               below. `b[i] > a[i]` ∀ `i ∈ 1:n`, except possibly `b[end]` as discussed below. 
               `a[1]` may be zero to denote a solid (non-annular) polygon as the first "ring".
                It is possible to let the outermost ring to extend completely to the unit cell boundary.  
                This is specified by setting `b[end]` < 0, in which case `-b[end]` is interpreted as the 
                number of edges along the shorter of the `s1` and `s2` lattice vectors.
- `sides`:  The number (>= 3) of polygon sides.
- `ntri`:  The desired total number of triangles distributed among all the annular regions. This is a guide, the actual number 
           will likely be different.
    
$(optional_kwargs)
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to locate the initial
           vertex of the polygonal rings.  The default is to locate the vertex on the
           positive x-axis.

"""
function polyring(;s1::Vector, s2::Vector, a::Vector{<:Real}, b::Vector{<:Real},
                  sides::Int, ntri::Int, units::PSSFSSLength,
                  orient::Real=0.0, kwarg...)::RWGSheet
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    @testpos(sides); @testpos(ntri)
    (length(s1) == length(s2) == 2) || throw(ArgumentError("s1 and s2 must have length 2"))
    
                  
    length(a) ≠ length(b) && throw(ArgumentError("length(a) !== length(b)"))
    nring = length(a)
    for i in 1:nring
        if i < nring || (i == nring && b[nring] > 0)
            a[i] ≥ b[i] && throw(ArgumentError("a[$i] ≥ b[$i]"))
        end
    end 
    for i in 1:nring-1
        b[i] ≥ a[i+1] && throw(ArgumentError("b[$i] ≥ a[$(i+1)]"))
        a[i+i] - a[i] ≤ 0 && throw(ArgumentError("Elements of a must be strictly increasing"))
        b[i+1] - b[i] ≤ 0 && i < nring-1 &&
            throw(ArgumentError("All but final element of b must be strictly increasing"))
    end

    
    if b[nring] < 0 
        fillcell = true  # outer ring extends all the way to unit cell boundaries.
        ns1 = Int(-b[nring])  # number of edges along s1 direction.
        norms1, norms2 = (norm(s1), norm(s2))
        if norms1 ≤ norms2
            ns2 = round(Int, ns1 * norms2 / norms1) # number of edges along s2 direction
        else
            ns2 = ns1  # given value is to be associated with shorter edge.
            ns1 = round(Int, ns2 * norms1 / norms2)
        end
        #ncvrt = [2*(ns1+ns2), sides] # Num of vertices in outer, inner boundary
    else
        fillcell = false # outer ring has finite width.
        #ncvrt = [sides, sides]
    end
    
    ρ₀ = 0.5 * (s1 + s2) # calculate center of polygon.
    α = 360/sides
    
    # compute area of each ring and total area of all rings:
    area_factor = sides/2 * sind(α)
    area = [(b[i]^2 - a[i]^2) * area_factor for i in 1:nring]
    if fillcell  # need to recompute outer ring's area:
        area[nring] = zhatcross(s1) ⋅ s2 - area_factor * a[nring]^2
    end
    areat = sum(area) # total area of all rings.
    areatri = areat/ntri # Desired area of a single triangle

    ρ = Array{SV2}(undef, 0)
    e1 = Array{Cint}(undef, 0)
    e2 = Array{Cint}(undef, 0)
    segmarkers = Array{Cint}(undef, 0)
    holes = Array{SV2}(undef, 0)
    boundary = 0
    node = 0
    for iring in 1:nring
        if a[iring] == 0
            #  solid polygon:
            boundary += 1
            for i = 1:sides
                node += 1
                push!(e1, node)
                push!(e2, node+1)
                push!(segmarkers, boundary)
                push!(ρ, ρ₀ + b[iring] * SV2([reverse(sincosd(orient+(i-1)*α))...]))
            end 
            e2[end] -= sides
        else 
            if iring == nring && fillcell
                # annulus with regular polygon as inner boundary and unit cell as outer:
                # outer boundary first:
                boundary += 1
                i1 = 1
                i2 = ns1 + 1
                nodesave = node + 1 # initial node of outer boundary
                for i in i1:i2  # bottom edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node+1)
                    push!(segmarkers, boundary)
                    push!(ρ, (i-i1)/ns1 * s1)
                end
                i1 = ns1 + 2
                i2 = ns1 + ns2 + 1 
                for i in i1:i2 # right edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node+1)
                    push!(segmarkers, boundary)
                    push!(ρ, s1 + (i-i1+1)/ns2 * s2)
                end
                i1 = ns1 + ns2 + 2 
                i2 = 2*ns1 + ns2 + 1 
                for i in i1:i2 # top edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node+1)
                    push!(segmarkers, boundary)
                    push!(ρ, s1 + s2 - (i-i1+1)/ns1 * s1)
                end
                i1 = 2*ns1 + ns2 + 2
                i2 = 2*ns1 + 2*ns2 
                for i in i1:i2 # left edge
                    node += 1
                    push!(e1, node)
                    push!(e2, node+1)
                    push!(segmarkers, boundary)
                    push!(ρ, s2 - (i-i1+1)/ns2 * s2)
                end
                e2[end] = nodesave
                          
                # now inner boundary:
                boundary += 1
                i1 = 2*(ns1 + ns2) + 1
                i2 = i1 + sides - 1
                nodesave = node+1 # first node of inner boundary
                for i in i1:i2
                    node += 1
                    push!(e1, node)
                    push!(e2, node+1)
                    push!(segmarkers, boundary)
                    push!(ρ,  ρ₀ + a[iring] * SV2([reverse(sincosd(orient+(i-i1)*α))...]))
                end
                e2[end] = nodesave
            else
                # regular polygonal annulus:
                for r in (b[iring], a[iring])
                    boundary += 1
                    nodesave = node + 1
                    for i in 1:sides
                        node += 1
                        push!(e1, node)
                        push!(e2, node+1)
                        push!(segmarkers, boundary)
                        push!(ρ, ρ₀ + r * SV2([reverse(sincosd(orient+(i-1)*α))...]))
                    end
                    e2[end] = nodesave
                end
            end
        end
    end

    # Calculation coordinates of "hole" points
    ρhole = Array{SV2}(undef, 0)
    a[1] > 0 && push!(ρhole, ρ₀)
    
    for i in 1:nring-1
        unitvector = SV2([reverse(sincosd(orient))...])
        r = 0.5 * (b[i] + a[i+1])
        push!(ρhole, ρ₀ + r * unitvector)
    end
    
    # Set up call to meshsub
    points = convert(Matrix{Cdouble}, hcat(ρ...))
    segments = convert(Matrix{Cint}, transpose(hcat(e1, e2)))
    if isempty(ρhole)
        holes = Array{Cdouble}(undef, 2,0)
    else
        holes = convert(Matrix{Cdouble}, hcat(ρhole...))
    end
    sheet = meshsub(points=points, seglist=segments, segmarkers=segmarkers,
                    holes=holes, area=areatri, ntri=ntri)
    
    # Set the face sheet resistance values.
    Rsheet = kwargs[:Rsheet]
    sheet.fr .= Rsheet  # Broadcast value to entire array.

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.,0.]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "polyring"
    sheet.ξη_check = fillcell
    sheet.units = units
    sheet.s₁ = SV2(s1)
    sheet.s₂ = SV2(s2)
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function polyring

"""
    meander(;a::Real, b::Real, h::Real, w1::Real, w2::Real, ntri::Int,
                  units::PSSFSSLength, kwarg...) --> sheet::RWGSheet

# Description:
Return a variable of type `RWGSheet` that contains the triangulation for 
a meanderline strip.  The returned `sheet` has the components `s₁`, `s₂`, 
`β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, and `fr` properly initialized.  
Geometrical parameters are shown in the following diagram:
 
      - - - - - - - - - - - - - - - - - - - - - - - - -             ^
     |                                                |             |
     |                                                |             |
     |                                                |             |
     |                                                |             |
     |                                                |             |
     |            <-------- a/2 ------->              |             |
     |               (center-to-center)               |             |
     |                                                |             |
     |          ----------------------------          |  ^    ^     b
     |          |                          |          |  w2   |     |
     |          |                          |          |  |    |     |
     |          | -----------------------  |          |  v    |     |
     |          | |                     |  |          |             |
     |       -->| |<--w1           w1-->|  |<--       |       h     |
     ------------ |                     |  ------------  ^          |
     |            |                     |             |  w2   |     |
     |            |                     |             |  |    |     |
     ------------ - - - - - - - - - - - ---------------  v    v     v
 
     <-------------------- a ------------------------->
 
 
`a` and `b` are unit cell dimensions.  `w1` and `w2` are the widths
   of the vertical and horizontal strips, resp. `h` is the total
   height of the meander.


# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `a`,`b`,`h`,`w1`, `w2`: Geometrical parameters as defined above.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles distributed among all the annular regions. 
           This is a guide, the actual number will likely be different.
    
$(optional_kwargs)
 """
function meander(;a::Real, b::Real, h::Real, w1::Real, w2::Real, ntri::Int,
                  units::PSSFSSLength, kwarg...)::RWGSheet
                 
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(a); @testpos(b); @testpos(h); @testpos(w1); @testpos(w2); @testpos(ntri) 

    ρ = Array{SV2}(undef, 0)
    e1 = Array{Cint}(undef, 0)
    e2 = Array{Cint}(undef, 0)
    segmarkers = Array{Cint}(undef, 0)
    node = 0

    # Calculate chopping increment using formulas from notes dated
    t1 = w2 * (a + 2w1) / (2 * (h - 2w2)^2)
    t2 = w1 / (h - 2w2)
    fny2 = sqrt(ntri/(4*(t1+t2)))
    ny2 = round(Int, fny2)
    ny2 ≤ 0 &&  (ny2 = 1)
    ny1 = round(Int, fny2 * w2 / (h - 2w2))
    ny1 ≤ 0 && (ny1 = 1)
    nx1 = round(Int, fny2 * (a - 2w1) / (4 * (h - 2w2)))
    nx1 ≤ 0 && (nx1 = 1)
    nx2 = round(Int, fny2 * w1 / (h - 2w2))
    nx2 ≤ 0 && (nx2 = 1)
    ntotal = 4*(2ny1*(nx1+nx2) + nx2*ny2) # Actual # of triangles
    
    # Triangulate first section:
    yoffset = (b - h) / 2
    Lx = (a/2 - w1) / 2
    Ly = h - 2w2
    ρbl = SV2([0.0, yoffset])
    ρtr = ρbl + SV2([Lx, w2])
    sh1 = recttri(ρbl,ρtr,nx1,ny1)
    # Triangulate third section:
    ρbl = SV2([Lx, yoffset])
    ρtr = ρbl + SV2([w1, w2])
    sh2 = recttri(ρbl,ρtr,nx2,ny1)
    # Combine them:
    sh3 = combine(sh1, sh2, 'x', Lx)
    # Add to section 5, store result in sh2:
    ρbl = SV2([Lx, yoffset+w2])
    ρtr = ρbl + SV2([w1, Ly])
    sh1 = recttri(ρbl, ρtr, nx2, ny2)
    sh2 = combine(sh3, sh1, 'y', ρbl[2])
    #  Add to section 7, store result in sh3
    ρbl = SV2([Lx, yoffset+w2+Ly])
    ρtr = ρbl + SV2([w1, w2])
    sh1 = recttri(ρbl, ρtr, nx2, ny1)
    sh3 = combine(sh2, sh1, 'y', ρbl[2])
    #  Add to sections 9 and 10, store result in sh2
    ρbl = SV2([Lx+w1, yoffset+w2+Ly])
    ρtr = ρbl + SV2([2*Lx, w2])
    sh1 = recttri(ρbl, ρtr, 2nx1, ny1)
    sh2 = combine(sh3, sh1, 'x', ρbl[1])
    #  Add to section 8, store result in sh3
    ρbl = SV2([ρtr[1], ρbl[2]])
    ρtr = ρbl + SV2([w1, w2])
    sh1 = recttri(ρbl, ρtr, nx2, ny1)
    sh3 = combine(sh2, sh1, 'x', ρbl[1])
    # Add to section 6, store result in sh2:
    ρbl = SV2([ρbl[1], yoffset+w2])
    ρtr = ρbl + SV2([w1, Ly])
    sh1 = recttri(ρbl, ρtr, nx2, ny2)
    sh2 = combine(sh3, sh1, 'y', ρtr[2])
    #  Add to section 4, store result in sh3
    ρbl = SV2([ρbl[1], yoffset])
    ρtr = ρbl + SV2([w1, w2])
    sh1 = recttri(ρbl, ρtr, nx2, ny1)
    sh3 = combine(sh2, sh1, 'y', ρtr[2])
    # Add to section 2, store result in sh2:
    ρbl = SV2([ρtr[1], yoffset])
    ρtr = ρbl + SV2([Lx, w2])
    sh1 = recttri(ρbl, ρtr, nx1, ny1)
    sheet = combine(sh3, sh1, 'x', ρbl[1])


    # Set the face sheet resistance values.
    sheet.fr = zeros(size(sheet.fv,2))
    Rsheet = kwargs[:Rsheet]
    sheet.fr .= Rsheet  # Broadcast value to entire array.

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    sheet.style = "meander"
    sheet.units = units
    sheet.s₁ = SV2([a, 0.0])
    sheet.s₂ = SV2([0.00, b])
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.,0.]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.ξη_check = true
    
    return sheet

end # function

"""
    loadedcross(;s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
    !
 
# Description:

Create a variable of type `RWGSheet` that
contains the triangulation for a "loaded cross" type of geometry.
The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, 
and `fr` properly initialized.


The following (very poor) "ascii art" attempts to show
the definitions of the geometrical parameters `L1`, `L2` and `w`.
Note that the structure is supposed to be symmetrical wrt reflections
about its horizontal and vertical centerlines, and wrt reflections through a line oriented
at a 45 degree angle wrt the x-axis.


     ^                 ----------------
     |                 |  _________   |
     |                 |  |       |   |
     |                 |  |       |   |
     |                 |  |    -->|   |<--- W
     |                 |  |       |   |
     |                 |  |       |   |
     |     ------------   |       |   -------------
     |     |  |-----------|       |------------|  |
     |     |  |                                |  |
     L1    |  |                                |  |
     |     |  |                                |  |
     |     |  |                                |  |
     |     |  ------------          ------------  |
     |     |-----------   |        |  ------------|
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |        |  |
     |                 |  |________|  |
     |                 |              |
     V                 ----------------
    
                       <---- L2 ------>
    
# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `P`: The period, i.e. the side length of the square unit cell.
- `L1`,`L2`,`w`: Geometrical parameters as defined above.  Note that it is permissible
   to specify `w ≥ L2/2` in which case a solid (i.e., singly-connected) cross will be 
   generated.  In that case the `L2` dimension will be used for the width of the cross pieces.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
           the actual number will likely be different.
    
$(optional_kwargs)
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to locate the initial
           vertex of the polygonal rings.  The default is to locate the vertex on the
           positive x-axis.
"""
function loadedcross(;s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real,
                     ntri::Int, orient::Real=0.0, units::PSSFSSLength, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    @testpos(L1); @testpos(L2); @testpos(w); @testpos(ntri)
    (length(s1) == length(s2) == 2) || throw(ArgumentError("s1 and s2 must be 2-vectors"))
    
    # Initialization:
    nv = (2w < L2 ? 24 : 12) # Total number of vertices
    ρ₀ = 0.5 * (s1 + s2) # Calculate center of polygon.
    ρ = Array{SV2}(undef, nv)
    e1 = Array{Cint}(undef, nv)
    e2 = Array{Cint}(undef, nv)
    segmarkers = Array{Cint}(undef, nv)
    holes = Array{Cdouble}(undef, 2,0)

    # Set up the (outer) polygon geometry:
    ρ[1] = SV2([L2/2, L2/2])
    ρ[2] = SV2([L1/2, ρ[1][2]]) 
    ρ[3] = SV2([ρ[2][1], -ρ[2][2]])
    ρ[4] = SV2([ρ[1][1], ρ[3][2]])
    ρ[5] = SV2([ρ[4][1], -L1/2])
    ρ[6] = SV2([-ρ[5][1], ρ[5][2]])
    ρ[7] = SV2([ρ[6][1], ρ[4][2]])
    ρ[8] = SV2([-ρ[3][1], ρ[7][2]])
    ρ[9] = SV2([ρ[8][1], ρ[1][2]])
    ρ[10] = SV2([ρ[7][1], ρ[9][2]])
    ρ[11] = SV2([ρ[10][1], -ρ[6][2]])
    ρ[12] = SV2([ρ[1][1], ρ[11][2]])
    e1[1:12] = 1:12
    e2[1:11] = 2:12
    e2[12] = 1
    segmarkers[1:12] .= 1
    areat =  2 * (L1 - L2) * L2 # total area for solid cross
    
    if 2w < L2
        #  Set up inner boundary for annulus:
        ρ[13] = ρ[1] .- w
        ρ[14] = ρ[2] .- w
        ρ[15] = ρ[3] + [-w,w]
        ρ[16] = ρ[4] + [-w,w]
        ρ[17] = ρ[5] + [-w,w]
        ρ[18] = ρ[6] .+ w
        ρ[19] = ρ[7] .+ w
        ρ[20] = ρ[8] .+ w
        ρ[21] = ρ[9] + [w,-w]
        ρ[22] = -ρ[16]
        ρ[23] = -ρ[17]
        ρ[24] = -ρ[18]
        e1[13:24] = 13:24
        e2[13:23] = 14:24
        e2[24] = 13
        segmarkers[13:24] .= 2
        holes = [holes ρ₀]
        areat -= 2 * (L1 - L2) * (L2 - 2w)  # Subract inner void
    end      

    if orient ≠ 0
        s,c = sincosd(orient)
        rotmat = SA[c -s;s c]
        for n in eachindex(ρ)
            ρ[n] = rotmat * ρ[n]
        end
    end

    ρ = [t + ρ₀ for t in ρ] # Center on the unit cell
    
    # Set up call to meshsub
    areatri = areat / ntri
    points = convert(Matrix{Cdouble}, hcat(ρ...))
    segments = convert(Matrix{Cint}, transpose(hcat(e1, e2)))
    sheet = meshsub(points=points, seglist=segments, segmarkers=segmarkers,
                    holes=holes, area=areatri, ntri=ntri)
    
    # Set the face sheet resistance values.
    Rsheet = kwargs[:Rsheet]
    sheet.fr .= Rsheet  # Broadcast value to entire array.

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.,0.]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "loadedcross"
    sheet.ξη_check = false
    sheet.units = units
    sheet.s₁ = s1
    sheet.s₂ = s2
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet
                                   
end # function




"""
    jerusalemcross(;P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that
contains the triangulation for a "loaded cross" type of geometry.
The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, 
and `fr` properly initialized.


The following "ascii art" attempts to show
the definitions of the geometrical parameters `P`, `L1`, `L2`, `A`, `B`, and `w`.
Note that the structure is supposed to be symmetrical wrt reflections
about its horizontal and vertical centerlines, and wrt reflections through a line oriented
at a 45 degree angle wrt the x-axis.


    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓ 
    ┃                                                       ┃ _______
    ┃               ┌────────────────────────┐              ┃    ↑
    ┃               │ ┌───────────────────┐  │              ┃    │
    ┃               │ └───────┐    ┌──────┘  │              ┃    │
    ┃               └──────┐  │    │ ┌───────┘              ┃    │
    ┃                      │  │    │ │                      ┃    │
    ┃  ┌───────┐           │  │    │ │            ┌──────┐  ┃    │
    ┃  │  ┌─┐  │           │  │    │ │            │ ┌──┐ │  ┃    │
    ┃  │  │ │  │           │  │   →│ │← w         │ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  │ │  └───────────┘  │    │ └────────────┘ │  │ │  ┃    │
    ┃  │  │ └─────────────────┘    └────────────────┘  │ │  ┃    
    ┃  │  │                                            │ │  ┃   L1 
    ┃  │  │ ┌─────────────────┐    ┌────────────────┐  │ │  ┃  
    ┃  │  │ │  ┌───────────┐  │    │ ┌────────────┐ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  │ │  │           │  │    │ │            │ │  │ │  ┃    │
    ┃  │  └─┘  │          →│  │    │ │← L2     B →│ └──┘ │← ┃    │
    ┃  └───────┘           │  │    │ │            └──────┘  ┃    │
    ┃                      │  │    │ │                      ┃    │
    ┃               ┌──────┘  │    │ └───────┐              ┃    │
    ┃               │ ┌───────┘    └──────┐  │              ┃    │
    ┃               │ └───────────────────┘  │              ┃    │
    ┃               └────────────────────────┘              ┃ ___↓___
    ┃               |<───────── A ──────────>|              ┃
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛ 
    |<─────────────────────── P ───────────────────────────>|
                        
    
    
# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `P`: The period, i.e. the side length of the square unit cell.
- `L1`,`L2`, `A`, `B`, `w`: Geometrical parameters as defined above.  Note that it is permissible
   to specify `w ≥ L2/2` and/or `w ≥ B/2` in which case the respective region will
   be filled in solidly with triangles.  If both conditions hold, then the entire structure will be
   filled in (i.e., singly-connected).  In that case the `L2` and `B` dimensions will be used 
   for the respective widths of the arms, and `w` will not be used.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
           the actual number will likely be different.
    
$(optional_kwargs)
"""
function jerusalemcross(;P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real,
                     ntri::Int, units::PSSFSSLength, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    @testpos(A); @testpos(B); @testpos(L1); @testpos(L2); @testpos(w); @testpos(ntri); @testpos(P)
    

    areaouter =  4 * (A*B + ((L1 - L2)/2 - B) * L2) + L2^2   # outer area for solid cross

    # Total number of vertices and holes and total area:
    if 2w < L2 && 2w < B
        nv = 28 + 28
        nholes = 1
        areat = areaouter - 4*((A-2w)*(B-2w) + (L2-2w)*((L1-L2)/2 - B)) - (L2-2w)^2
    elseif 2w ≥ L2 && 2w ≥ B
        nv = 28
        nholes = 0
        areat = areaouter
    elseif 2w < L2 && 2w ≥ B
        nv = 28 + 12
        nholes = 1
        areat = areaouter - 4*(L2-2w)*((L1-L2)/2 - B) - (L2-2w)^2
    elseif 2w ≥ L2 && 2w < B
        nv = 28 + 4*4
        nholes = 4
        areat = areaouter - 4*(A-2w)*(B-2w) - (L2-2w)^2
    end

    s1 = Cdouble[P, 0.0]; s2 = Cdouble[0.0, P]
    r0 = 0.5 * (s1 + s2) # Calculate center of polygon.
    r = zeros(Cdouble, 2,nv)
    e1 = Array{Cint}(undef, nv)
    e2 = Array{Cint}(undef, nv)
    segmarkers = Array{Cint}(undef, nv)
    holes = Array{Cdouble}(undef, 2,nholes)

    # Set up the (outer) polygon geometry:
    r[:,1] = [L1/2, A/2]
    r[:,2] = [L1/2-B, A/2]
    r[:,3] = [L1/2-B, L2/2]
    r[:,4] = [L2/2, L2/2]
    r[:,5] = [L2/2, L1/2-B]
    r[:,6] = [A/2, L1/2-B]
    r[:,7] = [A/2, L1/2]
    r[:,8:14] = reverse([-1 0; 0 1] * r[:,1:7], dims=2)
    r[:,15:28] = reverse([1 0; 0 -1] * r[:,1:14], dims=2)
    e1[1:28] = 1:28
    e2[1:27] = 2:28
    e2[28] = 1
    segmarkers[1:28] .= 1
    
    if 2w < L2 && 2w < B
        #  Set up inner boundary for annulus:
        r[:,29] = r[:,1] + [-w,-w]
        r[:,30] = r[:,2] + [w,-w]
        r[:,31] = r[:,3] + [w,-w]
        r[:,32] = r[:,4] + [-w,-w]
        r[:,33] = r[:,5] + [-w,w]
        r[:,34] = r[:,6] + [-w,w]
        r[:,35] = r[:,7] + [-w,-w]
        r[:,36:42] = reverse([-1 0; 0 1] * r[:,29:35], dims=2)
        r[:,43:56] = reverse([1 0; 0 -1] * r[:,29:42], dims=2)
        e1[29:56] = 29:56
        e2[29:55] = 30:56
        e2[56] = 29
        segmarkers[29:56] .= 2
        holes[:,1] = [0,0] 
    elseif 2w < L2 && 2w ≥ B
        r[:,29] = [L1/2-B, L2/2-w]
        r[:,30] = [L2/2-w, L2/2-w]
        r[:,31] = [L2/2-w, L1/2-B]
        r[:,32:34] = reverse([-1 0; 0 1] * r[:,29:31], dims=2)
        r[:,35:40] = reverse([1 0; 0 -1] * r[:,29:34], dims=2)
        e1[29:40] = 29:40
        e2[29:39] = 30:40
        e2[40] = 29
        segmarkers[29:40] .= 2
        holes[:,1] = [0,0]
    elseif 2w ≥ L2 && 2w < B
        r[:,29] = [L1/2-B+w, -A/2+w]
        r[:,30] = [L1/2-w, -A/2+w]
        r[:,31] = [L1/2-w, A/2-w]
        r[:,32] = [L1/2-B+w, A/2-w]
        e1[29:32] = 29:32
        e2[29:31] = 30:32
        e2[32] = 29
        segmarkers[29:32] .= 2
        holes[:,1] = [(L1 - B)/2, 0.0]

        r[:,33:36] = [0 1;1 0] * r[:,29:32]
        e1[33:36] = 33:36
        e2[33:35] = 34:36
        e2[36] = 33
        segmarkers[33:36] .= 3
        holes[:,2] = [0.0, (L1 - B)/2]

        r[:,37:40] = [-1 0;0 1] * r[:,29:32]
        e1[37:40] = 37:40
        e2[37:39] = 38:40
        e2[40] = 37
        segmarkers[33:36] .= 4
        holes[:,3] = -holes[:,1]

        r[:,41:44] = [0 1;1 0] * r[:,37:40]
        e1[41:44] = 41:44
        e2[41:43] = 42:44
        e2[44] = 41
        segmarkers[41:44] .= 5
        holes[:,4] = -holes[:,2]
        
    end      



    r .+= r0 # Center on unit cell
    isempty(holes) || (holes .+= r0)
    
    # Set up call to meshsub
    areatri = areat / ntri
    points = r
    segments = convert(Matrix{Cint}, transpose(hcat(e1, e2)))
    sheet = meshsub(points=points, seglist=segments, segmarkers=segmarkers,
                    holes=holes, area=areatri, ntri=ntri)
    
    # Set the face sheet resistance values.
    Rsheet = kwargs[:Rsheet]
    sheet.fr .= Rsheet  # Broadcast value to entire array.

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.,0.]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "jerusalemcross"
    sheet.ξη_check = false
    sheet.units = units
    sheet.s₁ = s1
    sheet.s₂ = s2
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet
                                   
end # function


end # module
