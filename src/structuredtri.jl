"""
    loadedcross_structured(;s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that contains the triangulation for a "loaded cross" 
type of geometry, using a structured mesh. The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, 
`ρ`, `e1`, `e2`, `fv`, `fe`, and `fr` properly initialized.


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
- `s1` and `s2`:  2-vectors containing the unit cell lattice vectors.
- `L1`,`L2`,`w`: Geometrical parameters as defined above.  Note that it is permissible
   to specify `w ≥ L2/2` in which case a solid (i.e., singly-connected) cross will be 
   generated.  In that case the `L2` dimension will be used for the width of the cross pieces.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
           the actual number will likely be different.
    
$(optional_kwargs)
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees used to locate the initial
           vertex of the loaded cross.  The default is to locate the vertex on the
           positive x-axis.
"""
function loadedcross_structured(; s1::Vector{<:Real}, s2::Vector{<:Real}, L1::Real, L2::Real, w::Real,
    ntri::Int, orient::Real=0.0, units::PSSFSSLength, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(L1)
    @testpos(L2)
    @testpos(w)
    @testpos(ntri)
    (length(s1) == length(s2) == 2) || throw(ArgumentError("s1 and s2 must be 2-vectors"))

    # Initialization:
    ρ₀ = 0.5 * (s1 + s2) # Calculate center of polygon.
    L1o2, L2o2 = (L1, L2) ./ 2
    if 2w > L2
        xunique = [-L1o2, -L2o2, L2o2, L1o2]
        areat = (2L1 - L2) * L2 # total area for solid cross
    else
        xunique = [-L1o2, -L1o2 + w, -L2o2, -L2o2 + w, L2o2 - w, L2o2, L1o2 - w, L1o2]
        areat = (2L1 - L2) * L2  - (2*(L1-2w) - (L2-2w)) * (L2-2w) # total area for loaded cross
    end
    yunique = copy(xunique)

    function is_inside(x::Real, y::Real)
        # predicate to determine if a point is within the region to be triangulated
        x, y = abs.((x,y)) # Due to left/right and up/down symmetry
        y > x && ((x, y) = (y, x)) # Due to symmetry about line x = y
        (x > L1o2 || y > L2o2) && return false
        if w ≥ L2o2 # solid
            return true
        else # "loaded"
            (L2/2 - w ≤ y || x ≥ L1/2 - w) && return true
        end
        return false
    end

    # Triangulate prior to rotating the orientation
    sheet = make_plaid_mesh(xunique, yunique, areat, ntri, is_inside)

    # Rotate, then center sheet on unit cell center
    s, c = sincosd(orient)
    rotmat = SA[c -s; s c]
    for n in eachindex(sheet.ρ)
        sheet.ρ[n] = rotmat * sheet.ρ[n] + ρ₀
    end


    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
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
    jerusalemcross_structured(;P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real, 
                 ntri::Int, units::PSSFSSLength, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that contains the triangulation for a 
"Jerusalem cross" type of geometry, using a structured mesh.
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
function jerusalemcross_structured(; P::Real, L1::Real, L2::Real, A::Real, B::Real, w::Real,
    ntri::Int, units::PSSFSSLength, orient::Real=0.0, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(A)
    @testpos(B)
    @testpos(L1)
    @testpos(L2)
    @testpos(w)
    @testpos(ntri)
    @testpos(P)

    s1 = SV2([P, 0])
    s2 = SV2([0, P])
    ρ₀ = 0.5 * (s1 + s2) # Calculate center of polygon.
    L1o2, L2o2, Ao2, Bo2 = (L1, L2, A, B) ./ 2
    areaouter = 4 * (A * B + ((L1 - L2) / 2 - B) * L2) + L2^2   # outer area for solid cross

    armsloaded = 2w < L2
    armsfilled = !armsloaded
    endsloaded = 2w < B
    endsfilled = !endsloaded
    allfilled = armsfilled && endsfilled
    allloaded = armsloaded && endsloaded
    # unique x vertices and total area:
    if allloaded
        xunique = [L2o2-w, L2o2, Ao2-w, Ao2, L1o2-B, L1o2-B+w, L1o2-w, L1o2]
        xunique = vcat(-1*reverse(xunique), xunique)
        areat = areaouter - 4 * ((A - 2w) * (B - 2w) + (L2 - 2w) * ((L1 - L2) / 2 - B)) - (L2 - 2w)^2
    elseif allfilled
        xunique = [L2o2, Ao2, L1o2-B, L1o2]
        xunique = vcat(-1*reverse(xunique), xunique)
        areat = areaouter
    elseif armsloaded && endsfilled
        xunique = [L2o2-w, L2o2, Ao2, L1o2-B, L1o2]
        xunique = vcat(-1*reverse(xunique), xunique)
        areat = areaouter - 4 * (L2 - 2w) * ((L1 - L2) / 2 - B) - (L2 - 2w)^2
    elseif armsfilled && endsloaded
        xunique = [L2o2, Ao2-w, Ao2, L1o2-B, L1o2-B+w, L1o2-w, L1o2]
        xunique = vcat(-1*reverse(xunique), xunique)
        areat = areaouter - 4 * (A - 2w) * (B - 2w) - (L2 - 2w)^2
    end
    yunique = copy(xunique)

    function is_inside(x::Real, y::Real)
        # predicate to determine if a point is within the region to be triangulated
        x, y = abs.((x,y)) # Due to left/right and up/down symmetry
        y > x && ((x, y) = (y, x)) # Due to symmetry about line x = y
        # (x,y) is now in the region 0 ≤ ϕ ≤ π/4
        x < L1o2 - B && y > Ao2 && return false
        x ≥ L1o2 && y > L2o2 && return false

        if allloaded
            L2o2-w ≤ x ≤ L1o2 - B && L2o2 - w ≤ y ≤ L2o2 && return true
            L1o2-B ≤ x ≤ L1o2-B+w && L2o2 - w ≤ y ≤ Ao2 && return true
            L1o2-B+w ≤ x ≤ L1o2-w && Ao2-w ≤ y ≤ Ao2 && return true
            L1o2-w ≤ x ≤ L1o2 && y ≤ Ao2 && return true
            return false
        elseif allfilled
            y ≤ L2o2 && return true
            x > L1o2 - B && y < Ao2 && return true
            return false
        elseif armsloaded && endsfilled
            L2o2-w ≤ x ≤ L1o2 - B && L2o2 - w ≤ y ≤ L2o2 && return true
            L1o2-B ≤ x ≤ L1o2 && y ≤ Ao2 && return true
            return false
        elseif armsfilled && endsloaded
            x ≤ L1o2-B+w && y ≤ L2o2 && return true
            L1o2-B ≤ x ≤ L1o2-B+w && L2o2 - w ≤ y ≤ Ao2 && return true
            L1o2-B+w ≤ x ≤ L1o2-w && Ao2-w ≤ y ≤ Ao2 && return true
            L1o2-w ≤ x ≤ L1o2 && y ≤ Ao2 && return true
            return false
        end
        error("This shouldn't happen")
    end

    # Triangulate prior to rotating the orientation
    sheet = make_plaid_mesh(xunique, yunique, areat, ntri, is_inside)

    # Rotate, then center sheet on unit cell center
    s, c = sincosd(orient)
    rotmat = SA[c -s; s c]
    for n in eachindex(sheet.ρ)
        sheet.ρ[n] = rotmat * sheet.ρ[n] + ρ₀
    end

    sheet.Zs = kwargs[:Zsheet]
    sheet.σ = kwargs[:σ]
    sheet.Rq = kwargs[:Rq]
    sheet.disttype = kwargs[:disttype]

    # Handle remaining optional arguments
    sheet.fufp = kwargs[:fufp]
    sheet.class = kwargs[:class]
    rotate!(sheet, kwargs[:rot])
    dxdy = SV2([kwargs[:dx], kwargs[:dy]])
    if dxdy ≠ [0.0, 0.0]
        sheet.ρ .= (dxdy + xy for xy in sheet.ρ)
    end

    sheet.style = "jerusalemcross"
    sheet.ξη_check = false
    sheet.units = units
    sheet.s₁ = SV2([P, 0])
    sheet.s₂ = SV2([0, P])
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function


"""
    make_plaid_mesh(xu, yu, area, ntri, is_inside) -> sheet::RWGSheet

Generate a structured, plaid triangular mesh from list of unique coordinates and predicate function

# Input Arguments
- `xu`, `yu`: Vectors of unique x and y coordinates for vertices of the geometry to be meshed.
- `area`:  The area of the geometry to be meshed.
- `ntri`:  The desired number of triangles for the area to be meshed.
- `is_inside`: A predicate function where `is_inside(x,y) -> tf::Bool` determines whether a point (x,y)
  is within the region to be meshed.

#  Return value:
sheet      A variable of type RWGSheet with fields ρ, e1, e2, fe, and fv properly initialized.
"""
function make_plaid_mesh(xu::AbstractVector, yu::AbstractVector, area, ntri, is_inside)::RWGSheet
    length(xu) == length(yu) || error("xu and yu not same length")
    xu, yu = sort.((xu, yu))
    bigarea = (xu[end] - xu[1]) * (yu[end] - yu[1]) # area of circumscribing rectangle
    bignsq = ceil(Int, bigarea / area * ntri/2) # desired number of squares to form in circumscribing rectangle
    s = sqrt(bigarea / bignsq) # ideal side length for squares used to tesselate the big area

    facevs = Tuple{Tuple{Int,Int}, Tuple{Int,Int}, Tuple{Int,Int}}[]
    edgevs = Tuple{Tuple{Int,Int}, Tuple{Int,Int}}[]
    xn = xu[begin:begin]
    yn = yu[begin:begin]
    for (tu, tn) in ((xu, xn), (yu, yn)), i in eachindex(tu)[begin+1:end]
        dt = tu[i] - tu[i-1]
        nt = max(1, round(Int, dt / s))
        append!(tn, tu[i-1] .+ collect((1:nt) * (dt / nt)))
        tn[end] = tu[i] # correct rounding errors
    end
    
    # xn and yn now contain the plaid vertex coordinates
    # Initialize vectors of face and edge indices into xn and yn:
    facevs = Tuple{Tuple{Int,Int}, Tuple{Int,Int}, Tuple{Int,Int}}[]
    edgevs = Tuple{Tuple{Int,Int}, Tuple{Int,Int}}[]

    # Add triangular faces and edges within the desired geometry.  Assumes that 
    # if the center of the rectangle is inside, then so are both triangles
    # partitioning that rectangle.
    for i in eachindex(xn)[begin+1:end]
        xcen = 0.5 * (xn[i] + xn[i-1])
        for j in eachindex(yn)[begin+1:end]
            ycen = 0.5 * (yn[j] + yn[j-1])
            is_inside(xcen, ycen) || continue
            topface = ((i,j-1), (i,j), (i-1,j))
            botface = ((i,j-1), (i-1,j), (i-1,j-1))
            push!(facevs, topface, botface)
            push!(edgevs, ((i-1,j-1), (i,j-1)),
                          ((i,j-1), (i,j)),
                          ((i,j), (i-1,j)),
                          ((i-1,j),(i-1,j-1)),
                          ((i,j-1), (i-1,j)))
        end
    end

    nodes = unique!([s[i] for s in facevs for i in 1:3]) # list of node (ix,iy) values
    nface = length(facevs)
    inodes = Dict(t => i for (i,t) in pairs(nodes)) # returns linear node index given (ix,iy)
    edgenodes = unique!([extrema((inodes[e[1]], inodes[e[2]])) for e in edgevs]) # List of edge (m,n) node indices with m < n
    iedges = Dict(t => i for (i,t) in pairs(edgenodes)) # returns linear edge index given (m,n) node indices with m < n
    e1 = first.(edgenodes)
    e2 = last.(edgenodes)
    fv = [inodes[f[k]] for k in 1:3, f in facevs] # face-vertex matrix
    fe = zeros(Int, 3, nface)
    previ = (3, 1, 2)
    nexti = (2, 3, 1)
    for (jf, f) in pairs(facevs), i in 1:3
        edgenodes = (inodes[f[nexti[i]]], inodes[f[previ[i]]])
        fe[i, jf] = iedges[extrema(edgenodes)]
    end

    ρ = [SV2([xn[ix], yn[iy]]) for (ix, iy) in nodes]

    sh = RWGSheet()
    sh.ρ = ρ
    sh.e1 = e1
    sh.e2 = e2
    sh.fv = fv
    sh.fe = fe
    return sh
end



