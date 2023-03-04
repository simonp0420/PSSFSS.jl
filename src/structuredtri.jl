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
        xrequired = [-L1o2, -L2o2, L2o2, L1o2]
        areat = (2L1 - L2) * L2 # total area for solid cross
    else
        xrequired = [-L1o2, -L1o2 + w, -L2o2, -L2o2 + w, L2o2 - w, L2o2, L1o2 - w, L1o2]
        areat = (2L1 - L2) * L2  - (2*(L1-2w) - (L2-2w)) * (L2-2w) # total area for loaded cross
    end
    yrequired = copy(xrequired)

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
    sheet = make_plaid_mesh(xrequired, yrequired, areat, ntri, is_inside)

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
        xrequired = [L2o2-w, L2o2, Ao2-w, Ao2, L1o2-B, L1o2-B+w, L1o2-w, L1o2]
        xrequired = vcat(-1*reverse(xrequired), xrequired)
        areat = areaouter - 4 * ((A - 2w) * (B - 2w) + (L2 - 2w) * ((L1 - L2) / 2 - B)) - (L2 - 2w)^2
    elseif allfilled
        xrequired = [L2o2, Ao2, L1o2-B, L1o2]
        xrequired = vcat(-1*reverse(xrequired), xrequired)
        areat = areaouter
    elseif armsloaded && endsfilled
        xrequired = [L2o2-w, L2o2, Ao2, L1o2-B, L1o2]
        xrequired = vcat(-1*reverse(xrequired), xrequired)
        areat = areaouter - 4 * (L2 - 2w) * ((L1 - L2) / 2 - B) - (L2 - 2w)^2
    elseif armsfilled && endsloaded
        xrequired = [L2o2, Ao2-w, Ao2, L1o2-B, L1o2-B+w, L1o2-w, L1o2]
        xrequired = vcat(-1*reverse(xrequired), xrequired)
        areat = areaouter - 4 * (A - 2w) * (B - 2w) - (L2 - 2w)^2
    end
    yrequired = copy(xrequired)

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
    sheet = make_plaid_mesh(xrequired, yrequired, areat, ntri, is_inside)

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


function polyring_structured(; s1::Vector, s2::Vector, a::Vector{<:Real}, b::Vector{<:Real},
    sides::Int, ntri::Int, units::PSSFSSLength,
    orient::Real=0.0, kwarg...)::RWGSheet
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = false)
    check_optional_kw_arguments!(kwargs)
    @testpos(sides)
    @testpos(ntri)
    (length(s1) == length(s2) == 2) || throw(ArgumentError("s1 and s2 must have length 2"))
    sides == 4 || throw(ArgumentError("Structured mesh requires exactly 4 sides"))

    length(a) ≠ length(b) && throw(ArgumentError("length(a) !== length(b)"))
    nring = length(a)
    for i in 1:nring
        if i < nring || (i == nring && b[nring] > 0)
            a[i] ≥ b[i] && throw(ArgumentError("a[$i] ≥ b[$i]"))
        end
    end
    for i in 1:nring-1
        b[i] ≥ a[i+1] && throw(ArgumentError("b[$i] ≥ a[$(i+1)]"))
        a[i+1] - a[i] ≤ 0 && throw(ArgumentError("Elements of a must be strictly increasing"))
        b[i+1] - b[i] ≤ 0 && i < nring - 1 &&
            throw(ArgumentError("All but final element of b must be strictly increasing"))
    end

    fillin = iszero(a[begin])
    if b[nring] < 0
        fillcell = true  # outer ring extends all the way to unit cell boundaries.
        iszero(s1 ⋅ s2) || throw(ArgumentError("structured mesh not possible for nonrectangular unit cell with b[end]<0"))
        testit = (orient - 45) / 90
        round(Int, testit) == testit || error("b[end]<0 requires orient=±45 for structured mesh")
    else
        fillcell = false # outer ring has finite width.
    end

    ρ₀ = 0.5 * (s1 + s2) # calculate center of polygon.

    # compute area of each ring and total area of all rings:
    area = [(b[i]^2 - a[i]^2) * 2 for i in 1:nring]
    if fillcell  # need to recompute outer ring's area:
        area[nring] = norm(s1) * norm(s2) - 2 * a[nring]^2
    end
    areat = sum(area) # total area of all rings.

    xrequired = vec([transpose(a) ; transpose(b)]) * inv(sqrt(2))
    fillin && deleteat!(xrequired, firstindex(xrequired))
    fillcell && (xrequired[end] = 0.5 * norm(s1))
    yrequired = copy(xrequired)
    fillcell && (yrequired[end] = 0.5 * norm(s2))
    xrequired = vcat(-1*reverse(xrequired), xrequired)
    yrequired = vcat(-1*reverse(yrequired), yrequired)

    function is_inside(x::Real, y::Real)
        # predicate to determine if a point is within the region to be triangulated
        x, y = abs.((x,y)) # Due to left/right and up/down symmetry
        y > x && ((x, y) = (y, x)) # Due to symmetry about line x = y
        # (x,y) is now in the region 0 ≤ ϕ ≤ π/4

        xroot2 = √2 * x

        for i in eachindex(a, b)
            a[i] ≤ xroot2 ≤ b[i] && return true
        end
        fillcell && xroot2 ≥ a[end] && return true
        return false
    end


    # Triangulate prior to rotating the orientation
    sheet = make_plaid_mesh(xrequired, yrequired, areat, ntri, is_inside)

    # Rotate, then center sheet on unit cell center
    s, c = sincosd(orient + 45)
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

    sheet.style = "polyring"
    sheet.ξη_check = fillcell
    sheet.units = units
    sheet.s₁ = SV2(s1)
    sheet.s₂ = SV2(s2)
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function polyring_structured

"""
    make_plaid_mesh(xr, yr, area, ntri, is_inside) -> sheet::RWGSheet

Generate a structured, plaid triangular mesh from list of required coordinates and predicate function

# Input Arguments
- `xr`, `yr`: Vectors of required x and y coordinates for vertices of the geometry to be meshed.
- `area`:  The area of the geometry to be meshed.
- `ntri`:  The desired number of triangles for the area to be meshed.
- `is_inside`: A predicate function where `is_inside(x,y) -> tf::Bool` determines whether a point (x,y)
  is within the region to be meshed.

#  Return value:
- `sheet`: A variable of type RWGSheet with fields ρ, e1, e2, fe, and fv properly initialized. The
  mesh results from a plaid rectangular tesselation containing at least the vertices in the Cartesian
  product of `xr` and `yr`, the latter supplemented with additional points to refine the mesh, and then
  converted to a triangular tesselation by adding a diagonal to each rectangle.
"""
function make_plaid_mesh(xr::AbstractVector, yr::AbstractVector, area, ntri, is_inside)::RWGSheet
    length(xr) == length(yr) || error("xr and yu not same length")
    xr, yr = sort.((xr, yr))
    bigarea = (xr[end] - xr[1]) * (yr[end] - yr[1]) # area of circumscribing rectangle
    bignsq = ceil(Int, bigarea / area * ntri/2) # desired number of squares to form in circumscribing rectangle
    s = sqrt(bigarea / bignsq) # ideal side length for squares used to tesselate the big area

    facevs = Tuple{Tuple{Int,Int}, Tuple{Int,Int}, Tuple{Int,Int}}[]
    edgevs = Tuple{Tuple{Int,Int}, Tuple{Int,Int}}[]
    xn = xr[begin:begin]
    yn = yr[begin:begin]
    for (tr, tn) in ((xr, xn), (yr, yn)), i in eachindex(tr)[begin+1:end]
        dt = tr[i] - tr[i-1]
        nt = max(1, round(Int, dt / s))
        append!(tn, tr[i-1] .+ collect((1:nt) * (dt / nt)))
        tn[end] = tr[i] # correct rounding errors
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
