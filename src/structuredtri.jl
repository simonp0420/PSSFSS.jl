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
    function manji(; L1, L2, L3, w, s1, s2, ntri, units, a=0, w2=0, CCW=false, orient=0, kwargs...)
 
# Description:

Create a variable of type `RWGSheet` that
contains the triangulation for a "manji" (Japanese for swastica shape) type of geometry.
The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, 
and `fr` properly initialized.


# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `L1`, `L2`, `L3`, `w`: Geometrical parameters defined in the diagram at
  ![https://simonp0420.github.io/PSSFSS.jl/stable/assets/manjidef.png](https://simonp0420.github.io/PSSFSS.jl/stable/assets/manjidef.png)
  Note that if `a` ≤ `w` then the center square shown in the figure will not be present.  Similarly, if 
  `L3` ≤ `2*w` then the bent portions of the arms will consist of a single strip of width `L3` 
  (without any gap in the middle).
- `s1` and `s2`:  2-vectors containing the unit cell lattice vectors.
- `units`:  Length units (`mm`, `cm`, `inch`, or `mil`)
- `ntri`:  The desired total number of triangles.  This is a guide/request, 
  the actual number will likely be different.
    
## Optional arguments:
- `a::Real=0`: A geometrical parameter defined in the above referenced diagram.  If `a` ≤ `w`
  then the center square shown in that figure will be absent, and the arms will continue uninterrupted
  to the center of the structure.
- `w2::Real=0`: The width of the square ring border between adjacent unit cells, as shown in 
  the above-referenced diagram.  If `w2` ≤ 0 the square ring will not be included in the triangulation.
  Note that `w2 > 0` is only allowed for square unit cells.
- `CCW::Bool=false`: By default, the chiral structure has a clockwise sense.  If 
  `CCW` is `true`, the structure will be counter-clockwise.
- `orient::Real=0.0`:  Counterclockwise rotation angle in degrees applied to the structure within the
  unrotated unit cell.  Note that the outer square ring present when `w2 > 0` will not be rotated by
  use of a nonzero `orient` value.
$(optional_kwargs)
"""
function manji(; L1::Real, L2::Real, L3::Real, a::Real=0, w::Real, w2::Real=0,
    s1::AbstractVector{<:Real}, s2=AbstractVector{<:Real}, CCW::Bool=false, 
    ntri::Int, units::PSSFSSLength, orient::Real=0, kwarg...)
    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    check_optional_kw_arguments!(kwargs)
    @testpos(L1)
    @testpos(L2)
    @testpos(L3)
    @testpos(w)
    @testpos(ntri)
    length(s1) == length(s2) == 2 || error("s1 and s2 must be 2-vectors")
    w2half = w2 / 2
    if w2 > 0 
        s1norm, s2norm = norm.((s1, s2))
        squnitcell = abs(s1 ⋅ s2) / (s1norm * s2norm) < 1e-10 && s1norm ≈ s2norm
        squnitcell || error("w2 > 0 not allowed unless unit cell is square")
        P = s1norm
        Po2 = P / 2
    end
    a < 2L2 || error("a must be less than 2*L2")
    ρ₀ = 0.5 * (s1 + s2) # Calculate center of polygon.
    wo2, ao2 = (w, a) ./ 2
    armsfolded = 2w < L3
    centersquare = a > w
    # Triangulation of the upper arm:
    # unique x and y vertices for arm, total area, and armarea:
    ymin = max(ao2, wo2)
    ymino2 = ymin / 2
    if armsfolded
        xrequired = [-wo2, wo2, L1-3wo2, L1-wo2]
        yrequired = [ymin, L2, L2+w, L2+L3-w, L2+L3]
        areaarm = ((L2 - ymin) + 2*L1 + (L3 - 2*w)) * w
    else
        xrequired = [-wo2, wo2, L1-wo2]
        yrequired = [ymin, L2, L2+L3]
        areaarm = (L2 - ymin) * w + L1 * L3
    end
    areat = ymin^2 + 4 * areaarm
    if w2 > 0 
        #push!(xrequired, Po2-w2half, Po2)
        areat += 4 * (P - w2half) * w2half
    end

    function arm_inside(x::Real, y::Real)
        # predicate to determine if a point is inside upper arm
        if ymin ≤ y < L2
            abs(x) ≤ wo2 && return true
        elseif L2 ≤ y ≤ L2 + w
            -wo2 ≤ x < L1 - wo2 && return true
        elseif L2 + w ≤ y ≤ L2 + L3 - w
            if armsfolded
               L1 - 3*wo2 ≤ x ≤ L1 - wo2 && return true
            else
                -wo2 ≤ x ≤ L1 - wo2 && return true
            end
        elseif L2 + L3 - w ≤ y ≤ L2 + L3
            -wo2 < x < L1 - wo2 && return true
        end
        return false
    end

    # Triangulate upper arm
    ntriarm = round(Int, ntri * areaarm / areat)
    arm = make_plaid_mesh(xrequired, yrequired, areaarm, ntriarm, arm_inside)

    # Triangulate center square
    xrsquare = [x for (x,y) in arm.ρ if y == ymin]
    if centersquare
        xrsquare = sort!(push!(xrsquare, -ao2, ao2))
    end
    areasquare = ymin^2
    ntrisquare = round(ntri * areasquare / areat)
    square = make_plaid_mesh(xrsquare, xrsquare, areasquare, ntrisquare, (x,y) -> true)

    # Combine rotated arms and square
    sheet = square
    rotmat = SA[0.0 -1.0; 1.0 0.0]
    for (coord, val) in zip(('y', 'x', 'y', 'x'), (ymin, -ymin, -ymin, ymin))
        sheet = combine(sheet, arm, coord, val)
        for n in eachindex(arm.ρ)
            arm.ρ[n] = rotmat * arm.ρ[n]
        end
    end

    if CCW
        flipmat = SA[-1.0 0.0; 0.0 1.0]
        for n in eachindex(sheet.ρ)
            sheet.ρ[n] = flipmat * sheet.ρ[n]
        end
    end

    if w2 > 0
        # Add outer ring
        xr = [-Po2, -Po2+w2half, Po2-w2half, Po2]
        yr = xr
        oring_inside(x,y) = Po2-w2half ≤ abs(x) ≤ Po2 || Po2-w2half ≤ abs(y) ≤ Po2 
        arearing = 4 * (P * w2half - w2half^2)
        ntriring = round(Int, ntri * arearing / areat)
        ring = make_plaid_mesh(xr, yr, arearing, ntriring, oring_inside)
        sheet = combine(sheet, ring, ' ', Inf)
    end

    # Rotate, then center sheet on unit cell center
    s, c = sincosd(orient)
    rotmat = SA[c -s; s c]
    for n in eachindex(sheet.ρ)
        x, y = sheet.ρ[n]
        if w2 > 0 && (Po2 - w2half ≤ abs(x) ≤ Po2 || Po2 - w2half ≤ abs(y) ≤ Po2)
                sheet.ρ[n] = sheet.ρ[n] + ρ₀ # Don't rotate points in outer loop
        else
            sheet.ρ[n] = rotmat * sheet.ρ[n] + ρ₀
        end
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

    sheet.style = "manji"
    sheet.ξη_check = w2 > 0
    sheet.units = units
    sheet.s₁ = s1
    sheet.s₂ = s2
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function manji

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
    xr, yr = sort.((xr, yr))
    bigarea = (xr[end] - xr[begin]) * (yr[end] - yr[begin]) # area of circumscribing rectangle
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
