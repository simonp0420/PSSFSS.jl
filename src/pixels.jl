# Note: this file is included by Elements.jl

"""
    SymPixels

A struct containing the pixel parameters of a symmetric pixel mesh element. Fields:
    - `halfnint`: Half the number of pixels across the interior (non-rim) portion of
      the unit cell. Must be a positive integer.
    - `nrim`: Width of the rim in pixels. Must be a nonnegative integer.
"""
struct SymPixels{T1<:Integer, T2<:Integer}
    halfnint::T1
    nrim::T2      # Width of the rim in pixels ≥ 0
    function SymPixels(halfint::T1, nrim::T2) where {T1<:Integer, T2<:Integer}
        halfint < 1 && throw(ArgumentError("SymPixels first argument must be ≥ 1"))
        nrim < 0 && throw(ArgumentError("SymPixels second argument must be ≥ 0"))
        return new{T1, T2}(halfint, nrim)
    end
end

SymPixels(; halfnint::Integer, nrim::Integer=0) = SymPixels(halfnint, nrim) # keyword constructor

ucsidelen(s::SymPixels) = 2 * (s.halfnint + s.nrim) # Total number of pixels along side of unit cell

"""
    patternveclen(s::SymPixels) -> vlen::Integer


Return the length of the 1/0 pattern vector needed to define the unit cell metallization pattern
for a symmetric pixel mesh.  This is the number of pixels in the irreducible zone for the mesh.
"""
patternveclen(s::SymPixels) = s.halfnint * (1 + s.halfnint) ÷ 2

"""
    sympixmat(s::SymPixels, v::AbstractVector) -> mat::Matrix{Bool}

Compute the metalization pixel pattern matrix for a symmetric pixel mesh.

The matrix will be square of dimension `ucsidelen(s)`.  `true` elements correspond
to metalized pixels in the unit cell.  Element (i,j) of the matrix corresponds
to to the pixel in the unit cell located at the i'th x location and j'th y location,
counting from the lower-left corner of the unit cell.

## Input Arguments
- `s::SymPixels`: Defines the pixel parameters of the mesh.
- `v::AbstractVector`: This is the metallization pattern vector. It must contain only ones and
  zeros, and its length should be equal to the number of pixels in the
  irreducible zone of the mesh.  A greater length is also allowed, in which case a warning is
  issued and the extra elements in `v` are ignored.

"""
function sympixmat(s::SymPixels, v::AbstractVector)
    niz = patternveclen(s)
    length(v) ≥ niz || throw(ArgumentError("v is not long enough"))
    all(x -> isone(x) || iszero(x), v) || throw(ArgumentError("v must consist only of 1s and 0s"))
    length(v) > niz && @warn "v is longer than necessary.  Ignoring extra entries"

    totsidelen = ucsidelen(s)
    mat = falses(totsidelen, totsidelen)

    # Fill the rim
    nrim = s.nrim
    if nrim > 0
        mat[1:nrim, :] .= true
        mat[:, 1:nrim] .= true
        mat[:, end-nrim+1:end] .= true
        mat[end-nrim+1:end, :] .= true
    end

    nir = patternveclen(s) # Number of pixels in irreducible zone
    halfnint = s.halfnint
    nint = 2 * halfnint

    # Work on the interior
    mint = @view mat[1+nrim:end-nrim, 1+nrim:end-nrim] # interior matrix

    # Fill upper right quadrant
    mintur = @view mint[1:halfnint, halfnint+1:end] # upper-right quadrant of interior matrix
    for (i,j,k) in zip(1:nir, ulindices(halfnint), lrindices(halfnint))
        mintur[j] = mintur[k] = v[i]
    end
    # Copy it into the upper left quadrant
    mintul = @view mint[1:halfnint, halfnint:-1:1] # Reverse the columns
    for col in 1:s.halfnint
        mintul[:, col] .= @view mintur[:, col]
    end
    # Copy into the bottom half:
    mintbot = @view mint[nint:-1:halfnint+1, :] # reverse the rows
    for row in 1:halfnint
        mintbot[row, :] .= @view mint[row, :]
    end

    return mat
end

"""
    ulindices(n::Integer)

Return a CartesianIndex generator that iterates over the indices in the upper left half
of a square matrix of order n.  The upper left half is bounded by (and includes) the
antidiagonal of the matrix.  The elements are iterated in column major order.
"""
function ulindices(n::Integer)
    n > 0 || throw(ArgumentError("uladindices requires a a positive argument"))
    return (CartesianIndex(i,j) for (j,k) in zip(1:n, n:-1:1) for i in 1:n if i ≤ k)
end

"""
    lrindices(n::Integer)

Return a CartesianIndex generator that iterates over the indices in the lower right half
of a square matrix of order n.  The lower right half is bounded by (and includes) the
antidiagonal of the matrix.  The elements are iterated in row major order, but with the
elements of each row reversed. This places the elements in 1:1 correspondence with their
mirror images in the upper left half.
"""
function lrindices(n::Integer)
    n > 0 || throw(ArgumentError("lradindices requires a a positive argument"))
    return (CartesianIndex(n-j+1, n-i+1) for (j,k) in zip(1:n, n:-1:1) for i in 1:n if i ≤ k)
end

const keepstart = first(findfirst("- `dx", optional_kwargs)) - 1
const pix_optional_kwargs =
"""
- `class::Char='M'`  Specify the class, either `'J'` or `'M'`. If `'J'`, the unknowns are electric surface
  currents in the areas corresponding to `1` values of the pixels.  If `'M'`,  the unknowns are magnetic surface
  currents in the areas corresponding to `0` values of the pixels.  It is known that using `'J'` can result in
  grossly incorrect results for some geometries where adjacent metallic pixels intersect at only a single point.
  Therefore, use of only `'M'` is strongly recommended for most `pixels` elements and all `sympixels` elements.
""" *
optional_kwargs[keepstart:end]


"""
    function sympixels(; P, nrim, halfnint, patternvec, units, kwargs...)

# Description:

Create a variable of type `RWGSheet` that contains the triangulation for a symmetrically pixelated
square unit cell. The pattern of metallic pixels has C4 symmetry, as well as left-right and up-down mirror symmetry,
as well as each quadrant exhibiting antidiagonal mirror symmetry. The returned value has fields `s₁`, `s₂`, `β₁`,
`β₂`, `ρ`, `e1`, `e2`, `fv`, `fe`, and `fr` properly initialized.  The pixels included in the triangulation are
determined by the `patternvec` input vector as described below.


# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `P::Real > 0`: The side length of the square unit cell, specified in units defined by the `units` keyword argument.
- `nrim::Integer`: The width of the solid metallic rim placed just inside the unit cell boundary, in units of pixels.
  A value of `0` implies that there is no rim.
- `halfnint`: Half the number of pixels in the side length for the interior (strictly inside the rim) square pixelated
  region of the unit cell.
- `patternvec::AbstractVector{<:Integer}`: A vector of length `halfnint*(halfnint+1)÷2`, consisting solely of 1's and 0's.
  The elements of this vector are mapped to pixels in the irreducible zone of the unit cell as shown in the
  diagram at ![https://simonp0420.github.io/PSSFSS.jl/stable/assets/sympixel_with_irzone_numbering.png]\
  (https://simonp0420.github.io/PSSFSS.jl/stable/assets/sympixel_with_irzone_numbering.png).
  Within the irreducible zone, pixels corresponding to a value of `1` (or `true`)
  are taken to be areas of metallization, while `0` or `false` values are metal-free (void) areas.  This holds for either
  `J` or `M` as the `class` value (see the `class` argument below for important limitations).
- `units`:  Length units for `P` (either `mm`, `cm`, `inch`, or `mil`).

## Optional arguments:
- `pdiv::Int = 1`: The number of "chops" or subdivisions applied to each square pixel side when forming the triangulation.
  A value of `1` (the default) means that the pixels included in the triangulation (`1` or `true` values for `class='J'`,
  `0` or `false` values for `class='M'`) are not subdivided any further, except for a single diagonal across each square
  pixel to form triangles. A value of `n>1` means that each square pixel is first divided into `n×n` square
  subpixels, after which a single diagonal edge (if `quad=false`) is added to each subpixel to form triangles.
- `quad::Bool=false`:  If `true`, each subpixel (or pixel, if `pdiv` is 1) is divided into four triangles by adding
  two diagonals.  If `false` (the default), only a single diagonal is added to each square to produce two triangles.
- `sym::Bool = false`: If true, the diagonals added to the squares will exhibit the same left-right and up-down
  mirror symmetry as the collection of `true` (`false`) pixel locations.  If `false` and `quad=false`, then all added diagonals
  across the unit cell will have the same orientation.  `sym` has no effect (i.e. is redundant) if `quad=true` since
  in that case the two added diagonals already ensure mirror symmetry of the triangulation.
$(pix_optional_kwargs)
"""
function sympixels(; P::Real, nrim::Integer, halfnint::Integer, class::Char='M',
    patternvec::AbstractVector{<:Integer}, units, pdiv::Integer=1, quad::Bool=false,
    kwargs...)::RWGSheet

    @testpos(P)
    @testpos(halfnint)
    @testpos(pdiv)
    nrim ≥ 0 || throw(ArgumentError("nrim must be nonnegative"))
    length(patternvec) ≥ halfnint * (halfnint + 1) ÷ 2 || throw(ArgumentError("patternvec is not long enough"))
    all(x -> isone(x) || iszero(x), patternvec) || throw(ArgumentError("patternvec must consist of ones and zeros"))

    s = SymPixels(halfnint, nrim)
    patternmat = sympixmat(s, patternvec)

    sheet = pixels(; P, patternmat, units, pdiv, quad, class, kwargs...)
    return sheet

end

"""
    pixels(; P, patternmat, units, kwargs...)

# Description:

Create a variable of type `RWGSheet` that contains the triangulation for an arbitrarily pixelated
square unit cell. The returned value has fields `s₁`, `s₂`, `β₁`, `β₂`, `ρ`, `e1`,
`e2`, `fv`, `fe`, and `fr` properly initialized.


# Arguments:

All arguments are keyword arguments which can be entered in any order.

## Required arguments:
- `P::Real > 0`: The side length of the square unit cell, specified in units defined by the `units` keyword argument.
- `patternmat::AbstractMatrix{<:Integer}`: A square matrix, consisting solely of 1's and 0's. The matrix entries control
  the metallization pattern in the unit cell, with a 1 or `true` value denoting a metallized pixel, and a 0 or `false` value
  indicating no metallization.  The `(i,j)` entry corresponds to the pixel centered at `(x,y) = ((j-1/2)d, P-(i-1/2)d)`, where
  `d = P / size(patternmat, 1)`.
- `units`:  Length units for `P` (either `mm`, `cm`, `inch`, or `mil`).

## Optional arguments:
- `pdiv::Int = 1`: The number of "chops" or subdivisions applied to each square pixel side when forming the geometry triangulation.
  A value of `1` (the default) means that the pixels are not subdivided any further, except for a single diagonal across each
  pixel to form triangles.  A value of `2` would subdivide each pixel into 4 squares.  A diagonal edge would be added to
  each of the resulting squares to form triangles.
- `quad::Bool=false`:  If `true`, each subpixel (or pixel, if `pdiv` is 1) is divided into four triangles by adding
  two diagonals.  If `false` (the default), only a single diagonal is added to produce two triangles.
- `sym::Bool = false`: A `true` value states that the pixel matrix has vertical and horizontal mirror symmetry, and
  consists of an even number of rows (and columns).  If `true`, the diagonals added to the squares to form
  triangles will exhibit the same left-right and up-down mirror symmetry as the collection of `true` (`false`) pixel
  locations.  If `sym=false` and `quad=false`, then all added diagonals across the unit cell will have the same
  orientation.  `sym` has no effect (i.e. is redundant) if `quad=true` since
  in that case the two added diagonals already ensure mirror symmetry of the triangulation.
$(pix_optional_kwargs)
"""
function pixels(; P::Real, patternmat::AbstractMatrix{<:Integer}, pdiv::Integer=1, class::Char='M',
    sym::Bool=false, quad::Bool=false, units, kwarg...)

    kwargs = Dict{Symbol,Any}(kwarg)
    haskey(kwargs, :fufp) || (kwargs[:fufp] = true)
    kwargs[:class] = class
    check_optional_kw_arguments!(kwargs)
    @testpos(P)
    @testpos(pdiv)
    all(x -> isone(x) || iszero(x), patternmat) || throw(ArgumentError("patternmat must consist of ones and zeros"))
    size(patternmat,1) == size(patternmat,2) || throw(ArgumentError("patternmat must be a square matrix"))

    if kwargs[:class] == 'M'
        patternmat = (!).(patternmat) # Triangulate the non-metal regions
    else
        @warn """

        `class='J'` detected for `sympixels` or `pixels` element.
        This is is known to produce incorrect results for certain geometries containing metallic islands that intersect in a single point.
        It is strongly recommended to use `class='M'` for most `sympixels` and `pixels` elements.
        """
    end

    s1 = [P, 0.0]
    s2 = [0.0, P]

    d = P / size(patternmat, 1) # Pixel side length
    npside = round(Int, P / d) # pixels on side of unit cell

    xrequired = range(start = 0.0, stop = P, length = 1 + pdiv * npside) |> collect
    yrequired = copy(xrequired)

    function is_inside(x::Real, y::Real)
        # predicate to determine if a point is within the region to be triangulated
        (x < 0 || x > P || y < 0 || y > P) && return false
        # column index into patternmat:
        if iszero(x)
            j = 1
        elseif x == P
            j = npside
        else
            j = 1 + trunc(Int, x / d)
        end
        # row index into patternmat:
        if iszero(y)
            i = npside
        elseif y == P
            i = 1
        else
            i = 1 + trunc(Int, (P - y) / d)
        end
        return isone(patternmat[i,j]) ? true : false
    end

    # Triangulate
    npixtri = count(isone, patternmat) # Number of triangulated pixels
    areat = npixtri * d^2  # Total area to triangulate
    ntri = npixtri * pdiv^2 * 2 # Number of triangles expected
    ntri_requested = ntri ÷ 4 # Much smaller to ensure no extra triangles are added by make_plaid_mesh
    if sym && !quad
        xrequired .-= xrequired[(length(xrequired)+1) ÷ 2]
        yrequired .-= yrequired[(length(yrequired)+1) ÷ 2]
        #xrequired = xrequired[((length(xrequired)+1) ÷ 2):end]
        #yrequired = yrequired[((length(xrequired)+1) ÷ 2):end]
        is_inside_sym(x,y) = is_inside(x + P/2, y + P/2)
        sheet = make_sym_plaid_mesh(xrequired, yrequired, areat, ntri_requested, is_inside_sym)
    else
        sheet = make_plaid_mesh(xrequired, yrequired, areat, ntri_requested, is_inside, quad)
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

    sheet.style = "pixels"
    sheet.ξη_check = true
    sheet.units = units
    sheet.s₁ = SV2(s1)
    sheet.s₂ = SV2(s2)
    sheet.β₁, sheet.β₂ = s₁s₂2β₁β₂(sheet.s₁, sheet.s₂)

    return sheet

end # function pixels
