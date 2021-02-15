module Rings
export Ring

"""
    Ring(r::Integer)

An iterator of (m,n) pairs comprising the r'th summation ring.
"""
struct Ring
    r::Int
    Ring(r) = r < 0 ? error("ring must be non-negative") : new(r)
end

function Base.iterate(ring::Ring)
    r = ring.r
    #topbot = ((m,n) for m in (r,-r) for n in -r:r)
    #leftright = ((m,n) for m in 1-r:r-1 for n in (-r,r))
    m = r
    n = -r
    topdone = botdone = leftdone = rightdone = (r == 0)
    item = (m,n)
    state = (m, n, topdone, botdone, leftdone, rightdone)
    return (item, state)
end

function Base.iterate(ring::Ring, state)
    r = ring.r
    (m, n, topdone, botdone, leftdone, rightdone) = state
    rightdone && return nothing
    if leftdone
        # Doing right
        (m,n) == (r-1,-r) && ((m,n) = (-r,r)) # Just finished left
        m == r-2 && (rightdone = true) # last iteration
        m += 1
        return (m,n), (m, n, topdone, botdone, leftdone, rightdone)
    end
    if botdone
        # Doing left
        (m,n) == (-r,r) && (n = -r) # Just finished bottom
        m == r-2 && (leftdone = true)
        m += 1
        return (m,n), (m, n, topdone, botdone, leftdone, rightdone)
    end
    if topdone
        # Doing bottom
        (m,n) == (r,r) && ((m,n) = (-r, -(r+1))) # Just finished top
        n == r-1 && (botdone = true)
        n += 1
        return (m,n), (m, n, topdone, botdone, leftdone, rightdone)
    end
    # Doing top
    n == r-1 && (topdone = true)
    n += 1
    return (m,n), (m, n, topdone, botdone, leftdone, rightdone)
end

Base.IteratorSize(::Type{Ring}) = Base.HasLength()
Base.IteratorEltype(::Type{Ring}) = Base.HasEltype()
Base.eltype(::Type{Ring}) = Tuple{Int,Int}
Base.length(ring::Ring) = ring.r == 0 ? 1 : 8*ring.r


end # module
