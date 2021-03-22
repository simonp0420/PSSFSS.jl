module Meshsub

export meshsub

using Printf: @sprintf
using Triangulate: TriangulateIO, triangulate
using ..Sheets: RWGSheet, MV2


"""
    meshsub(;points, seglist, area, 
            segmarkers=Int[], holes=Array{Float64}(undef,2,0),  switches="") -> sheet::RWGSheet

Triangulate polygonal region(s) and create a sheet object with the triangulation fields properly initialized.

#  Required Input Keyword Arguments:

- `points::Matrix{<:Real}` with size(points,1) == 2. Input list of boundary segment vertex coordinates.
- `seglist::Matrix{<:Integer}` with size(seglist,1) == 2. List of boundary segments for the 
               regions(s) to be triangulated.  seglist[:,n] contains the point 
               indices for the n'th boundary segment. 
- `area::Real` Desired maximum area of any triangle in the generated mesh.
- `ntri::Int`  Desired number of triangles in the generated mesh.

#  Optional Input Keyword Arguments:

- `segmarkers::Vector{<:Integer}` with length(segmarkers) == size(seglist,2). Each boundary 
               polygon is numbered starting with 1. These boundary numbers are entered 
               in segmarkers to denote which boundary a particular segment belongs to.
               Required if more than a single boundary is specified (i.e. for multiply connected regions).
- `holes::Matrix{<:Real} with size(holes,1) == 2.  Each column `n` contains the coordinates
               of a point within the n'th *hole* (i.e., region which is not to be triangulated). There
               must be a boundary polygon defined for the boundary of each hole.  Required if more
               than a single boundary is specified (i.e., for multiply connected regions).
- `switches::String` A set of switches to pass to `Triangulation.triangulate` instead of the default
               `"pDa\$(area)q20.0Qe"`.  If you pass your own switches, be sure to include `"e"`.

#  Return value:
    sheet      A variable of type RWGSheet with fields ρ, e1, e2, fe, and fv properly initialized.

"""    
function meshsub(;points::Matrix{<:Real}, seglist::Matrix{<:Integer},
                 segmarkers=Int[],
                 holes=Array{Float64}(undef,2,0),
                 area::Real=0.0,
                 ntri::Int,
                 switches::String="")::RWGSheet

    size(points,1) == 2 || throw(ArgumentError("points first dimension must be 2"))
    size(seglist,1) == 2 || throw(ArgumentError("seglist first dimension must be 2"))
    isempty(segmarkers) || (length(segmarkers) == size(seglist,2)) ||
        throw(ArgumentError("wrong length(segmarkers)"))
    isempty(holes) || (size(holes,1) == 2) ||
        throw(ArgumentError("holes first dimension must be 2"))
    minangstr = "30.0" # Angle quality requirement
    if isempty(switches)
        iter = 14 # Number of triangulation iterations to reach desired ntri
        astr = @sprintf("%.14f", area)
        switches = "pDa$(astr)q$(minangstr)Qe"
    else
        iter = 1
    end

    # Set up for call to triangulate:
    triin = TriangulateIO()
    triin.pointlist=Matrix{Cdouble}(points)
    triin.segmentlist=Matrix{Cint}(seglist)
    isempty(segmarkers) || (triin.segmentmarkerlist = Vector{Cint}(segmarkers))
    isempty(holes) || (triin.holelist = Matrix{Cdouble}(holes))

    areaold = area
    triout = deepcopy(triin) # Establish scope
    for k in 1:iter
        # Perform triangulation:
        (triout, vorout)=triangulate(switches, triin)
        ntrinew = size(triout.trianglelist,2)
        (0.9ntri ≤ ntrinew ≤ 1.1ntri) && break
        correction = (ntrinew/ntri)^(1/3)
        areanew = areaold * correction
        astr = @sprintf("%.14f", areanew)
        switches = "pDa$(astr)q$(minangstr)Qe"
        areaold = areanew
    end

    # Set up the sheet data structures:
    sh = RWGSheet()
    plist = triout.pointlist
    sh.ρ = [MV2(plist[:,i]) for i in 1:size(plist, 2)]
    sh.e1 = triout.edgelist[1,:]
    sh.e2 = triout.edgelist[2,:]
    sh.fv = copy(triout.trianglelist)
    # Set up sh.fe so that for any triangle, edge 1 connects nodes 1 and 3,
    # edge 2 connects nodes 1 and 2, and edge 3 connects nodes 2 and 3:
    sh.fe = similar(sh.fv) # Preallocation
    for tri in 1:size(sh.fv,2)
        n1, n2, n3 = sh.fv[:,tri]
        for (side, m1, m2) in [(1,n2,n3), (2,n1,n3), (3,n1,n2)]
                found = false
            for i in 1:length(sh.e1)
                if (sh.e1[i] == m1 && sh.e2[i] == m2) || (sh.e1[i] == m2 && sh.e2[i] == m1)
                    sh.fe[side,tri] = i
                    found = true
                    break
                end
            end
            !found && error("triangle edge not found connecting vertices $m1 and $m2 for face $tri")
        end
    end
    sh.fr = zeros(Float64, size(sh.fv,2)) # Face resistance vector
    return sh               
end # function

end # module
