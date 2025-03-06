#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Checkerboard Metasurface
# Checkerboard-style metasurfaces have been exensively studied, for example in [nakata2013plane](@cite)
# and [kempa2010percolation](@cite).  They exhibit some very unusual properties, which we will demonstrate here.
#
# We consider a series of checkerboad-like geometries.  The square unit cell has dimension ``P = 5~\text{mm}``.
# The side length for a PEC square sheet rotated 45° and perfectly inscribed in the unit cell is ``L_0 = P / \sqrt{2}``.
# We will analyze a series of squares at 1 GHz with side lengths ``L = L_0 + δ``.  The function `computed_and_plot`
# defined below will take a given value of ``δ`` as its input, then perform four distinct analyses based on two 
# complementary triangulations.  Let's run it first for ``δ = -0.5 \text{mm}`` and look
# at the resulting triangulations:
using PSSFSS, Plots, PrettyTables
using Printf: @sprintf

units = mm
P = 5  # period in x and y direction
L0 = P / √2 # Side length for self-complementary square
Nx = Ny = 10

function compute_and_plot(δ; P=P, units=units, L0=L0, Nx=Nx, Ny=Ny)

    Px = Py = P
    Lx = Ly = L0 - abs(δ)
    orient = 45
    islekwds = (; Lx, Ly, Px, Py, units, Nx, Ny, orient)
    holekwds = (; units, s1=[P,0], s2=[0,P], sides=4, 
                  a = iszero(δ) ? [P/2] : [(L0-abs(δ)) / √2], 
                  b=[-2*Nx], ntri=2*Nx*Ny, structuredtri=false)
    if δ ≤ 0
        redj = rectstrip(; class='J', islekwds...)
        redm = rectstrip(; class='M', islekwds...)
        bluej = polyring(; class='J', holekwds...)
        bluem = polyring(; class='M', holekwds...)
    else
        redj = polyring(; class='J', holekwds...)
        redm = polyring(; class='M', holekwds...)
        bluej = rectstrip(; class='J', islekwds...)
        bluem = rectstrip(; class='M', islekwds...)
    end

    p1 = plot(redj, unitcell=true, faces=true, fillcolor=:red, title="Unit Cell") 
    p2 = plot(redj, faces=true, edges=false, fillcolor=:red, rep=(3,3), title="3×3 Array")
    p3 = plot(bluej, unitcell=true, faces=true, fillcolor=:blue, title="Unit Cell") 
    p4 = plot(bluej, faces=true, edges=false, fillcolor=:blue, rep=(3,3), title="3×3 Array")
    pl = plot(p1,p2,p3,p4, layout=(2,2), size=(550,600), plot_title=" Triangulations for δ = $δ")

    flist = 1  # Analysis frequency in GHz
    steering = (θ=0, ϕ=0) # Normal incidence
    logfile = resultfile = devnull  # Suppress creation of output files
    showprogress = false  # Suppress screen output

    redjres = analyze([Layer(), redj, Layer()], flist, steering; logfile, resultfile, showprogress)
    bluejres = analyze([Layer(), bluej, Layer()], flist, steering; logfile, resultfile, showprogress)
    redmres = analyze([Layer(), redm, Layer()], flist, steering; logfile, resultfile, showprogress)
    bluemres = analyze([Layer(), bluem, Layer()], flist, steering; logfile, resultfile, showprogress)


    s11rj = extract_result(redjres, @outputs s11(v,v)) |> only
    s21rj = extract_result(redjres, @outputs s21(v,v)) |> only
    s11rm = extract_result(redmres, @outputs s11(v,v)) |> only
    s21rm = extract_result(redmres, @outputs s21(v,v)) |> only
    s11bj = extract_result(bluejres, @outputs s11(v,v)) |> only
    s21bj = extract_result(bluejres, @outputs s21(v,v)) |> only
    s11bm = extract_result(bluemres, @outputs s11(v,v)) |> only
    s21bm = extract_result(bluemres, @outputs s21(v,v)) |> only

    return pl, (; s11rj, s21rj, s11rm, s21rm, s11bj, s21bj, s11bm, s21bm)
end

pl, _ = compute_and_plot(-0.5)
pl
# The function creates a "red" triangulation occupying the triangle of side length ``L_0 + δ``, and a "blue" triangulation, 
# occupying the complementary portion of the unit cell. Since for this case the red square
# side length is 0.5 mm shorter than the critical length ``L_0``, it lies strictly inside the unit cell.  So if we choose to
# use the red triangulation to model electric surface current, then we can consider the red regions to be "islands" of 
# metal in otherwise empty space.  We could also use the blue triangulation to model magnetic surface current, which again would 
# lead to the conclusion that the small untriangulated squares are conducting patches or "islands" of metalization.  
# Either of these two choices, when analyzed with PSSFSS, should yield the same values for computed reflection or transmission 
# coefficients (within modeling accuracy).  
#
# A different approach would be to choose the red triangulation for representing magnetic surface current, in which case 
# the small red squares would represent "holes" in an otherwise solid metallic sheet. The same "hole" interpretation would 
# follow from using the blue triangulation to represent electric surface current. In fact, for this case, the blue region 
# in the full plane can be regarded as the union of an infinite number of metallic squares of dimension ``L_0 + δ``.  
# So positive values of ``δ`` can be handled by using ``-δ`` and reversing the roles of the red and blue triangulations.
# This is exactly what is done in the function `compute_and_plot` above.
#-
# We'll now exercise the function for the set of ``δ`` values ``\{ -0.2, -0.05, 0, 0.05, 0.2 \}``, observing both the plotted
# triangulations and the resulting scattering parameter predictions for each of the four modeling choices outlined above.
δs = [-0.2, -0.05, 0, 0.05, 0.2] # Departure in mm from self-complementary square side length 
s11rj, s21rj, s11rm, s21rm, s11bj, s21bj, s11bm, s21bm = (zeros(ComplexF64, length(δs)) for _ in 1:8)
pls = [] #hide
for (i, δ) in pairs(δs)
    plt, r = compute_and_plot(δ)
    s11rj[i] = r.s11rj
    s21rj[i] = r.s21rj
    s11rm[i] = r.s11rm
    s21rm[i] = r.s21rm
    s11bj[i] = r.s11bj
    s21bj[i] = r.s21bj
    s11bm[i] = r.s11bm
    s21bm[i] = r.s21bm
    push!(pls, plt) #hide
end
#-
i = 1   #hide
pls[i]   #hide
#-
i += 1   #hide
pls[i]   #hide
#-
i += 1   #hide
pls[i]   #hide
#-
i += 1   #hide
pls[i]   #hide
#-
i += 1   #hide
pls[i]   #hide
#-
# Let the letter "J" denote use of a triangulation to represent electric surface current and let "M" denote magnetic surface current.
# So, e.g. "Blue M" means that the blue triangulation is used to represent magnetic current.  We can make the following observations about
# the above plots:
#
# 1. The red and blue triangulations alway occupy complementary regions of the unit cell. Then for a given choice
#    of ``δ``, PSSFSS analysis of "Blue J" and "Red M" should result identical scattering parameters (apart from modeling errors).
#    Likewise analysis of "Red J" and "Blue M" should result in the same scattering parameters.
# 2. It is well known from Babinet's principle [nakata2013plane](@cite) [tan2012babinet](@cite) that reflection coefficients
#    of thin perforated screens are equal to the negative of the transmission coefficients of the complementary structure, provided the
#    structures exhibit sufficient rotational symmetry as is the case here.  "Red J" and "Blue J" form such a 
#    complementary pair, as do "Red M" and "Blue M".
# 3. All of the ``δ`` values are quite small compared to ``L_0 \approx 3.5355 \text{mm}``.  In particular, the plots for 
#    ``δ = \pm 0.05`` are almost indistinguishable from the plot for ``δ = 0``.  Therefore, we expect the scattering 
#    parameters of all of these cases to be nearly equal.
# 4. For the case of ``δ = 0``, the red and blue regions are both self-complementary. From our earlier considerations, we then
#    expect equal reflection and transmission coefficients for any of the choices "Red J", "Blue J", "Red M", or "Blue M".  In 
#    fact all four of these cases should yield the same values.
#
# The following code generates a summary table showing how well these expectations are satisfied:
mat = hcat(s11rj, s11bm, -s21rm, -s21bj) |> transpose
row_labels = ["S₁₁: Red J ", "S₁₁: Blue M", "-S₂₁: Red M ", "-S₂₁: Blue J"]
header = (["δ = $δ mm" for δ in δs], ["islands or holes" for δ in δs])
header[2][findfirst(iszero, δs)] = "SC (Self-Complementary)"
formatters = (v,i,j) -> imag(v) > 0 ? @sprintf("%7.4f + %6.4fim", real(v), imag(v)) : 
                                      @sprintf("%7.4f - %6.4fim", real(v), -imag(v))
highlighters = Highlighter((data, i, j) -> (j == findfirst(iszero,δs)), crayon"red bold")
pretty_table(mat; header, row_labels, alignment=:c, formatters, highlighters)
#-
# For ``δ < 0`` the results seem reasonable.  In this regime of electrically small metal islands (for Red J and Blue M), 
# one wouldn't expect much reflection, and that is what is observed.  For ``δ > 0``, the geometry (for Red J and Blue M)
# consists of a metal plate with small holes, so one would expect almost total reflection, as observed.  Now consider 
# how well our above observations are aligned with the computed results...
#
# From observations 1 and 2, the numbers in any one column should all be approximately equal, which they are
# except for the center ``δ = 0`` column (highlighted in red).  From observation 3 we expect all the entries 
# in any one row to be nearly equal, but this is not true at all.  In any one row there is a violent jump in 
# amplitude at or near ``δ = 0``.  Finally, from observation 4 we expect all of the numeric entries for ``δ = 0``
# to be approximately equal--which they are most definitely not.  What is going on here?
#
# As discussed in the previously cited references, the idealized model being analyzed for ``δ = 0`` is unphysical. 
# As ``\delta \rightarrow 0`` the response of the structure does not approach a limit, since the violent jump in 
# reflection coefficient occurs for arbitrarily small positive and negative excursions from ``L_0``.  Any physically
# realizable structure must exhibit continuous dependence on its physical parameters. If one attempts
# to approximate this ideal surface in the real world one must use finite thickness and conductivity, so that neighboring
# unit cells will intersect all along the thickness of the metal, rather than at a single point, thus destroying the 
# self-complementary property.  Similarly, finite losses in any real metal preclude self-complementarity.
#
# The same problem observed here for an ideal model of a self-coplementary surface arises when analyzing pixelated 
# structures such as shown in the next example, where adjacent metal pixels 
# can intersect at only a single point. In that case, using a "Blue M" approach is known to agree with measurements.