#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Meanderline/Strip-Based CPSS
# This example comes from the same authors as the previous example.  The paper is
# A. Ericsson and D. Sjöberg, "Design and Analysis of a Multilayer Meander Line
# Circular Polarization Selective Structure", IEEE Trans. Antennas Propagat., 
# Vol. 65, No. 8, Aug 2017, pp. 4089-4101.
# The design is similar to that of the previous example except that here, the two ``\pm 45^\circ``
# rotated meanderlines are replaced with rectangular strips.  
# This allows us to employ the `diagstrip` element and the `orient` keyword for the 
# `meander` elements to maintain the same, square unit cell for all sheets. By doing this
# we allow PSSFSS to rigorously account for the inter-sheet coupling using multiple 
# high-order modes in the generalized scattering matrix (GSM) formulation. 
#
#-
# We begin by computing the skin depth and sheet resistance for the 
# copper traces.  The conductivity and thickness are as stated in the paper:

## Compute skin depth and sheet resistance:
using PSSFSS.Constants: μ₀ # free-space permeability [H/m]
f = (10:0.1:20) * 1e9 # frequencies in Hz
σ = 58e6 # conductivity of metalization [S/m]
t = 18e-6 # metalization thickness [m]
Δ = sqrt.(2 ./ (2π*f*σ*μ₀)) # skin depth [m]
@show extrema(t./Δ)
#-
Rs = 1 ./ (σ * Δ)
@show extrema(Rs)

# We see that the metal is many skin depths thick (effectively infinitely thick) so that we can use
# the thick metal surface sheet resistance formula.  Since the latter varies with frequency, we approximate
# it over the band 10-20 GHz by a value near its mean: 0.032 Ω/□.

# Here is the script that analyzes the design from the referenced paper:

using PSSFSS
P = 5.2 # side length of unit square
d1 = 2.61 # Inner layer thickness
d2 = 3.81 # Outer layer thickness
h0 = 2.44 # Inner meanderline dimension (using paper's definition of h)
h2 = 2.83 # Outer meanderline dimension (using paper's definition of h)
w0x = 0.46 # Inner meanderline line thickness of traces running along x
w0y = 0.58 # Inner meanderline line thickness of traces running along y
w1 = 0.21 # Rectangular strips width
w2x = 0.25   # Outer meanderline line thickness of traces running along x
w2y = 0.17 # Outer meanderline line thickness of traces running along y

outer(orient) = meander(a=P, b=P, w1=w2y, w2=w2x, h=h2+w2x, units=mm, ntri=600, 
                        Rsheet=0.032, orient=orient)
inner = meander(a=P, b=P, w1=w0y, w2=w0x, h=h0+w0x, units=mm, ntri=600, Rsheet=0.032)
strip(orient) = diagstrip(P=P, w=w1, units=mm, Nl=60, Nw=4, orient=orient, Rsheet=0.032)

substrate = Layer(width=0.127mm, epsr=2.17, tandel=0.0009)
foam(w) = Layer(width=w, epsr=1.043, tandel=0.0017)
sheets = [outer(-90), strip(-45), inner, strip(45), outer(90)]
strata = [
    Layer()
    substrate
    sheets[1]
    foam(d2 * 1mm)
    substrate
    sheets[2]
    foam(d1 * 1mm)
    sheets[3]
    substrate
    foam(d1 * 1mm)
    substrate
    sheets[4]
    foam(d2 * 1mm)
    sheets[5]
    substrate
    Layer() ]
steering = (θ=0, ϕ=0)
flist = 10:0.1:20

results = analyze(strata, flist, steering, logfile=devnull, 
                  resultfile=devnull, showprogress=false)
#md nothing # hide

# The PSSFSS run took about 85 seconds on my machine.  Here are plots of the five sheets:

using Plots
default()
ps = []
for k in 1:5
    push!(ps, plot(sheets[k], unitcell=true, title="Sheet $k", linecolor=:red))
end
plot(ps..., layout=5)
#md savefig("cpssb1.png"); nothing  # hide
#-
#md # ![](cpssb1.png)
#-
# Notice that for all 5 sheets, the unit cell is a square of constant side length and is unrotated. 
# We can see from the log file (of a previous run where it was not suppressed) that this allows
# PSSFSS to use additional modes in the GSM cascading procedure:
#-
# ```
# Starting PSSFSS analysis on 2021-05-26 at 09:54:02.902
# 
# 
# Dielectric layer information... 
# 
#  Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
#  ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
#      1    0.0000  mm    1.00 0.0000    1.00 0.0000     2  1208.3     0.0     0.0  1208.3
#      2    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
#  ==================  Sheet   1  ========================  1208.3     0.0     0.0  1208.3
#      3    3.8100  mm    1.04 0.0017    1.00 0.0000    10  1208.3     0.0     0.0  1208.3
#      4    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
#  ==================  Sheet   2  ========================  1208.3     0.0     0.0  1208.3
#      5    2.6100  mm    1.04 0.0017    1.00 0.0000    18  1208.3     0.0     0.0  1208.3
#  ==================  Sheet   3  ========================  1208.3     0.0     0.0  1208.3
#      6    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      7    2.6100  mm    1.04 0.0017    1.00 0.0000    18  1208.3     0.0     0.0  1208.3
#      8    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
#  ==================  Sheet   4  ========================  1208.3     0.0     0.0  1208.3
#      9    3.8100  mm    1.04 0.0017    1.00 0.0000    10  1208.3     0.0     0.0  1208.3
#  ==================  Sheet   5  ========================  1208.3     0.0     0.0  1208.3
#     10    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
#     11    0.0000  mm    1.00 0.0000    1.00 0.0000     2  1208.3     0.0     0.0  1208.3
# ...
# ```
#-
# Layers 3 and 9 were assigned 10 modes each.  Layers 5 and 7, being thinner were assigned 
# 18 modes each. The numbers of modes are determined automatically by PSSFSS to ensure 
# accurate cascading.  
#-
# Here are comparison plots of PSSFSS versus highly converged CST predictions digitized from 
# plots presented in the paper:

using Plots, DelimitedFiles
RLl = -extract_result(results, @outputs s11db(l,l))
AR11l = extract_result(results, @outputs ar11db(l))
IL21r = -extract_result(results, @outputs s21db(r,r))
AR21r = extract_result(results, @outputs ar21db(r))

default(lw=2, xlabel="Frequency (GHz)", xlim=(10,20), xtick=10:2:20,
        framestyle=:box, gridalpha=0.3)

plot(flist,RLl,title="LHCP → LHCP Return Loss", label="PSSFSS",
         ylabel="Return Loss (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_rllhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
#md savefig("cpssb2.png"); nothing  # hide
#-
#md # ![](cpssb2.png)
#-
plot(flist,AR11l,title="LHCP → LHCP Reflected Axial Ratio", label="PSSFSS",
         ylabel="Axial Ratio (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_arlhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
#md savefig("cpssb3.png"); nothing  # hide
#-
#md # ![](cpssb3.png)
#-
plot(flist,AR21r,title="RHCP → RHCP Transmitted Axial Ratio", label="PSSFSS",
     ylabel="Axial Ratio (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_arrhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
#md savefig("cpssb4.png"); nothing  # hide
#-
#md # ![](cpssb4.png)
#-
plot(flist,IL21r,title="RHCP → RHCP Insertion Loss", label="PSSFSS",
         ylabel="Insertion Loss (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_ilrhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
#md savefig("cpssb5.png"); nothing  # hide
#-
#md # ![](cpssb5.png)
#-

# Differences between the PSSFSS and CST predictions are attributed to the fact that the 
# metalization thickness of 18 μm was included in the CST model but cannot be accommodated by PSSFSS.
