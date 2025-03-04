#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Meanderline-Based CPSS
# A "CPSS" is a circular polarization selective structure, i.e., a structure that reacts differently
# to the two senses of circular polarization.  
# We'll first look at analyzing a design presented in the literature, and then proceed to optimize another 
# design using PSSFSS as the analysis engine inside the optimization objective function.
# ### Sjöberg and Ericsson Design
# This example comes from [sjoberg2014multi](@cite).
# The authors describe an ingenious structure consisting of 5 progressively rotated meanderline sheets, which
# acts as a circular polarization selective surface: it passes LHCP (almost) without attenuation or 
# reflection, and reflects RHCP (without changing its sense!) almost without attenuation or transmission.
#

# Here is the script that analyzes their design:

using PSSFSS
## Define convenience functions for sheets:
outer(rot) = meander(a=3.97, b=3.97, w1=0.13, w2=0.13, h=2.53+0.13, units=mm, ntri=600, rot=rot)
inner(rot) = meander(a=3.97*√2, b=3.97/√2, w1=0.1, w2=0.1, h=0.14+0.1, units=mm, ntri=600, rot=rot)
center(rot) = meander(a=3.97, b=3.97, w1=0.34, w2=0.34, h=2.51+0.34, units=mm, ntri=600, rot=rot)
## Note our definition of `h` differs from that of the reference by the width of the strip.
t1 = 4mm # Outer layers thickness
t2 = 2.45mm # Inner layers thickness
substrate = Layer(width=0.1mm, epsr=2.6)
foam(w) = Layer(width=w, epsr=1.05) # Foam layer convenience function
rot0 = 0 # rotation of first sheet
strata = [
    Layer()
    outer(rot0)
    substrate
    foam(t1)
    inner(rot0 - 45)
    substrate
    foam(t2)
    center(rot0 - 2*45)
    substrate
    foam(t2)
    inner(rot0 - 3*45)
    substrate
    foam(t1)
    outer(rot0 - 4*45)
    substrate
    Layer() ]
steering = (θ=0, ϕ=0)
flist = 10:0.1:20
##
results = analyze(strata, flist, steering, showprogress=false, 
                  resultfile=devnull, logfile=devnull); 

# Here are plots of the five meanderline sheets:

using Plots
plot(outer(rot0), unitcell=true, title="Sheet1")
#-
plot(inner(rot0-45), unitcell=true, title="Sheet2")
#-
plot(center(rot0-2*45), unitcell=true, title="Sheet3 (Center)")
#-
plot(inner(rot0-3*45), unitcell=true, title="Sheet4")
#-
plot(outer(rot0-4*45), unitcell=true, title="Sheet5")

# Notice that not only are the meanders rotated, but so too are the unit cell rectangles.
# This is because we used the generic `rot` keyword argument that rotates the entire unit
# cell and its contents. `rot` can be used for any FSS or PSS element type.  As a consequence of
# the different rotations applied to each unit cell, interactions between sheets due to higher-order
# modes cannot be accounted for; only the dominant ``m=n=0`` TE and TM modes are used in cascading
# the individual sheet scattering matrices.  This approximation is adequate for sheets that are 
# sufficiently separated.  We can see from the log file (saved from a previous run where it was not
# disabled) that only 2 modes are used to model the interactions between sheets:
#
# ```
# Starting PSSFSS 1.0.0 analysis on 2022-09-12 at 16:35:20.914
# Julia Version 1.8.1
# Commit afb6c60d69 (2022-09-06 15:09 UTC)
# Platform Info:
#   OS: Windows (x86_64-w64-mingw32)
#   CPU: 8 × Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
#   Threads: 8 on 8 virtual cores
#   BLAS: LBTConfig([ILP64] libopenblas64_.dll)
# 
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 1 and 2
#    Setting #modes in dividing layer 3 to 2
# ******************* Warning ***********************
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 2 and 3
#    Setting #modes in dividing layer 5 to 2
# ******************* Warning ***********************
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 3 and 4
#    Setting #modes in dividing layer 7 to 2
# ******************* Warning ***********************
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 4 and 5
#    Setting #modes in dividing layer 9 to 2
# ******************* Warning ***********************
# 
# 
# Dielectric layer information... 
# 
#  Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
#  ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
#      1    0.0000  mm    1.00 0.0000    1.00 0.0000     2  1582.7    -0.0    -0.0  1582.7
#  ==================  Sheet   1  ========================  1582.7    -0.0    -0.0  1582.7
#      2    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      3    4.0000  mm    1.05 0.0000    1.00 0.0000     2  1582.7    -0.0    -0.0  1582.7
#  ==================  Sheet   2  ========================   791.3  -791.3  1582.7  1582.7
#      4    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      5    2.4500  mm    1.05 0.0000    1.00 0.0000     2   791.3  -791.3  1582.7  1582.7
#  ==================  Sheet   3  ========================     0.0 -1582.7  1582.7     0.0
#      6    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      7    2.4500  mm    1.05 0.0000    1.00 0.0000     2     0.0 -1582.7  1582.7     0.0
#  ==================  Sheet   4  ========================  -791.3  -791.3  1582.7 -1582.7
#      8    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      9    4.0000  mm    1.05 0.0000    1.00 0.0000     2  -791.3  -791.3  1582.7 -1582.7
#  ==================  Sheet   5  ======================== -1582.7     0.0     0.0 -1582.7
#     10    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#     11    0.0000  mm    1.00 0.0000    1.00 0.0000     2 -1582.7     0.0     0.0 -1582.7
# 
# ...
# ```
# 
# Note that PSSFSS prints warnings to the log file where it is forced to set the number of layer
# modes to 2 because of unequal unit cells.  Also, in the dielectric layer list it can be seen
# that these layers are assigned 2 modes each.  The thin layers adjacent to sheets are assigned 0
# modes because these sheets are incorporated into so-called "GSM blocks" or "Gblocks" wherein
# the presence of the thin layer is accounted for using the stratified medium Green's functions.
# Analyzing this multilayer structure at 101 frequencies required about 22 seconds on my machine.


# Here is the script that compares PSSFSS predicted performance with very
# high accuracy predictions from CST and COMSOL that were digitized from figures in the paper.

using Plots, DelimitedFiles
RL11rr = -extract_result(results, @outputs s11db(r,r))
AR11r = extract_result(results, @outputs ar11db(r))
IL21L = -extract_result(results, @outputs s21db(L,L))
AR21L = extract_result(results, @outputs ar21db(L))

RLgoal, ILgoal, ARgoal  = ([0.5, 0.5], [0.5, 0.5], [0.75, 0.75])
foptlimits = [12, 18]
default(lw=2, xtick=10:20, xlabel="Frequency (GHz)", ylabel="Amplitude (dB)", gridalpha=0.3)

p = plot(flist,RL11rr,title="RHCP → RHCP Return Loss", label="PSSFSS", 
         ylim=(0,2), ytick=0:0.25:2)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_rl.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_rl.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, RLgoal, color=:black, lw=4, label="Goal")
#md savefig("cpssa1.png"); nothing  # hide
#-
#md # ![](cpssa1.png)
#-
p = plot(flist,AR11r,title="RHCP → RHCP Reflected Axial Ratio", label="PSSFSS",
         xlim=(10,20), ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_reflected.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_reflected.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ARgoal, color=:black, lw=4, label="Goal")
#md savefig("cpssa2.png"); nothing  # hide
#-
#md # ![](cpssa2.png)
#-          
p = plot(flist,IL21L,title="LHCP → LHCP Insertion Loss", label="PSSFSS",
         xlim=(10,20), ylim=(0,2), ytick=0:0.25:2)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_il.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_il.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ILgoal, color=:black, lw=4, label="Goal")
#md savefig("cpssa3.png"); nothing  # hide
#-
#md # ![](cpssa3.png)
#-
p = plot(flist,AR21L,title="LHCP → LHCP Transmitted Axial Ratio", label="PSSFSS",
         xlim=(10,20), ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_transmitted.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_transmitted.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ARgoal, color=:black, lw=4, label="Goal")
#md savefig("cpssa4.png"); nothing  # hide
#-
#md # ![](cpssa4.png)


# The PSSFSS results generally track well with the high-accuracy solutions, but are less accurate
# especially at the high end of the band, possibly because in PSSFSS metallization thickness is 
# neglected and cascading is performed for this structure using only the two principal Floquet modes.  
# As previosly discussed, this is necessary because the rotated meanderlines are achieved by rotating 
# the entire unit cell, and the unit cell for sheets 2 and 4 are not square.  Since the periodicity of 
# the sheets in the structure varies from sheet to sheet, higher order Floquet modes common to neighboring
# sheets cannot be defined, so we are forced to use only the dominant (0,0) modes which are independent of
# the periodicity.  This limitation is removed in the next example.
# Meanwhile, it is of interest to note that their high-accuracy runs
# required 10 hours for CST and 19 hours for COMSOL on large engineering workstations versus about 22 
# seconds for PSSFSS on my desktop machine.

