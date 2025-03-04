#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Split Ring-Based CPSS
# This circular polarization selective surface (CPSS) example comes from the paper
# [wu2022circular](@cite) by Wu et al.
# The design consists of three sequentially rotated split rings separated by dielectric
# layers. Since the unit cells for all three rings are identical, PSSFSS can rigorously
# account for multiple scattering between the individual sheets using multiple 
# high-order modes in the generalized scattering matrix (GSM) formulation. 
#
#-
# We begin by defining the three `splitring` sheets:
using PSSFSS 
b = [3.8, 4.18, 3.8] # outer radius of each ring
a = b - [1, 1.1, 1]  # inner radius of each ring
gw = [3.1, 1.0, 3.1] # gap widths
gc = [90, 45, 0]     # gap centers
s1 = [10, 0]; s2 = [0, 10] # lattice vectors
units = mm; sides = 42; ntri = 900
sheets = [splitring(;units, sides, ntri, a=[a[i]], b=[b[i]], 
          s1, s2, gapwidth=gw[i], gapcenter=gc[i])  for i in 1:3]

# We generate a plot of the three sheets:
using Plots
default() #hide
ps = []
for i in 1:3
    push!(ps, plot(sheets[i], unitcell=true, lc=:red, title="Sheet $i", size=(400,400)))
end
p = plot(ps..., layout = (1,3), size=(900,300), margin=5Plots.mm)
#md savefig("wu2022_sheets.png"); nothing  #hide 
#-
#md # ![](wu2022_sheets.png)
#-
# Next we define the dielectric layers: `F4B` and `prepreg` (the latter is the bonding agent),
# then set up and run the PSSFSS analysis:
F4B = Layer(ϵᵣ=2.55, tanδ=0.002, width=3mm)
prepreg = Layer(ϵᵣ=3.71, width=0.07mm)
strata = [
    Layer()
    sheets[1]
    F4B
    prepreg
    sheets[2]
    F4B
    prepreg
    sheets[3]
    Layer()
    ]
freqs = 8:0.05:12
steering = (θ = 0, ϕ = 0)
results = analyze(strata, freqs, steering, logfile=devnull, resultfile=devnull, showprogress=false)
#md nothing #hide 
#-
# PSSFSS analysis of this 3-sheet structure at 81 frequencies required 28 seconds on my machine.
# As seen from the portion of the log file below (from a previous run where the log file was not discarded),
# PSSFSS chose 42 modes in layers 2 and 4 to ensure acccurate cascading of the GSMs.
# ```
# Starting PSSFSS 1.2.1 analysis on 2022-12-01 at 09:48:54.807
# Julia Version 1.8.3
# Commit 0434deb161 (2022-11-14 20:14 UTC)
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
# Dielectric layer information... 
# 
#  Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
#  ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
#      1    0.0000  mm    1.00 0.0000    1.00 0.0000     2   628.3    -0.0    -0.0   628.3
#  ==================  Sheet   1  ========================   628.3    -0.0    -0.0   628.3
#      2    3.0000  mm    2.55 0.0020    1.00 0.0000    42   628.3    -0.0    -0.0   628.3
#      3    0.0700  mm    3.71 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#  ==================  Sheet   2  ========================   628.3    -0.0    -0.0   628.3
#      4    3.0000  mm    2.55 0.0020    1.00 0.0000    42   628.3    -0.0    -0.0   628.3
#      5    0.0700  mm    3.71 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#  ==================  Sheet   3  ========================   628.3    -0.0    -0.0   628.3
#      6    0.0000  mm    1.00 0.0000    1.00 0.0000     2   628.3    -0.0    -0.0   628.3
# ```
#-
# The circular polarization reflection and transmission amplitudes are now extracted from the PSSFSS
# results and are plotted along with digitized results from the reference.  We first plot the case
# where the excitation is a LHCP polarized plane wave traveling in the positive $z$ direction, incident upon Region 1:
using DelimitedFiles
default(lw=2)
(s11ll,s11rl,s21ll, s21rl) = eachcol(extract_result(results, @outputs s11db(L,L) s11db(R,L) s21db(L,L) s21db(R,L)))
p = plot(xlabel="Frequency (GHz)", ylabel="Amplitude (dB)", xminorticks=2, yminorticks=2, framestyle=:box,
    xtick=8:12, xlim=(8, 12), ytick = -30:5:0, ylim=(-20,0), legend=:top, gridalpha=0.3)
plot!(p, freqs, s11ll, lc=:black, label = "PSSFSS S11(L,L)")
plot!(p, freqs, s11rl, lc=:red, label = "PSSFSS S11(R,L)")
plot!(p, freqs, s21rl, lc=:blue, label = "PSSFSS S21(R,L)")
plot!(p, freqs, s21ll, lc=:green, label = "PSSFSS S21(L,L)")
data = readdlm("../src/assets/rll_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:black, ls=:dash, label = "Wu S11(L,L)")
data = readdlm("../src/assets/rrl_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:red, ls=:dash, label = "Wu S11(R,L)")
data = readdlm("../src/assets/trl_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:blue, ls=:dash, label = "Wu S21(R,L)")
data = readdlm("../src/assets/tll_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:green, ls=:dash, label = "Wu S21(L,L)")
#md savefig("wu2022_fig2a_compare.png"); nothing  #hide 
#-
#md # ![](wu2022_fig2a_compare.png)
#-
# And then the case where the excitation is a RHCP polarized plane wave:
(s11lr,s11rr,s21lr, s21rr) = eachcol(extract_result(results, @outputs s11db(L,R) s11db(R,R) s21db(L,R) s21db(R,R)))
p = plot(xlabel="Frequency (GHz)", ylabel="Amplitude (dB)", xminorticks=2, yminorticks=2, framestyle=:box,
    xtick=8:12, xlim=(8, 12), ytick = -30:5:0, ylim=(-20,0), legend=:top, gridalpha=0.3)
plot!(p, freqs, s11lr, lc=:black, label = "PSSFSS S11(L,R)")
plot!(p, freqs, s11rr, lc=:red, label = "PSSFSS S11(R,R)")
plot!(p, freqs, s21rr, lc=:blue, label = "PSSFSS S21(R,R)")
plot!(p, freqs, s21lr, lc=:green, label = "PSSFSS S21(L,R)")
data = readdlm("../src/assets/rlr_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:black, ls=:dash, label = "Wu S11(L,R)")
data = readdlm("../src/assets/rrr_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:red, ls=:dash, label = "Wu S11(R,R)")
data = readdlm("../src/assets/trr_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:blue, ls=:dash, label = "Wu S21(R,R)")
data = readdlm("../src/assets/tlr_wu_digitized.csv", ',')
plot!(p, data[:,1], data[:,2], lc=:green, ls=:dash, label = "Wu S21(L,R)")
#md savefig("wu2022_fig2b_compare.png"); nothing  #hide 
#-
#md # ![](wu2022_fig2b_compare.png)
#-
# The agreement between Wu et al and PSSFSS is generally quite good, with larger differences at smaller 
# amplitudes.  This is attributed to the fact that conductor thickness was included in the reference but
# can not yet be accommodated by PSSFSS.
