#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Meanderline-Based CPSS
# A "CPSS" is a circular polarization selective structure, i.e., a structure that reacts differently
# to the two senses of circular polarization.  
# We'll first look at analyzing a design presented in the literature, and then proceed to optimize another 
# design using PSSFSS as the analysis engine inside the optimization objective function.
# ### Sjöberg and Ericsson Design
# This example comes from the paper D. Sjöberg and A. Ericsson, "A multi layer meander line circular 
# polarization selective structure (MLML-CPSS)," The 8th European Conference on Antennas and Propagation 
# (EuCAP 2014), 2014, pp. 464-468, doi: 10.1109/EuCAP.2014.6901792.
#
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
# Note our definition of `h` differs from that of the reference by the width of the strip.

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

results = analyze(strata, flist, steering, showprogress=false); 

# Here is the script that compares PSSFSS predicted performance with very
# high accuracy predictions from CST and COMSOL that were digitized from figures in the paper.

using Plots, DelimitedFiles
RL11rr = -extract_result(results, @outputs s11db(r,r))
AR11r = extract_result(results, @outputs ar11db(r))
IL21L = -extract_result(results, @outputs s21db(L,L))
AR21L = extract_result(results, @outputs ar21db(L))

default(lw=2, xlim=(10,20), xtick=10:20, ylim=(0,3), ytick=0:0.5:3, gridalpha=0.3)
p = plot(flist,RL11rr,title="RHCP → RHCP Return Loss", label="PSSFSS")
cst = readdlm("../src/assets/cpss_cst_fine_digitized_rl.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_rl.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
#-
p = plot(flist,AR11r,title="RHCP → RHCP Reflected Axial Ratio", label="PSSFSS")
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_reflected.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_reflected.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
#-          
p = plot(flist,IL21L,title="LHCP → LHCP Insertion Loss", label="PSSFSS")
cst = readdlm("../src/assets/cpss_cst_fine_digitized_il.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_il.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
#-
p = plot(flist,AR21L,title="LHCP → LHCP Transmitted Axial Ratio", label="PSSFSS")
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_transmitted.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_transmitted.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")


# The PSSFSS results generally track well with the high-accuracy solutions, but are less accurate
# especially at the high end of the band, presumably because cascading is performed in PSSFSS
# for this structure using only the two principal Floquet modes.  This is necessary because the 
# rotated meanderlines are achieved by rotating the entire unit cell, and the unit cell for sheets
# 2 and 4 are not square.  Since the periodicity of the sheets in the structure varies from sheet
# to sheet, higher order Floquet modes common to neighboring sheets cannot be defined, so we are forced
# to use only the dominant (0,0) modes which are independent of the periodicity.  This is a limitation
# which could be removed in the future using the same technique employed in the paper to enable full-wave
# analysis with the commercial tools.  Meanwhile, it is of interest to note that their high-accuracy runs
# required 10 hours for CST and 19 hours for COMSOL on large engineering workstations.  The PSSFSS run 
# took about 60 seconds on my desktop machine.


# ### Design Based on PSSFSS Optimization with CMAES
# Here we use PSSFSS in conjunction with the CMAES optimizer from the 
# [CMAEvolutionStrategy](https://github.com/jbrea/CMAEvolutionStrategy.jl) package.  I've used CMAES
# in the past with good success on some tough optimization problems.  Here is the code that defines 
# the objective function:

# ```julia
# using PSSFSS
# using Dates: now
# 
# let bestf = typemax(Float64)
#     function objective(x)
#         period, wo, ho, wi, hi, wc, hc, t1, t2 = x
#         ao = bo = ai = bi = ac = bc = period
#         ai *= √2
#         bi /= √2
#         # Ensure physically realizable or return large value:
#         (bo > ho > 2.1*wo && bi > hi > 2.1*wi && bc > hc > 2.1*wc) || (return 5000.0)
# 
#         outer(rot) = meander(a=ao, b=bo, w1=wo, w2=wo, h=ho, units=mm, ntri=400, rot=rot)
#         inner(rot) = meander(a=ai, b=bi, w1=wi, w2=wi, h=hi, units=mm, ntri=400, rot=rot)
#         center(rot) = meander(a=ac, b=bc, w1=wc, w2=wc, h=hc, units=mm, ntri=400, rot=rot)
# 
#         substrate = Layer(width=0.1mm, epsr=2.6)
#         foam(w) = Layer(width=w, epsr=1.05)
#         rot0 = 0
# 
#         strata = [
#                 Layer()
#                 outer(rot0)
#                 substrate
#                 foam(t1*1mm)
#                 inner(rot0 - 45)
#                 substrate
#                 foam(t2*1mm)
#                 center(rot0 - 2*45)
#                 substrate
#                 foam(t2*1mm)
#                 inner(rot0 - 3*45)
#                 substrate
#                 foam(t1*1mm)
#                 outer(rot0 - 4*45)
#                 substrate
#                 Layer() ]
#         steering = (θ=0, ϕ=0)
#         flist = 11:0.25:19
#         results = analyze(strata, flist, steering, showprogress=false)
#         s11rr, s21ll, ar11db, ar21db = eachcol(extract_result(results, 
#                        @outputs s11db(R,R) s21db(L,L) ar11db(R) ar21db(L)))
#         RL = -s11rr
#         IL = -s21ll
#         obj = maximum(vcat(RL,IL,ar11db,ar21db))
#         if obj < bestf
#             bestf = obj
#             open("optimization_best.log", "a") do fid
#                 xround = round.(x, digits=4)
#                 println(fid, round(obj,digits=4), " at x = ", xround, "  #", now())
#             end
#         end
#         return obj
#     end
# end
# ```

# We optimize at 33 frequencies between 11 and 19 GHz.  The actual frequency range of interest is
# 12 to 18 GHz; the wider optimization band provides some safety margin.  Each objective function
# evaluation takes about 24  seconds on my machine.
# As you can see from the code above, each successive sheet in the structure is rotated an additional
# 45 degrees relative to its predecessor.
# The objective is defined as the largest value of RHCP reflected return loss, LHCP insertion loss, or 
# reflected or transmitted axial ratio that
# occurs at any of the analysis frequencies (i.e. we are setting up for "minimax" optimization). Also,
# the `let` block allows the objective function to maintain persistent state in the
# variable `bestf` which is initialized to the largest 64-bit floating point value. Each time a set 
# of inputs results in a lower objective function value, `bestf` is updated with this value and 
# the inputs and objective function value are
# written to the file "optimization_best.log", along with a date/time stamp.  This allows the user
# to monitor the optimization and to terminate the it prematurely, if desired, without losing the 
# best result achieved so far.

# Here is the code for running the optimization:

# ```julia
# using CMAEvolutionStrategy
# #  x = [period,  wo,  ho,  wi,   hi,  wc,   hc,  t1,  t2]
# xmin = [3.0,     0.1, 0.1, 0.1,  0.1, 0.1,  0.1, 1.5, 1.5]
# xmax = [5.5,     0.35,4.0, 0.35, 4.0, 0.35, 4.0, 6.0, 6.0]
# x0 = 0.5 * (xmin + xmax)
# popsize = 2*(4 + floor(Int, 3*log(length(x0))))
# isfile("optimization_best.log") && rm("optimization_best.log")
# opt = minimize(objective, x0, 1.0;
#            lower = xmin,
#            upper = xmax,
#            maxfevals = 9000,
#            xtol = 1e-4,
#            ftol = 1e-5,
#            popsize=popsize)
# ```

# Note that I set the population size to twice the normal default value.  Based
# on previous experience, using 2 to 3 times the default population size helps the 
# optimizer better on tough objective functions like the present one.
# I let the optimizer run for 6 hours, during which time it reduced the objective function
# value from 11.88 dB to 0.86 dB.  It was then interrupted due to a file system error.  I
# restarted it after setting the starting value to the current best and reducing the 
# "sigma" value (the third argument to `minimize`, which controls the algorithms
#  exploratory tendency) to 0.2 (slightly greater than the value it had achieved during the aborted run)
# from its default value of 1.  After about 13 hours the program terminated normally due to 
# insufficient changes in the `x` variable.  The final value of objective function was 0.65 dB.

# Here is a look at the final portion of the file "optimization_best.log":
# ```
# 0.6535 at x = [3.0968, 0.1025, 2.1601, 0.1003, 0.95, 0.3377, 2.3584, 4.3813, 2.2974]  #2021-05-09T23:02:40.230
# 0.6533 at x = [3.0991, 0.1028, 2.162, 0.1002, 0.944, 0.3372, 2.3562, 4.3801, 2.3007]  #2021-05-09T23:27:38.132
# 0.6532 at x = [3.0985, 0.1029, 2.1652, 0.1003, 0.9414, 0.337, 2.3547, 4.3766, 2.3036]  #2021-05-09T23:46:19.068
# 0.6531 at x = [3.0998, 0.1028, 2.164, 0.1001, 0.9443, 0.3378, 2.3558, 4.3783, 2.3028]  #2021-05-09T23:49:13.202
# 0.6529 at x = [3.0988, 0.1028, 2.1652, 0.1002, 0.9422, 0.3372, 2.3545, 4.3765, 2.3039]  #2021-05-09T23:50:58.064
# 0.6529 at x = [3.0985, 0.1029, 2.1687, 0.1002, 0.9407, 0.3368, 2.3526, 4.3714, 2.3073]  #2021-05-10T00:11:36.314
# 0.6529 at x = [3.0976, 0.1028, 2.1688, 0.1004, 0.9443, 0.3374, 2.3532, 4.3716, 2.3062]  #2021-05-10T00:16:16.205
# 0.6528 at x = [3.0991, 0.1028, 2.1679, 0.1002, 0.9444, 0.3376, 2.3539, 4.3735, 2.3056]  #2021-05-10T00:20:04.181
# 0.6527 at x = [3.0987, 0.1027, 2.1664, 0.1001, 0.9443, 0.3375, 2.3544, 4.3751, 2.3045]  #2021-05-10T00:33:08.779
# 0.6527 at x = [3.0972, 0.1029, 2.1692, 0.1002, 0.9362, 0.336, 2.3512, 4.3707, 2.3083]  #2021-05-10T00:39:14.306
# 0.6527 at x = [3.0985, 0.1028, 2.1663, 0.1001, 0.9426, 0.3372, 2.3539, 4.3749, 2.3049]  #2021-05-10T00:40:59.223
# 0.6527 at x = [3.0974, 0.1028, 2.1688, 0.1002, 0.9424, 0.3371, 2.3527, 4.3718, 2.3069]  #2021-05-10T00:51:08.383
# 0.6526 at x = [3.0982, 0.1028, 2.1678, 0.1001, 0.9422, 0.3371, 2.3531, 4.373, 2.3061]  #2021-05-10T00:56:58.562
# 0.6526 at x = [3.0979, 0.1028, 2.1682, 0.1001, 0.9406, 0.3369, 2.3529, 4.3725, 2.3066]  #2021-05-10T00:59:51.875
# 0.6525 at x = [3.0974, 0.1028, 2.1675, 0.1001, 0.9411, 0.3369, 2.353, 4.3735, 2.306]  #2021-05-10T01:03:39.402
# ```

# The performance of this design is shown below:

# ![](../src/assets/cpss_cmaesopt_ar_refl.png)

# ![](../src/assets/cpss_cmaesopt_ar_trans.png)

# ![](../src/assets/cpss_cmaesopt_il_trans.png)

# ![](../src/assets/cpss_cmaesopt_rl_refl.png)


# It would probably be possible to do improve the performance somewhat over the 12-18 GHz band by tapering the 
# requirements of the objective function at the band edges.  Also, in a serious design effort, several additional
# runs of the optimizer should be attempted, since results vary for stochastic algorithms like CMAES.
