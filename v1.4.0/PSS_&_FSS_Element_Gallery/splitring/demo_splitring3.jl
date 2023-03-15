using PSSFSS, Plots

units = mm
s1 = [4, 0]
s2 = [0, 4]
gapwidth = [0.6]
gapcenter = [0]
sides = 4
a = [1.626]
b = [2.475]
sheet = splitring(; units, sides, orient = 45, a, b, ntri=400, s1, s2, gapwidth, gapcenter)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

