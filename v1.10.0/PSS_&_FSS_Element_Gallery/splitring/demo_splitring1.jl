using PSSFSS, Plots

units = mm
a = [1.9, 5.7]
b = a + [1.4, 1.3]
gapwidth = [(1.5, 1.5), (1.5, 1.5)]
gapcenter = [(45,225), (45,225)]
s1 = [15, 0]
s2 = [0, 15]
sheet = splitring(; units=mm, sides=45, ntri=700, a, b, s1, s2, gapwidth, gapcenter)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
