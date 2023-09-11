using PSSFSS, Plots

a = [4.7, 6.8]
b = [5.9, 8]
s1 = [19.5, 0]
s2 = [0, 19.5]
gapangle = [26, 85.6]
gapcenter = [90, 0]
sheet = splitring(; units=mm, sides=42, ntri=800, a, b, s1, s2, gapangle, gapcenter)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
