using PSSFSS, Plots
P = 0.55
L2 = 0.95P
s1 = P * [1, 0]; s2 = P * [0, 1]; sides = 45; units = cm
ntri = 1400; w = w2 = 0.03; g = 0.02; rc = 0.05
b=[0.1, 0.15, 0.2]
sheet = sinuous(; arms=2, b, w, rc, g, w2, L2, sides, ntri, units, s1, s2)
p1 = plot(sheet, linecolor=:red, size=(400,400), unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(3,3))
plot(p1, p2, layout = (1,2), size=(800,400), margin=10Plots.pt)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
