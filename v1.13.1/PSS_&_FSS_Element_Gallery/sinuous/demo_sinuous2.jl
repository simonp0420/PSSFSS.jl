using PSSFSS, Plots
s1 = [1, 0]; s2 = [0, 1]
b = [0.12, 0.2, 0.3]
sides = 50; ntri = 2800; units = cm
sheet = sinuous(; arms=4, b, w=0.03, rc=0.05, s1, s2,
                   L2=0.95, w2=0.03, c2=0.12, g=0.04, sides, ntri, units)


p1 = plot(sheet, linecolor=:red, size=(400,400), unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(3,3))
plot(p1, p2, layout = (1,2), size=(800,400), margin=10Plots.pt)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
