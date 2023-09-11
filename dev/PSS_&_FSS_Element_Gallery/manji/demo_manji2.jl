using PSSFSS, Plots
L1, L2, L3, L4 = 0.481, 0.22, 0.24, 1.1
w = 0.06
a = 0.15
w2 = 0.03
s1=[1.2, 0]; s2=[0, 1.2]
sheet = manji(; s1, s2, units=cm, L1, L2, L3, L4, w, w2, a, CCW=true, ntri=1000)
p1 = plot(sheet, linecolor=:red, size=(400,400), unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(3,3))
plot(p1, p2, layout = (1,2), size=(800,400), margin=10Plots.pt)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
