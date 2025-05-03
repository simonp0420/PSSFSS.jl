using PSSFSS, Plots
L1 = 0.481
L2 = 0.22
w = 0.06
L3 = 0.24
P = 1.15
orient = 30
sheet = manji(; s1=[P, 0], s2=[P/2, P*√3/2], units=cm, L1, L2, L3, w, orient, ntri=1000)
p1 = plot(sheet, linecolor=:red, size=(400,400), unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(3,3))
plot(p1, p2, layout = (1,2), size=(800,400), margin=10Plots.pt)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
