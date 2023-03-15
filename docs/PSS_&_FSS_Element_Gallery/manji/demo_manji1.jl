# ---
# title: Manji (Clockwise)
# cover: "assets/demo_manji1.png"
# description: "Manji (clockwise) in a square lattice, with folded arms, center square, and outer ring, created by the manji function"
# ---

using PSSFSS, Plots
L1 = 0.481
L2 = 0.22
w = 0.06
L3 = 0.24 
a = 0.15
w2 = 0.04
sheet = manji(; s1=[1.1, 0], s2=[0, 1.1], units=cm, L1, L2, L3, w, w2, a, ntri=1000)
p1 = plot(sheet, linecolor=:red, size=(400,400), unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(3,3))
plot(sheet, axis=false, xlabel="", ylabel="", grid=false, linecolor=:green, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_manji1.png") #src
plot(p1, p2, layout = (1,2), size=(800,400), margin=10Plots.pt)
