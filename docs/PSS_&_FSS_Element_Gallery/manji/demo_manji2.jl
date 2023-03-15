# ---
# title: Manji (Counter-clockwise)
# cover: "assets/demo_manji2.png"
# description: "Manji (counter-clockwise) in a square lattice, with folded arms, center square, and no outer ring, created by the manji function"
# ---

using PSSFSS, Plots
L1 = 0.481
L2 = 0.22
w = 0.06
L3 = 0.24 
a = 0.15
w2 = 0.0
sheet = manji(; s1=[1.1, 0], s2=[0, 1.1], units=cm, L1, L2, L3, w, w2, a, CCW=true, ntri=1000)
p1 = plot(sheet, linecolor=:red, size=(400,400), unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(3,3))
plot(sheet, axis=false, xlabel="", ylabel="", grid=false, linecolor=:green, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_manji2.png") #src
plot(p1, p2, layout = (1,2), size=(800,400), margin=10Plots.pt)
