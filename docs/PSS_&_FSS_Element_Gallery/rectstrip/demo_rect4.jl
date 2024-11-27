# ---
# title: Rotated Rectangular Strip
# cover: "assets/demo_rect4.png"
# description: "Rotated rectangular strips created by the rectstrip function"
# ---

# This example uses the `rot` keyword and should be compared to the  
# [`diagstrip` example](../diagstrip/demo_diagstrip1.html) and the 
# [`rectstrip` example](./demo_rect5.html) that uses the `orient` keyword.


using Plots, PSSFSS
strip = rectstrip(Nx=10, Ny=2, Px=1, Py=1, Lx=0.5, Ly=0.1, rot=45, units=cm)
p1 = plot(strip, unitcell=true, linecolor=:red)
plot(strip, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:green, size=(300,300), rep=(4,3)) #src
savefig("assets/demo_rect4.png") #src
p2 = plot(strip, rep=(3,3), linecolor=:blue)
plot(p1, p2, layout=(1,2), size=(600,300))
