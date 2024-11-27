# ---
# title: Oriented Rectangular Strip
# cover: "assets/demo_rect5.png"
# description: "Oriented rectangular strips created by the rectstrip function"
# ---

# This example uses the `orient` keyword and should be compared to the  
# [`diagstrip` example](../diagstrip/demo_diagstrip1.html) and the 
# [`rectstrip` example](./demo_rect4.html) that uses the `rot` keyword.


using Plots, PSSFSS
strip = rectstrip(Nx=10, Ny=2, Px=1, Py=1, Lx=0.5, Ly=0.1, orient=45, units=cm)
p1 = plot(strip, unitcell=true, linecolor=:red)
plot(strip, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:green, size=(400,300), rep=(4,3)) #src
savefig("assets/demo_rect5.png") #src
p2 = plot(strip, rep=(3,3), linecolor=:blue)
plot(p1, p2, layout=(1,2), size=(600,300))
