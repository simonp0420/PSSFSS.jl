using Plots, PSSFSS
patch = rectstrip(Nx=10, Ny=10, Px=1, Py=1, Lx=0.5, Ly=0.5, units=cm)
p1 = plot(patch, unitcell=true, linecolor=:red)
p2 = plot(patch, rep=(3,3), linecolor=:blue)
plot(p1, p2, layout=(1,2))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
