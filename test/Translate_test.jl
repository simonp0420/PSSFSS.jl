# Tests that translating geometry does not break the simulation

using PSSFSS
using Test
using Unitful
using PSSFSS.Sheets:translate!

# Test that the simulation runs normally for a variety of translations and mesh densities

function run_sim(ntri, shift)
    period = 20.0 # um
    l1 = 0.8*period
    l2 = 0.1*period

    flist = [5u"THz"]
    steering = (phi=0, theta=0)
    sheet = loadedcross(s1=[period,0.0], s2=[0.0,period], L1=l1, L2=l2, w=l2, units=μm, ntri=ntri, class='M', dx=shift*period, dy=shift*period)
    strata = [Layer(), sheet, Layer()]
    try
        results = analyze(strata, flist, steering, showprogress=false,resultfile=devnull, logfile=devnull)
        return true
    catch e 
        return false
    end
end

@testset "translate errors" for ntri in [100, 200, 400], shift in [0.0, 0.1, 0.25, 0.45, 0.50]
    @test run_sim(ntri, shift)
end

# Compare simulations using a translated mesh vs a mesh using a phase shift for each type of element. 
# Shifts will be small to avoid the bug fixed by changing to a phase shift
period = 20.0
flist = 2u"THz":2u"THz":10u"THz"
dx = 0.1
dy = 0.1
l1 = 0.5*period
l2 = 0.3*period
ntri = 100
steering = (phi=0, theta=0)
sheet_ps = loadedcross(s1=[period,0.0], s2=[0.0,period], L1=l1, L2=l2, w=l2, units=μm, ntri=ntri, class='M', dx=dx*period, dy=dy*period)
sheet_t = loadedcross(s1=[period,0.0], s2=[0.0,period], L1=l1, L2=l2, w=l2, units=μm, ntri=ntri, class='M')
translate!(sheet_t, dx*period, dy*period)

strata_ps = [Layer(), sheet_ps, Layer()]
strata_t = [Layer(), sheet_t, Layer()]

results_ps = analyze(strata_ps, flist, steering, showprogress=false,resultfile=devnull, logfile=devnull)
results_t = analyze(strata_t, flist, steering, showprogress=false,resultfile=devnull, logfile=devnull)

data_ps = extract_result(results_ps, @outputs s21mag(v,v))
data_t = extract_result(results_t, @outputs s21mag(v,v))

tol = 1e-6

@testset "phase shift vs translation" begin
    @test sum([data_ps[i] - data_t[i] for i in eachindex(flist)]) <= tol
end
ntri=400
sheet_list = []
name_list = []
loaded_sheet = loadedcross(s1=[period,0.0], s2=[0.0,period], L1=0.6*period, L2=0.2*period, w=0.05*period, units=μm, ntri=ntri, class='M')
push!(sheet_list, loaded_sheet)
push!(name_list, "Loaded Cross")
jerusalem_sheet = jerusalemcross(P=period, L1=0.8*period, L2=0.1*period, A=0.4*period, B=0.1*period, w=0.2*period, units=μm, ntri=ntri)
push!(sheet_list, jerusalem_sheet)
push!(name_list, "Jerusalem Cross")
diagstrip_sheet = diagstrip(P=period, w=period*0.5, orient=45, units=μm, Nl=10, Nw=10, class='M')
push!(sheet_list, diagstrip_sheet)
push!(name_list, "Diagstrip")
rectstrip_sheet = rectstrip(Lx=period*0.5, Ly=period*0.5, Nx=10, Ny=10, Px=period, Py=period, units=μm)
push!(sheet_list, rectstrip_sheet)
push!(name_list, "Rectstrip")
meander_sheet = meander(a=period, b=period, w1=period*0.1, w2=period*0.1, ntri=ntri, units=μm, h = period*0.5, class='M')
push!(sheet_list, meander_sheet)
push!(name_list, "Meander")
polyring_sheet = polyring(units=μm, s1=[period,0.0], s2=[0.0,period], a=[0.2*period], b=[0.4*period], sides=8, ntri=ntri)
push!(sheet_list, polyring_sheet)
push!(name_list, "Polyring")
splitring_sheet = splitring(s1=[period,0.0], s2=[0.0,period], units=μm, a=[0.2*period], b=[0.4*period], sides=16, ntri=ntri, gapcenter=0, gapwidth=period*0.1)
push!(sheet_list, splitring_sheet)
push!(name_list, "Split ring")
sinuous_sheet = sinuous(s1=[period,0.0], s2=[0.0,period], units=μm, L2=period*0.9, w=0.075*period, b=[period*0.3, period*0.45], rc=period*0.15, arms=4, sides=20, ntri=ntri, g=period*0.1)
push!(sheet_list, sinuous_sheet)
push!(name_list, "Sinuous")


function compare_results(sheet1, dx, dy, test, debug)
    flist1 = 2u"THz":2u"THz":10u"THz"
    tol = 1e-6
    steering1 = (phi=0, theta=0)

    # Set values for dx and dy for phase shift
    sheet_ps = deepcopy(sheet1)
    sheet_ps.dx = dx
    sheet_ps.dy = dy

    # Translate sheet
    sheet_t = deepcopy(sheet1)
    if test
        translate!(sheet_t, dx, dy)
    end

    strata_ps = [Layer(), sheet1, Layer(width=10μm), sheet_ps, Layer()]
    strata_t = [Layer(),sheet1, Layer(width=10μm),  sheet_t, Layer()]

    try
        results_ps = analyze(strata_ps, flist1, steering1, showprogress=false,resultfile=devnull, logfile=devnull)
        results_t = analyze(strata_t, flist1, steering1, showprogress=false,resultfile=devnull, logfile=devnull)
        data_ps = extract_result(results_ps, @outputs s21(v,v))
        data_t = extract_result(results_t, @outputs s21(v,v))
        # return true if the different between the two simulations is less than the desired tolerance
        total = sum([abs((data_ps[i] - data_t[i]))/abs(data_ps[i]) for i in eachindex(flist)]) / size(flist,1)
        #println("\t", total)
        return total <= tol
    catch e 
        println("Simulation Failed")
        return false
    end
end

@testset "translate vs phase shift" begin
    for (i, sheet) in enumerate(sheet_list)
        @testset "$(name_list[i])" begin
            @test compare_results(sheet, (period*rand(Float32)*0.2), (period*rand(Float32)*0.2), true, false)
            @test !compare_results(sheet, (period*rand(Float32)*0.2), (period*rand(Float32)*0.2), false, false)
        end
    end
end