import StatsBase as stb
include("Particles.jl")
include("Setup.jl")
include("Stat.jl")

STEP = UInt(0)

#Cell unit
# cu = 20.0
# cu½ = cu / 2
v₀ = 0.95
N = 125
n = 3000
nstart = 0

T = 1.5
v₀ = sqrt(3T)
ρ = 0.3 # N / cu^3
begin
    cu = cbrt(N / ρ)
    cu½ = cu / 2
end
σ_int = pi * 2.5^2
l_fr = 1 / (ρ * σ_int)

# vs = genv(100000; T = T)
# Vssplit, Vscount = freqinintervals(getindex.(vs, 2), 300)
# bar(Vssplit[1:300], Vscount[1:300], bar_width = Vssplit[2] - Vssplit[1])

# energy = sen()
# states, energy = sf()
# virial = svi()
# n = 2
# nstart = 0
# @time simulateall();
states, kin, virial = simulateall(0);
"""
begin
    data_type = Tuple{@Tnv,@Tr,@Tr}
    data = Vector{data_type}(undef, 500)
    states = Vector{@Tnv}(undef, length(data))
    kin = Vector{@Tr}(undef, length(data))
    virial = Vector{@Tr}(undef, length(data))
end
"""
r = iterations!(iterverle!, (r[1], r[2]), 500, process, data, forces, forcessum)
# r_prev, r = genwithprev(gencryst, N, cu; v₀ = v₀);
# energy = cfs(countenergy, states)
# plotstates([r_prev, r], frames = 2, duration = 1, new = "test new")
# # countv(r₁::@Tnv, r₂::@Tnv)::@Tnv() = countv.(r₁, r₂)

pot = potenergy.(states);
energy = pot + kin;

plot(energy[1:end], label = "Kinetic + potential energy", yguide = "E(tᵢ)", xguide = "tᵢ")
plot(pot, label = "Potential energy", widen = true)
plot(kin[1:end], label = "Kinetic energy", widen = true)

average(Vs[end])

begin
    plot(energy[1:end], label = "Kinetic + potential energy", yguide = "E(tᵢ)", xguide = "tᵢ")
    plot!(pot, label = "Potential energy", widen = true)
    plot!(kin[1:end], label = "Kinetic energy", widen = true) #, ylims = (1191.25, 1191.75))
    # plot!(virial[1:end], label = "Internal virial", widen = true, legend = :bottomright)#, ylims=(minimum(virial), 0))
    png(plotname(N, n, cu, T, new = "energies"))
    println(unstablerate(energy))
end

# Vssplit, Vscount = freqinintervals(Vs, n)
# bar()
# # allVs3 = round.(allVs, sigdigits = 3)
# # # # plotstart = 1
plotstates(states, frames = 100, new = join(("new cryst", n, cu, v₀), ' '))
# plot(energy[1:end], label = "Kinetic + potential energy", yguide = "E(tᵢ)", xguide = "tᵢ")#, ylims = (1232.54, 1232.545))
# # histogram(allVs, bins = 30, normalize = :pdf)
# V3fmap = stb.countmap(allVs)
begin
    Vs = cfs(countvabs2s, states)
    allVs = reduce(vcat, Vs)
    Vssplit, Vscount = freqinintervals(allVs, n)#; stop = 300)
end
bar(Vssplit[1:300], Vscount[1:300], bar_width = Vssplit[2] - Vssplit[1])

linemaxwell(V, fᵥ) = log(fᵥ) - log(V)
linVsmap = Dict(V => linemaxwell(V, fᵥ) for (V, fᵥ) in zip(Vssplit, Vscount))
# bar(linVsmap, bar_width = 0.1)
scatter(linVsmap, bar_width = Vssplit[2] - Vssplit[1], title = "Linearized Maxwell's distribution of velocities",
    label = "computed distribution")
# function lindistrV(Vs::Vector, n::Int; kwargs...)
#     Vssplit, Vscount = freqinintervals(Vs, n; kwargs...)
#     linVsmap = Dict(V => linemaxwell(V, fᵥ) for (V, fᵥ) in zip(Vssplit, Vscount))
#     scatter(linVsmap, bar_width = Vssplit[2] - Vssplit[1], title = "Linearized Maxwell's distribution of velocities",
#         label = "computed distribution")
#     # b, k = mnk(filter(isfinite ∘ last, linVsmap))
#     # plot!(Vssplit, k .* Vssplit .+ b, xlabel = "v²", ylabel = "ln[fᵥ(v)] - ln[v²]", label = "linear regression")
# end

# bar(Vssplit[1:30], Vscount[1:30], bar_width = Vssplit[2] - Vssplit[1])
# lindistrV(Vs[500], 50)
# savefig("Maxwell Linearized")

# n = 5000
# Vssplit, Vscount = freqinintervals(allVs, n; stop = 1e5)
# linVsmap = Dict(V => linemaxwell(V, fᵥ) for (V, fᵥ) in zip(Vssplit, Vscount))
# scatter(linVsmap)  # , bar_width = Vssplit[2] - Vssplit[1])
# b, k = mnk(filter(isfinite ∘ last, linVsmap))
# plot!(Vssplit, k .* Vssplit .+ b)

# png("Kin 2000 50 5")
# # plotstates(states, frames = 400, new="2withstart")

# # vir2 = cfs(countvirial, states)
# # plot(vir2, widen=true)



# # 3/2NkT = K

# function pressure(K::@Tr, W::@Tr, V::@Tr)
#     1 / 3V * (2K + W)
# end

# W = average(virial)
# K = average(kin)
# V = (2cu½) ^ 3
# P = pressure(K, W, V)
# press = pressure.(kin, virial, V)
# plot(press[1:end-40], label="Pressure", widen=true)#, ylims=(3.06e-5, 3.063e-5), xguide="tᵢ", yguide="P(tᵢ)")
# # png("Pressure 2000 50 5")

# temperature(K::@Tr, N::Number)::@Tr() = 2K / 3N
# temp = temperature.(kin, N)
# T = average(temp)
# σ = 3.465e-10
# ϵk = 116
# M=39.94
# k = 1.38e-23