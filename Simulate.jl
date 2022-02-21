import StatsBase as stb
include("Particles.jl")
include("Setup.jl")
include("Stat.jl")

# Half the cell length
celllen = 50.0
cu = 2celllen
v₀ = 5.0
N = 100
n = 1000
nstart = 101

n_conc = N / cu^3
σ_int = pi * 2.5^2
l_fr = 1 / (n_conc * σ_int)

# energy = sen()
# states, energy = sf()
# virial = svi()
n = 400
nstart = 0
@time simulateall();
states, kin, virial = simulateall();

energy = cfs(countenergy, states)
# # # pot = potenergy.(states)

# # countv(r₁::@Tnv, r₂::@Tnv)::@Tnv() = countv.(r₁, r₂)
Vs = cfs(countvabs2s, states);
allVs = reduce(vcat, Vs)

average(Vs[end])
# allVs3 = round.(allVs, sigdigits = 3)
# # # plotstart = 1
# plotstates(states, frames = 100, new = join(("all", n, celllen, v₀), ' '))
plot(energy[1:end], label = "Kinetic + potential energy", yguide = "E(tᵢ)", xguide = "tᵢ")#, ylims = (1232.54, 1232.545))
# histogram(allVs, bins = 30, normalize = :pdf)
# V3fmap = stb.countmap(allVs)
Vssplit, Vscount = freqinintervals(allVs, 30; stop = 300)

linemaxwell(V, fᵥ) = log(fᵥ) - log(V)
linVsmap = Dict(V => linemaxwell(V, fᵥ) for (V, fᵥ) in zip(Vssplit, Vscount))
# bar(linV3map, bar_width = 0.1)

function lindistrV(Vs::Vector, n::Int; kwargs...)
    Vssplit, Vscount = freqinintervals(Vs, n; kwargs...)
    linVsmap = Dict(V => linemaxwell(V, fᵥ) for (V, fᵥ) in zip(Vssplit, Vscount))
    scatter(linVsmap, bar_width = Vssplit[2] - Vssplit[1], title = "Linearized Maxwell's distribution of velocities",
        label = "computed distribution")
    # b, k = mnk(filter(isfinite ∘ last, linVsmap))
    # plot!(Vssplit, k .* Vssplit .+ b, xlabel = "v²", ylabel = "ln[fᵥ(v)] - ln[v²]", label = "linear regression")
end

bar(Vssplit[1:30], Vscount[1:30], bar_width = Vssplit[2] - Vssplit[1])
lindistrV(Vs[500], 50)
savefig("Maxwell Linearized")

n = 5000
Vssplit, Vscount = freqinintervals(allVs, n; stop = 1e5)
linVsmap = Dict(V => linemaxwell(V, fᵥ) for (V, fᵥ) in zip(Vssplit, Vscount))
scatter(linVsmap)  # , bar_width = Vssplit[2] - Vssplit[1])
b, k = mnk(filter(isfinite ∘ last, linVsmap))
plot!(Vssplit, k .* Vssplit .+ b)

png("Kin 2000 50 5")
# # plot(pot)
plot(kin[1:end], label = "Kinetic energy", widen=true, ylims=(1191.25, 1191.75))
# plot(virial[1:1500], label="Internal virial", widen=true)#, ylims=(minimum(virial), 0))
unstablerate(energy)
# # plotstates(states, frames = 400, new="2withstart")

# # vir2 = cfs(countvirial, states)
# # plot(vir2, widen=true)



# # 3/2NkT = K

# function pressure(K::@Tr, W::@Tr, V::@Tr)
#     1 / 3V * (2K + W)
# end

# W = average(virial)
# K = average(kin)
# V = (2celllen) ^ 3
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