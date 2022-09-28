using LinearAlgebra

macro Tr()
    Float64
end
macro Tv()
    Vector{@Tr}
end
macro Tnv()
    Vector{@Tv}
end

const dt = 1e-3
const dt² = dt^2

function vectorabs2(v::@Tv)::@Tr
    sum(x^2 for x in v)
end

function vectorabs(v::@Tv)::@Tr
    sqrt(sum(v .^ 2))
end


function dist(r₁::@Tv, r₂::@Tv)::@Tr
    sqrt(sum((r₂ - r₁) .^ 2))
end

const zerovector = zeros(@Tr, 3)
iszerovector(v::@Tv) = all(iszero, v)

nearest(x::@Tr)::@Tr() = nearest(x, cu½)
nearest(Δr::@Tv)::@Tv() = nearest.(Δr)
function nearest(Δx::@Tr, cu½::@Tr)::@Tr
    x = abs(Δx)
    x <= cu½ ? Δx : copysign((x + cu½) % 2cu½ - cu½, Δx)
end

internal(Δr::@Tv)::@Tv() = map(internal, Δr)
internal(x::@Tr)::@Tr() = internal(x, cu)
function internal(x::@Tr, cu::@Tr)::@Tr
    if x >= cu
        return x % cu
    elseif x < 0
        return x % cu + cu
    end
    return x
end

function φ2(r2)
    if r2 != 0
        r6 = r2^-3
        return 4 * r6 * (r6 - 1)
    end
    return 0
end
φ(r) = φ2(r^2)
ϕ = φ

function fᵣtrue(Δr::@Tv)::@Tv
    Δr2 = vectorabs2(Δr)
    if Δr2 != 0
        Δr_2 = 1 / Δr2
        Δr_4 = Δr_2^2
        return (24 * (2 * Δr_4 * Δr_2 - 1) * Δr_4^2) .* Δr
    end
    return zeros(@Tr, 3)
end

function f(r)
    r_6 = r^-6
    24 * (2r_6 - 1) * r_6 / r
end

function forcesum(rᵢ::@Tv, r::@Tnv)::@Tv()
    sum((fᵣtrue ∘ nearest)(rᵢ - rⱼ) for rⱼ in r)
end

const cutoff_add = f(2.5) / 2.5
function fᵣ(Δr::@Tv)::@Tv
    Δr2 = vectorabs2(Δr)
    if 0 < Δr2  # < 6.25
        Δr_2 = 1 / Δr2
        Δr_4 = Δr_2^2
        return (24 * (2 * Δr_4 * Δr_2 - 1) * Δr_4^2 - 0) .* Δr
    end
    return zeros(@Tr, 3)
end

function force(rᵢ::@Tv, rⱼ::@Tv)::@Tv()
    (fᵣ ∘ nearest)(rᵢ - rⱼ)
end

function fullforcessum!(r_cur::@Tnv, forces::Array{@Tv,2}, forcessum::@Tnv)
    N = length(r_cur)
    for (i, rᵢ) in enumerate(r_cur)
        forces[1:i, i] = [force(rᵢ, rⱼ) for rⱼ in view(r_cur, 1:i)]
    end
    for i = 1:N
        @views forces[(i+1):end, i] = -forces[i, (i+1):N]
    end
    for i = 1:N
        forcessum[i] = sum(filter(!iszerovector, forces[:, i]), init = zerovector)
    end
    nothing
end

function verlestep(rₜ::@Tv, rₜ₋₁::@Tv, aₜ::@Tv)::@Tv()
    2 .* rₜ .- rₜ₋₁ .+ aₜ .* dt²
end

function iterverle!(p::Number, r_prev::@Tnv, r_cur::@Tnv, r_next::@Tnv,
    forces::Array{@Tv,2}, forcessum::@Tnv)
    fullforcessum!(r_cur, forces, forcessum)
    for i in eachindex(r_cur)
        r_next[i] = verlestep(r_cur[i], r_prev[i], forcessum[i])
    end
    nothing
end


function iterations!(iter!::Function, r::Tuple{@Tnv,@Tnv,@Tnv}, STEP::UInt,n::Int,
    process::Function, savedata::Vector, iter_args...)
    println(n, " steps:")
    for i = 1:n
        iter!(i / n, r..., iter_args...)
        STEP += 1
        savedata[i] = process(r[2], r[3])
        r = r[2], r[3], r[1]
        if i % 100 == 0
            println(i)
        end
    end
    r
end



# energy computing
countv(r₁₂::Tuple{@Tv,@Tv})::@Tv() = countv(r₁₂...)
countv(r₁::@Tv, r₂::@Tv)::@Tv() = countv(r₂ - r₁)
countv(Δr::@Tv)::@Tv() = Δr / dt

countvabs2 = vectorabs2 ∘ countv
countvabs = vectorabs ∘ countv

countvabs2s(r_prev::@Tnv, r_cur::@Tnv)::@Tv() = broadcast(vectorabs2 ∘ countv, r_prev, r_cur)
countvabss(r_prev::@Tnv, r_cur::@Tnv)::@Tv() = broadcast(vectorabs ∘ countv, r_prev, r_cur)

function potenergy(r::@Tnv)::@Tr
    if isempty(r)
        return 0
    end
    sum((φ2 ∘ vectorabs2 ∘ nearest)(rᵢ - rⱼ) for rⱼ in r, rᵢ in r) / 2 / length(r)
end

function kinenergy(r_prev::@Tnv, r_cur::@Tnv)::@Tr
    if isempty(r_cur)
        return 0
    end
    sum(vectorabs2 ∘ countv, zip(r_prev, r_cur)) / 2 / length(r_cur)
end

function countvirial(r_prev::@Tnv, r_cur::@Tnv)::@Tr
    sum(dot(rᵢ, forcesum(rᵢ, r_cur)) for rᵢ in r_cur) / length(r_cur)
end

function countenergy(r_prev::@Tnv, r_cur::@Tnv)::@Tr
    potenergy(r_cur) + kinenergy(r_prev, r_cur)
end

# function temperature(K::@Tr)::@Tr

# end
