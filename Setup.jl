import Base.product

function gencrystpos(N::Int, cu::@Tr)::@Tnv()
    a3 = ceil(Int, cbrt(N))
    a2 = cld(N, (a3^2))
    V = a3 * a3 * a2
    r = reshape([collect(t) for t in product(
            range(cu / a3, cu, a3),
            range(cu / a3, cu, a3),
            range(cu / a2, cu, a2))], :)
    if V > N
        splice!(r, 1:div(V, V - N):V)
    end
    r
end

function genr(N::Int, cu::@Tr)
    [rand(@Tr, 3) .* cu for _ = 1:N]
end

function genv(N::Int; T::@Tr() = 1.5)::@Tnv()
    Δvₓ = sqrt(T)
    [randn(@Tr, 3) .* Δvₓ for _ = 1:N]
end

function genv_v(N::Int; Δv₀::@Tr() = 2.12)::@Tnv()
    Δvₓ = Δv₀ / sqrt(3)
    [randn(@Tr, 3) .* Δvₓ for _ = 1:N]
end

function gencryst(N::Int, cu::@Tr; kwargs...)
    return gencrystpos(N, cu), genv(N; kwargs...)
end

function genrand(N::Int, cu::@Tr; kwargs...)
    return genr(N, cu), genv(N; kwargs...)
end

function gen2(cu½::@Tr, v₀::@Tr)
    ([[cu½/2 + 0.1, cu½, cu½], [3cu½/2 - 0.1, cu½, cu½]]), ([[v₀, 0, 0], [-v₀, 0, 0]])
end


function getprev(r_cur::@Tnv, v::@Tnv)
    r_cur - v * dt
end

function genwithprev(genwithv::Function, args...; kwargs...)
    r_cur, v = genwithv(args...; kwargs...)
    r_prev = getprev(r_cur, v)
    r_prev, r_cur
end


function setup!(ch, data_type = Nothing)
    global r, forces, forcessum, STEP
    if ch == 0
        r = genwithprev(gencryst, N, cu; T = T)
    elseif ch == 1
        r = genwithprev(genrand, N, cu; T = T)
    elseif ch == 2
        r = genwithprev(gen2, cu½, v₀)
    end
    if ch != -1
        STEP = UInt(0)
        r = (r..., [zeros(@Tr, 3) for _ = 1:N])
    end

    forces = Array{@Tv,2}(undef, N, N)
    forcessum = @Tnv()(undef, N)

    if data_type != Nothing
        Vector{data_type}(undef, nstart), Vector{data_type}(undef, n)
    end
end

narg(n, args...) = copy(args[n])
firstarg(args...) = narg(1, args...)
secondarg(args...) = narg(2, args...)

function processall(r_prev::@Tnv, r_cur::@Tnv)::Tuple{@Tnv,@Tr,@Tr}
    kin = kinenergy(r_prev, r_cur)
    virial = countvirial(r_prev, r_cur)
    state = copy(r_cur)
    state, kin, virial
end


function iterstart!(p::@Tr, r_prev::@Tnv, r_cur::@Tnv, r_next::@Tnv,
    forces::Array{@Tv,2}, forcessum::@Tnv)
    fullforcessum!(r_cur, forces, forcessum)
    for i in eachindex(r_cur)
        r_next[i] = verlestep(r_cur[i], r_prev[i], forcessum[i] * p)
    end
    nothing
end


function simulate(process; init = 1, datatype = Nothing)
    global r, STEP
    startdata, data = setup!(init, datatype)

    r = iterations!(iterstart!, (r), STEP, nstart, process, startdata, forces, forcessum)
    r = iterations!(iterverle!, (r), STEP, n, process, data, forces, forcessum)

    [startdata; data]
end

sen() = simulate(countenergy, datatype = @Tr)
svi() = simulate(countvirial, datatype = @Tr)
sf() = simulate(secondarg, datatype = @Tnv)
sall(init) = simulate(processall, init = init, datatype = Tuple{@Tnv,@Tr,@Tr})

function simulateall(init = 1, scale = 1)
    data = sall(init)[1:scale:end]
    states = Vector{@Tnv}(undef, length(data))
    kin = Vector{@Tr}(undef, length(data))
    virial = Vector{@Tr}(undef, length(data))
    for (i, entry) in enumerate(data)
        states[i], kin[i], virial[i] = entry
    end
    states, kin, virial
end
