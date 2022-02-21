
function genrand(N::Int, celllen::@Tr, v₀::@Tr)::Tuple{@Tnv,@Tnv}
    ([rand(@Tr, 3) .* (2 * celllen) for _ = 1:N],
        [(rand(@Tr, 3) .- 0.5) .* (2 * v₀) for _ = 1:N])

end


function gen2(celllen::@Tr, v₀::@Tr)::Tuple{@Tnv,@Tnv}
    ([[celllen / 2 + 0.1, celllen, celllen], [3celllen / 2 - 0.1, celllen, celllen]]), ([[v₀, 0, 0], [-v₀, 0, 0]])
end

function getprev(r_cur::@Tnv, v::@Tnv)::@Tnv
    r_cur - v * dt
end

function genwithprev(genwithv::Function, args...)
    r_cur, v = genwithv(args...)
    r_prev = getprev(r_cur, v)
    r_prev, r_cur
end


function setup!(ch, data_type = Nothing)
    global r, forces, forcessum
    if ch == 2
        r = genwithprev(gen2, celllen, v₀)
    else
        r = genwithprev(genrand, N, celllen, v₀)
    end
    r = (r..., [zeros(@Tr, 3) for _ = 1:N])

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
    global r
    startdata, data = setup!(init, datatype)

    r = iterations!(iterstart!, (r), nstart, process, startdata, forces, forcessum)
    r = iterations!(iterverle!, (r), n, process, data, forces, forcessum)

    [startdata; data]
end

sen() = simulate(countenergy, datatype = @Tr)
svi() = simulate(countvirial, datatype = @Tr)
sf() = simulate(secondarg, datatype = @Tnv)
sall() = simulate(processall, datatype = Tuple{@Tnv,@Tr,@Tr})

function simulateall(scale = 1)
    data = sall()[1:scale:end]
    states = Vector{@Tnv}(undef, length(data))
    kin = Vector{@Tr}(undef, length(data))
    virial = Vector{@Tr}(undef, length(data))
    for (i, entry) in enumerate(data)
        states[i], kin[i], virial[i] = entry
    end
    states, kin, virial
end
