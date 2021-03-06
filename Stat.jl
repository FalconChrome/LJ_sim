using Plots
gr()

@userplot PointsPlot
@recipe function f(cp::PointsPlot)
    positions, boundaries = cp.args
    n = length(positions[1])
    linewidth --> 3
    xlims --> boundaries[1]
    ylims --> boundaries[2]
    zlims --> boundaries[3]
    aspect_ratio --> 1
    xguide --> "x"
    yguide --> "y"
    zguide --> "z"
    label --> false
    seriestype --> :scatter
    positions
end


function plotname(N, fps, n, acelllen, vā, new = nothing)
    filename_split = ["anim", N, "points", n, acelllen, vā, "fps", fps]
    if !isnothing(new)
        filename_split[end-1] = new
        insert!(filename_split, 2, 1)
        listdir = readdir()
        while join(filename_split, ' ') * ".gif" in listdir
            filename_split[6] += 1
        end
    end
    join(filename_split, ' ') * ".gif"
end


function plotstates(states; frames = 100, duration = 5, new = nothing)
    positions = [internal(states[i][j][k]) for j = 1:N, k = 1:3, i = 1:n]
    boundaries = fill((0, 2celllen), 3)

    fps = div(frames, duration)
    if fps > 30
        fps = 30
    end
    anim = @animate for i ā range(1, n, step = div(n, frames))
        print(i, '/', n, " of anima done")
        pointsplot(Tuple(positions[:, k, i] for k = 1:3), boundaries)
    end

    filename = plotname(N, fps, n, celllen, vā, new)
    gif(anim, filename, fps = fps)
end


function countfromstates(count::Function, states::Vector{@Tnv()})
    count.(states[1:end-1], states[2:end])
end

cfs = countfromstates


average(v::Vector) = length(v) > 0 ? sum(v) / length(v) : 0
average(f::Function, v::Vector) = length(v) > 0 ? sum(f(x) for x in v) / length(v) : 0

function stdd(values::Vector; average)
    if isempty(values)
        return 0
    end
    if isnothing(average)
        average = avg(values)
    end
    sqrt(sum((values .- average) .^ 2) / length(values))
end

function stdd(values::Vector, weights::Vector)
    if isempty(values)
        return 0
    end
    average = avg(values, weights)
    sqrt(sum((values .- average) .^ 2 .* weights) / sum(weights))
end


function unstablerate(values)
    a, b = extrema(values)
    (b - a) / (a + b)
end

mnk(data)::Tuple = mnk(keys(data), values(data))

function mnk(x, y)::Tuple
    sx = sum(x)
    sy = sum(y)
    n = length(x)
    k = (n * dot(x, y) - sx * sy) / (n * dot(x, x) - sx ^ 2)
    b = (sy - k * sx) / n
    b, k
end

function freqinintervals(v::Vector, n::Int;
    start::Number = -Inf, stop::Number = Inf)
    if !any(isfinite, (start, stop))
        start, stop = extrema(v)
    else
        if !isfinite(start)
            start = minimum(v)
        end
        if !isfinite(stop)
            stop = maximum(v)
        end
    end
    delimiters = range(start, stop, length = n + 1)
    s = step(delimiters)
    splitted = zeros(Int, n)
    for x in v
        if !(start <= x <= stop)
            continue
        end
        if x == stop
            splitted[end] += 1
            continue
        end
        index = findlast(t -> (x >= t), delimiters)
        if !isnothing(index)
            splitted[index] += 1
        end
    end
    intervals = [round(delimiters[i] + delimiters[i+1], sigdigits = round(Int, log10(n + 1) + 2)) / 2 for i = 1:n]
    intervals, splitted
end