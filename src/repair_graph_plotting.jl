get_mid(lay_x, lay_y, src, dst) = (;x=(1.3lay_x[src]+0.7lay_x[dst])/2, y=(1.3lay_y[src]+0.7lay_y[dst])/2 )

function plt(g::AbstractGraph, repairdays::Vector{Int},day=1;seed=1, curvature_scalar=0.1, markersize = 0.4)
    @assert length(repairdays) == ne(g) "Graph has $(ne(g)) edges, but $(length(repairdays)) repair days are provided"
    repairs = zeros(Int, ne(g), maximum(repairdays))
    for (i,d) in enumerate(repairdays)
        if d > 0
            repairs[i,d] = 1
        end
    end
    plt(g, nothing, repairs, day;seed, curvature_scalar, markersize)
end

function plt(g::AbstractGraph, flows::Union{Matrix{Int}, Nothing}=nothing, repairs::Union{Matrix{Int}, Nothing}=nothing, day=1;seed=1, curvature_scalar=0.1, markersize = 0.4)
    if isnothing(flows) && isnothing(repairs)
        flows = zeros(Int, ne(g),day)
        repairs = zeros(Int, ne(g),day)
    elseif isnothing(flows)
        flows = zeros(Int, size(repairs))
    elseif isnothing(repairs)
        repairs = zeros(Int, size(flows))
    end
    @show flows, repairs
    plt(g, JuMP.Containers.DenseAxisArray(flows, Tuple.(edges(g)), 1:size(flows,2)), JuMP.Containers.DenseAxisArray(repairs, Tuple.(edges(g)), 1:size(repairs,2)), day;seed, curvature_scalar, markersize)
end


function plt(g::AbstractGraph, flows::JuMP.Containers.DenseAxisArray, repairs::JuMP.Containers.DenseAxisArray, day=1;seed=1, curvature_scalar=0.1, markersize = 0.4)
    Random.seed!(seed)
    lay = spring(g;seed)
    lay_x = first.(lay)
    lay_y = last.(lay)
    ms = fill(markersize,nv(g))
    ms[1] = 2*markersize
    ms[end] = 2*markersize
    node_weights = ms
    gp = graphplot(g; x=lay_x, y=lay_y, names=1:nv(g), curvature_scalar=curvature_scalar, markersize,node_weights,
    edgecolor = (s,d,w)-> (repairs[(s,d),day] < 1) ? :black : :red, edgestyle= (s,d,w)-> (repairs[(s,d),day]) < 1 ? :solid : :dash, trim=true)
    for e in edges(g)
        if flows[(e.src, e.dst), day] > 0
            mid = get_mid(lay_x, lay_y, e.src, e.dst)
            annotate!(mid.x, mid.y, text(flows[(e.src, e.dst),day], :blue))
        end
    end
    gp
end