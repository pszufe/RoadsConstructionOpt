using Pkg; Pkg.activate(".")

using JuMP, Graphs, Ipopt, InvertedIndices, Juniper, Random
import JuMP.Containers.DenseAxisArray
using GraphRecipes, NetworkLayout
using Plots, Random, Statistics

using RoadsConstructionOpt

function plan_repairs(g::SimpleDiGraph;
     # set of days
     T = 1:3,
     # start node
     A=1,
     #destination node
     B=nv(g),
     # total flow from A to B (A and B  are  the same for all days and all cars)
     μ = 1_000,
     # maximum speed on each edge (the same od all edges)
     vmax = 100,
     # minimum speed on each edge
     vmin = 10,
    # set of edges
    E = Tuple.(edges(g)),
     # repair plan - here we assume that all edges are requiring a repair
     θ = DenseAxisArray(ones(Int, length(E)), E)

    )
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
    optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt)
    m = Model(optimizer)



    # maximum number of cars on each edge
    N = DenseAxisArray(fill(μ, length(E)), E)
    #edge distances (the same for all edges = 1)
    d = DenseAxisArray(fill(1, length(E)), E)
    # repair plan - binary variable x[e,u] = 1 if edge e is repaired on day u
    @variable(m, x[E, T], Bin)
    # flow on edge e on day u (number of cars)
    @variable(m, 0 <= n[E, T] <= μ) # , Int can be dropped when the traffic flows are high

    # excatly μ cars are leaving A
    @constraint(m, [u in T], sum(n[(A,i),u] for i in outneighbors(g,A)) == μ)
    # exactly μ cars are arriving at B
    @constraint(m, [u in T], sum(n[(i,B),u] for i in inneighbors(g,B)) == μ)

    #balancing car flows at each vetex
    @constraint(m, [u in T, i in vertices(g)[Not(A,B)]], sum(n[(j,i),u] for j in inneighbors(g,i)) == sum(n[(i,j),u] for j in outneighbors(g,i)))

    #Flow on each edge cannot exceed its capacity
    @constraint(m, [e in E, u in T], n[e,u] <= N[e])
    #the edge e needs to repaired θ times over the time horizon T
    @constraint(m, [e in E], sum(x[e,u] for u in T) == θ[e])
    # if edge e is repaired on day u, then the flow on this edge is 0
    @constraint(m, [e in E, u in T], n[e,u] <= (1-x[e,u])*μ )

    # total time in traffic spend by all cars over the time horizon T
    @NLobjective(m, Min, sum(d[e] * n[e,u]/(vmax - (n[e,u]/N[e])*(vmax-vmin)) for e in E, u in T))

    #running the optimization
    optimize!(m)

    #extracting the results
    n_flows = round.(Int, value.(n))
    x_repairs = round.(Int, value.(x))

    (;n_flows, x_repairs, obj=objective_value(m))
end





#adjacency matrix of a smaller graph
a = zeros(Bool, 6,6)
a[1,2:3] .= 1
a[2,4:5] .= 1
a[3,4:5] .= 1
a[4,6] = 1
a[5,6] = 1


#adjacency matrix of a larger graph
b = zeros(Bool, 10,10)
b[1,2:4] .= 1
b[4,5] = 1
b[3,5:6] .= 1
b[2,6] = 1
b[5,6:8] .= 1
b[6,[5,9]] .= 1
b[7,[8,10]] .= 1
b[8,[7,9,10]] .= 1
b[9,[8,10]] .= 1

c =b[1:8,1:8]
b[6,8] = 1
b[3,7] = 1

gs = Dict(:small=>SimpleDiGraph(a), :large=>SimpleDiGraph(b), :medium=>SimpleDiGraph(c))

#graphplot(SimpleDiGraph(c);names=1:8, curvature_scalar=0.1, markersize=0.4)


n_flows, x_repairs, obj = plan_repairs(gs[:small])

plt(gs[:small], n_flows, x_repairs, 1;seed=15, curvature_scalar=0.1, markersize=0.9)

elaps = Dict()
u = nothing
for gsize in [:small, :large, :medium]
    seed = 15
    ttt = @elapsed n_flows, x_repairs, obj = plan_repairs(gs[gsize])
    elaps[gsize] = ttt
    println("$gsize Total time in traffic ", obj)

    println("$gsize Flow on each edge on each day:\n ", n_flows)

    println("Repair plan for each edge for each day:\n ", x_repairs)
    for day in 1:3
        plt(gs[gsize], n_flows, x_repairs, day;seed, curvature_scalar=0.1, markersize=0.9)
        savefig("$(gsize)_$day.pdf")
        savefig("../5e923468bfa62e0001281297/pictures/$(gsize)_$day.pdf")
    end

    E = Tuple.(edges(gs[gsize]))
    θ = DenseAxisArray(zeros(Int, length(E)), E)
    @time n_flowsz, x_repairsz, objz = plan_repairs(gs[gsize];θ)
    plt(gs[gsize], n_flowsz, x_repairsz, 1;seed, curvature_scalar=0.1, markersize=0.9)
    savefig("toy$(gsize).pdf")
    savefig("../5e923468bfa62e0001281297/pictures/toy$gsize.pdf")
end


println(elaps)
