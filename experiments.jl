using Pkg
pkg"activate ."

using OpenStreetMapX
using Graphs
using Parameters
using NetworkLayout
using LinearAlgebra

using RoadsConstructionOpt


const p = ModelSettings(N=10)

pth = "winnipeg_downtownf.osm"
#pth = joinpath(dirname(pathof(OpenStreetMapX)),"..","test","data","reno_east3.osm")

map_data =  get_map_data(pth;use_cache = false, trim_to_connected_graph=true );
sim = get_sim(map_data, p)


@time stats = run_simulation(sim)

#plot_edge_load(map_data,stats)

roads = top_congested_roads(sim,stats.vehicle_load,20)

#plot_edge_load_removed(map_data, stats, roads) #Removed roads colored green

reference_time=stats.simulation_total_time
#get_solution(sim,roads,reference_time,5, Dict{Vector{Tuple{Int,Int},Float64}())


ooo = @time opt(f, 2, 100, 0.001,roads, sim, reference_time)


### Toy model
a = zeros(Bool, 6,6)
a[1,2:3] .= 1
a[2,4:5] .= 1
a[3,4:5] .= 1
a[4,6] = 1
a[5,6] = 1
#a = LinearAlgebra.symmetric(a)
g = SimpleDiGraph(a)
# creating artificial map data

function get_map_data(g::AbstractGraph)
    graph_layout = spring(g)
    enus =  ENU.(first.(graph_layout)*1000, last.(graph_layout)*1000)
    refLLATemp = LLA(0.0,0.0)
    llamin = LLA(ENU(minimum(getX.(enus)), minimum(getY.(enus))), refLLATemp)
    llamax = LLA(ENU(maximum(getX.(enus)), maximum(getY.(enus))), refLLATemp)
    md_bounds = Bounds{LLA}(llamin.lat, llamax.lat, llamin.lon, llamax.lon)
    refLLA = OpenStreetMapX.center(md_bounds)
    # recode enus
    enus =  ENU.(LLA.(enus, Ref(refLLATemp)), Ref(refLLA))

    md_nodes = Dict(1:nv(g) .=> enus )
    md_roadways = Way.(1:ne(g))
    for (i, e) in enumerate(edges(g))
        append!(md_roadways[i].nodes, [src(e), dst(e)])
    end
    md_intersections = Dict(1:nv(g) .=> Set.(1:nv(g)))
    md_v = Dict{Int64, Int64}(1:nv(g) .=> 1:nv(g))
    md_e = [(src(e), dst(e)) for e in edges(g)]
    md_w = adjacency_matrix(g) .* 1000
    md_n = collect(1:nv(g))
    md_class =ones(Int, ne(g))
    MapData(md_bounds, md_nodes, md_roadways, md_intersections, g, md_v, md_n, md_e, md_w, md_class)
end
md = get_map_data(g)

psmall = ModelSettings(N=100)
sim2 =  get_sim(md, psmall; f_get_nodes=RoadsConstructionOpt.get_nodesfirstlast)

@time stats2 = run_simulation(sim2)


reference_time=stats2.simulation_total_time

roads2 = tuple.(src.(edges(g)), dst.(edges(g)))

ooo = @time opt(f, 3, 100, 0.001,roads2, sim2, reference_time)
