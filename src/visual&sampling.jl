using StatsBase

using Pkg
#1. Using OSM tutorial and folium visualisation

Pkg.add("PyCall")
Pkg.add("Conda")
Pkg.add("OpenStreetMapX")
using Conda
Conda.runconda(`install folium -c conda-forge`)

mx = get_map_data(map_file_path, use_cache=false);
using Random 
Random.seed!(0)
node_ids = collect(keys(mx.nodes)) 
routes = Vector{Vector{Int}}()
visits = Dict{Int,Int}()
for i in 1:10000
    a,b = [point_to_nodes(generate_point_in_bounds(mx), mx) for _ in 1:2]
    route, route_time = OpenStreetMapX.shortest_route(mx,a,b)
    if route_time < Inf # when we select points neaer edges no route might be found
        push!(routes, route)
        for n in route
            visits[n] = get(visits, n,0)+1
        end 
    end
end                                   
println("We have generated ",length(routes)," non-empty routes")

#2. Visualization of most frequently visited intersections
using PyCall
flm = pyimport("folium")
matplotlib_cm = pyimport("matplotlib.cm")
matplotlib_colors = pyimport("matplotlib.colors")
cmap = matplotlib_cm.get_cmap("prism")

m = flm.Map()

max_visits= maximum(values(visits))
for k=keys(visits)
    visits[k] < 25 && continue  #skip nodes infrequently visited
    loc = LLA(mx.nodes[k],mx.bounds)
    info = "Node $(k) at ($(round(loc.lat,digits=4)), $(round(loc.lon,digits=4)))\n<br>"*
           "visits: $(visits[k])"
    flm.Circle(
      location=[loc.lat,loc.lon],
      popup=info,
      tooltip=info,
      radius=500*visits[k]/max_visits,
      color="crimson",
      weight=0.5,
      fill=true,
      fill_color="crimson"
   ).add_to(m)
end
MAP_BOUNDS = [(mx.bounds.min_y,mx.bounds.min_x),(mx.bounds.max_y,mx.bounds.max_x)]
flm.Rectangle(MAP_BOUNDS, color="black",weight=6).add_to(m)
m.fit_bounds(MAP_BOUNDS)
m


#3. Weighted sampling of intersections.

prob = Dict{Real,Real}() 

for k=keys(visits)
    prob[k]=visits[k]/sum(values(visits))
end

weights = collect(values(prob)) #probabilities 

intersections = collect(keys(visits)) #intersections in the map

StatsBase.direct_sample!(intersections,weights)[1:10] #Example - choosing 10 roads
