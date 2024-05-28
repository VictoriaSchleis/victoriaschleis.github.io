victoria@victoria-ThinkPad-L15-Gen-3:~$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.8.4 (2022-12-23)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> using Oscar
  ___   ____   ____    _    ____
 / _ \ / ___| / ___|  / \  |  _ \   |  Combining ANTIC, GAP, Polymake, Singular
| | | |\___ \| |     / _ \ | |_) |  |  Type "?Oscar" for more information
| |_| | ___) | |___ / ___ \|  _ <   |  Manual: https://docs.oscar-system.org
 \___/ |____/ \____/_/   \_\_| \_\  |  Version 1.0.0

julia> function cols(c)
           return collect(c[:,i] for i in 1:size(c)[2])
       end
cols (generic function with 1 method)

julia> function normalize_groundset(M::Matroid)
           nb = Vector{Vector{Int64}}([])
           for b in bases(M)
               push!(nb, [findfirst(i->i==j,M.groundset) for j in b])
           end
           return matroid_from_bases(nb,length(M.groundset))
       end
normalize_groundset (generic function with 1 method)

julia> function matroid_from_polytope(P::Polyhedron)
           verts = Vector{Vector{QQFieldElem}}(collect(vertices(P)))
           bases = Vector{Vector{Int64}}([])
           for v in verts
               v = Int.(v)
               push!(bases,findall(i -> i == 1,v))
           end
           return matroid_from_bases(bases, ambient_dim(P))
       end
matroid_from_polytope (generic function with 1 method)

julia> function add_loop(M::Matroid)
           return matroid_from_bases(bases(M),M.groundset[length(M.groundset)]+1)
       end
add_loop (generic function with 1 method)

julia> function add_coloop(m)
           bases_new = Vector{Vector{Int64}}([])
           newelem = length(m.groundset)+1
           for b in bases(m)
               b_new = push!(b,newelem)
               push!(bases_new, b_new)
           end
           return matroid_from_bases(bases_new,newelem)
       end
add_coloop (generic function with 1 method)

julia> function matroid_polytope(M::Matroid)
           b = bases(M)
           v = Vector{Vector{Int64}}([])
           for i in b
               push!(v,Vector{Int64}(map(j->(j in i),M.groundset)))
           end
           return convex_hull(transpose(reshape(collect(Iterators.flatten(v)), (length(v[1]),length(v)))))
       end
matroid_polytope (generic function with 1 method)

julia> function is_flag_matroid(M::Matroid,N::Matroid)
           return issubset(sort.(flats(M)),sort.(flats(N)))
       end
is_flag_matroid (generic function with 1 method)

julia> ## requires f to be a map from [n] to [m]
       function induced_matroid(M::Matroid,f::Dict)
           new_groundset = keys(f)
           im_groundset = unique!([get!(f,ele,0) for ele in new_groundset])
           new_rk = rank(M,im_groundset)
           b_m = bases(uniform_matroid(new_rk,length(new_groundset)))
           new_bases = Vector{Vector{Int64}}([])
           for b in b_m
               im_b = [get!(f,elem,0) for elem in b]
               if rank(M,im_b) == new_rk  
                   push!(new_bases,b)
               end
           end
           return matroid_from_bases(new_bases,new_groundset)
       end
induced_matroid (generic function with 1 method)

julia> function is_morphism_of_matroids(M::Matroid, N::Matroid, f::Dict)
           return is_flag_matroid(induced_matroid(N,f),M)
       end
is_morphism_of_matroids (generic function with 1 method)

julia> function is_strong_map(M::Matroid, N::Matroid, f::Dict)
           loopmap = Dict((length(M.groundset)+1) => (length(N.groundset)+1))
           pointed_f = merge(f,loopmap)
           return is_flag_matroid(induced_matroid(add_loop(N),pointed_f),add_loop(M))
       end
is_strong_map (generic function with 1 method)

julia> function find_all_morphisms(M::Matroid,N::Matroid)
           all_maps = all_maps_m_to_n(length(M.groundset),length(N.groundset),false)
           return filter((f) -> is_morphism_of_matroids(M,N,f) , all_maps)
       end
find_all_morphisms (generic function with 1 method)

julia> function all_maps_m_to_n(m::Int,n::Int,sm::Bool)
           n+=sm
           k = 1:n
           for i in 2:m
               k = Iterators.product(k,1:n)
           end
           c = collect(k)
           f = Iterators.flatten(c)
           for i in 2:m
               f = Iterators.flatten(f)
           end
           long_vec = collect(f)
           l = Int(length(long_vec)/m)
           mat_im = reshape(long_vec,m,l)
           if sm 
               add_row = reshape(collect(Iterators.flatten([n  for i in 1:l])), 1,l)
               m+=sm
               mat_im = vcat(mat_im,add_row)
           end
           images = cols(mat_im)
           return [Dict(zip(collect(1:m),img)) for img in images]
       end
all_maps_m_to_n (generic function with 1 method)

julia> function find_all_strong_maps(M::Matroid,N::Matroid)
           all_maps = all_maps_m_to_n(length(M.groundset),length(N.groundset),true)
           return filter((f) -> is_strong_map(M,N,f) , all_maps)
       end
find_all_strong_maps (generic function with 1 method)

julia> function filter_by_rank(r::Int,l::Vector{Dict{Int64,Int64}},N::Matroid)
           return filter((f)-> (rank(induced_matroid(add_loop(N),f))==r),l)
       end
filter_by_rank (generic function with 1 method)

julia> M = uniform_matroid(1,3)
Matroid of rank 1 on 3 elements

julia> flats(add_loop(M))
2-element Vector{Vector{Int64}}:
 [4]
 [1, 2, 3, 4]

julia> flats(add_coloop(M))
4-element Vector{Vector{Int64}}:
 []
 [4]
 [1, 2, 3]
 [1, 2, 3, 4]

julia> N = uniform_matroid(2,3)
Matroid of rank 2 on 3 elements

julia> flats(add_coloop(N))
10-element Vector{Vector{Int64}}:
 []
 [1]
 [2]
 [3]
 [4]
 [1, 2, 3]
 [1, 4]
 [2, 4]
 [3, 4]
 [1, 2, 3, 4]

julia> p14 = uniform_matroid(1,4)
Matroid of rank 1 on 4 elements

julia> p24 = uniform_matroid(2,4)
Matroid of rank 2 on 4 elements

julia> p14 = matroid_base_polytope(uniform_matroid(1,4))
Polyhedron in ambient dimension 4

julia> p24 = matroid_base_polytope(uniform_matroid(2,4))
Polyhedron in ambient dimension 4

julia> p34 = matroid_base_polytope(uniform_matroid(3,4))
Polyhedron in ambient dimension 4

julia> p13l = matroid_base_polytope(normalize_groundset(add_loop(uniform_matroid(1,3))))
Polyhedron in ambient dimension 4

julia> p23l = matroid_base_polytope(normalize_groundset(add_loop(uniform_matroid(2,3))))
Polyhedron in ambient dimension 4

julia> p13c = matroid_base_polytope(normalize_groundset(add_coloop(uniform_matroid(1,3))))
Polyhedron in ambient dimension 4

julia> p23c = matroid_base_polytope(normalize_groundset(add_coloop(uniform_matroid(2,3))))
Polyhedron in ambient dimension 4

julia> p33l = matroid_base_polytope(normalize_groundset(add_loop(uniform_matroid(3,3))))
Polyhedron in ambient dimension 4

julia> collect(vertic)
vertical_connectivity  vertical_direction     vertical_map           vertices               vertices_and_rays
julia> q = p14+p24+p34+p13c+p23c
Polyhedron in ambient dimension 4

julia> collect(vertices(q))
24-element Vector{PointVector{QQFieldElem}}:
 [5, 3, 1, 2]
 [5, 3, 0, 3]
 [5, 1, 3, 2]
 [5, 0, 3, 3]
 [5, 2, 0, 4]
 [5, 0, 2, 4]
 [3, 5, 1, 2]
 [3, 5, 0, 3]
 [1, 5, 3, 2]
 [0, 5, 3, 3]
 [2, 5, 0, 4]
 [0, 5, 2, 4]
 [3, 1, 5, 2]
 [3, 0, 5, 3]
 [1, 3, 5, 2]
 [0, 3, 5, 3]
 [2, 0, 5, 4]
 [0, 2, 5, 4]
 [4, 2, 0, 5]
 [4, 0, 2, 5]
 [2, 4, 0, 5]
 [0, 4, 2, 5]
 [2, 0, 4, 5]
 [0, 2, 4, 5]

julia> q = p13l+p23l+p33l+p13c+p23c
Polyhedron in ambient dimension 4

julia> collect(vertices(q))
6-element Vector{PointVector{QQFieldElem}}:
 [5, 3, 1, 2]
 [5, 1, 3, 2]
 [3, 5, 1, 2]
 [1, 5, 3, 2]
 [3, 1, 5, 2]
 [1, 3, 5, 2]

julia> boundary_lattice_points(q)
12-element SubObjectIterator{PointVector{ZZRingElem}}:
 [1, 3, 5, 2]
 [1, 4, 4, 2]
 [1, 5, 3, 2]
 [2, 2, 5, 2]
 [2, 5, 2, 2]
 [3, 1, 5, 2]
 [3, 5, 1, 2]
 [4, 1, 4, 2]
 [4, 4, 1, 2]
 [5, 1, 3, 2]
 [5, 2, 2, 2]
 [5, 3, 1, 2]

julia> q = p14+p24+p34+p13c+p23c
Polyhedron in ambient dimension 4

julia> boundary_lattice_points(q)
68-element SubObjectIterator{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 0, 4, 5]
 [2, 0, 5, 4]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 ⋮
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [4, 0, 2, 5]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 0, 5]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 3, 0, 3]
 [5, 3, 1, 2]

julia> q1 = boundary_lattice_points(q)
68-element SubObjectIterator{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 0, 4, 5]
 [2, 0, 5, 4]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 ⋮
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [4, 0, 2, 5]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 0, 5]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 3, 0, 3]
 [5, 3, 1, 2]


julia> q1 = collect(boundary_lattice_points(q))
68-element Vector{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 0, 4, 5]
 [2, 0, 5, 4]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 ⋮
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [4, 0, 2, 5]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 0, 5]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 3, 0, 3]
 [5, 3, 1, 2]

julia> q1 = boundary_lattice_points(q)
68-element SubObjectIterator{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 0, 4, 5]
 [2, 0, 5, 4]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 ⋮
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [4, 0, 2, 5]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 0, 5]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 3, 0, 3]
 [5, 3, 1, 2]


julia> filter!(i->sum(i)==11, collect(q1))
68-element Vector{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 0, 4, 5]
 [2, 0, 5, 4]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 ⋮
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [4, 0, 2, 5]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 0, 5]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 3, 0, 3]
 [5, 3, 1, 2]

julia> sort!.(q1)
68-element Vector{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 4, 5]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 3, 5]
 [0, 2, 4, 5]
 [0, 3, 4, 4]
 [0, 3, 4, 4]
 [0, 2, 4, 5]
 [0, 3, 3, 5]
 [1, 1, 4, 5]
 [1, 1, 4, 5]
 [1, 2, 3, 5]
 [1, 2, 3, 5]
 [1, 2, 3, 5]
 [1, 2, 3, 5]
 [1, 1, 4, 5]
 [1, 2, 4, 4]
 [1, 1, 4, 5]
 [1, 2, 3, 5]
 [1, 2, 3, 5]
 [0, 2, 4, 5]
 [0, 2, 4, 5]
 [1, 2, 3, 5]
 [1, 2, 3, 5]
 ⋮
 [0, 3, 4, 4]
 [2, 2, 3, 4]
 [0, 3, 3, 5]
 [1, 2, 3, 5]
 [0, 2, 4, 5]
 [0, 3, 4, 4]
 [0, 3, 4, 4]
 [1, 1, 4, 5]
 [1, 2, 4, 4]
 [0, 2, 4, 5]
 [2, 2, 3, 4]
 [0, 3, 4, 4]
 [2, 2, 3, 4]
 [0, 3, 4, 4]
 [1, 2, 4, 4]
 [0, 2, 4, 5]
 [0, 3, 3, 5]
 [1, 1, 4, 5]
 [1, 2, 3, 5]
 [1, 2, 3, 5]
 [0, 2, 4, 5]
 [1, 2, 3, 5]
 [2, 2, 2, 5]
 [0, 3, 3, 5]
 [1, 2, 3, 5]

julia> unique!(collect(q1))
68-element Vector{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 0, 4, 5]
 [2, 0, 5, 4]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 ⋮
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [4, 0, 2, 5]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 0, 5]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 3, 0, 3]
 [5, 3, 1, 2]

julia> unique!(collect(sort!.(q1)))
9-element Vector{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [1, 1, 4, 5]
 [1, 2, 3, 5]
 [1, 2, 4, 4]
 [2, 2, 2, 5]
 [2, 2, 3, 4]
 [2, 3, 3, 3]

julia> p = p14 +p24+p24+p34+p34
Polyhedron in ambient dimension 4

julia> q = collect(boundary_lattice_points(p))
88-element Vector{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 4, 5, 2]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [0, 5, 4, 2]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 5, 1]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [1, 5, 4, 1]
 [2, 0, 4, 5]
 ⋮
 [4, 3, 4, 0]
 [4, 4, 0, 3]
 [4, 4, 3, 0]
 [4, 5, 0, 2]
 [4, 5, 1, 1]
 [4, 5, 2, 0]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 0, 4, 2]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 1, 4, 1]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 2, 3, 1]
 [5, 2, 4, 0]
 [5, 3, 0, 3]
 [5, 3, 1, 2]
 [5, 3, 2, 1]
 [5, 3, 3, 0]
 [5, 4, 0, 2]
 [5, 4, 1, 1]
 [5, 4, 2, 0]

julia> v = collect(vertices(p))
24-element Vector{PointVector{QQFieldElem}}:
 [5, 4, 2, 0]
 [5, 4, 0, 2]
 [5, 2, 4, 0]
 [5, 0, 4, 2]
 [5, 2, 0, 4]
 [5, 0, 2, 4]
 [4, 5, 2, 0]
 [4, 5, 0, 2]
 [2, 5, 4, 0]
 [0, 5, 4, 2]
 [2, 5, 0, 4]
 [0, 5, 2, 4]
 [4, 2, 5, 0]
 [4, 0, 5, 2]
 [2, 4, 5, 0]
 [0, 4, 5, 2]
 [2, 0, 5, 4]
 [0, 2, 5, 4]
 [4, 2, 0, 5]
 [4, 0, 2, 5]
 [2, 4, 0, 5]
 [0, 4, 2, 5]
 [2, 0, 4, 5]
 [0, 2, 4, 5]

julia> q_not_v = filter(i->!(in(i,v)), q)
64-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 5, 1]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [1, 5, 4, 1]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 [2, 2, 2, 5]
 [2, 2, 5, 2]
 [2, 3, 1, 5]
 [2, 3, 5, 1]
 [2, 5, 1, 3]
 ⋮
 [3, 5, 1, 2]
 [3, 5, 2, 1]
 [3, 5, 3, 0]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 5, 1]
 [4, 3, 0, 4]
 [4, 3, 4, 0]
 [4, 4, 0, 3]
 [4, 4, 3, 0]
 [4, 5, 1, 1]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 1, 4, 1]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 2, 3, 1]
 [5, 3, 0, 3]
 [5, 3, 1, 2]
 [5, 3, 2, 1]
 [5, 3, 3, 0]
 [5, 4, 1, 1]

julia> transpose(q_not_v)
1×64 transpose(::Vector{PointVector{ZZRingElem}}) with eltype LinearAlgebra.Transpose{ZZRingElem, PointVector{ZZRingElem}}:
 [0 3 3 5]  [0 3 4 4]  [0 3 5 3]  [0 4 3 4]  [0 4 4 3]  [0 5 3 3]  [1 1 4 5]  [1 1 5 4]  [1 2 3 5]  …  [5 1 4 1]  [5 2 1 3]  [5 2 2 2]  [5 2 3 1]  [5 3 0 3]  [5 3 1 2]  [5 3 2 1]  [5 3 3 0]  [5 4 1 1]

julia> filter(i->in(5,i), q_not_v)
52-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 3, 5]
 [0, 3, 5, 3]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 5, 1]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [1, 5, 4, 1]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 [2, 2, 2, 5]
 [2, 2, 5, 2]
 [2, 3, 1, 5]
 [2, 3, 5, 1]
 [2, 5, 1, 3]
 [2, 5, 2, 2]
 [2, 5, 3, 1]
 [3, 0, 3, 5]
 ⋮
 [3, 1, 5, 2]
 [3, 2, 1, 5]
 [3, 2, 5, 1]
 [3, 3, 0, 5]
 [3, 3, 5, 0]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [3, 5, 2, 1]
 [3, 5, 3, 0]
 [4, 1, 1, 5]
 [4, 1, 5, 1]
 [4, 5, 1, 1]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 1, 4, 1]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 2, 3, 1]
 [5, 3, 0, 3]
 [5, 3, 1, 2]
 [5, 3, 2, 1]
 [5, 3, 3, 0]
 [5, 4, 1, 1]

julia> filter(i->!in(5,i), q_not_v)
12-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 4, 4]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [3, 0, 4, 4]
 [3, 4, 0, 4]
 [3, 4, 4, 0]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 3, 0, 4]
 [4, 3, 4, 0]
 [4, 4, 0, 3]
 [4, 4, 3, 0]

julia> filter(i->in(5,i)&in(4,i), q_not_v)
12-element Vector{PointVector{ZZRingElem}}:
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 4, 1, 5]
 [1, 4, 5, 1]
 [1, 5, 1, 4]
 [1, 5, 4, 1]
 [4, 1, 1, 5]
 [4, 1, 5, 1]
 [4, 5, 1, 1]
 [5, 1, 1, 4]
 [5, 1, 4, 1]
 [5, 4, 1, 1]

julia> filter(i->in(5,i)&in(3,i), q_not_v)
36-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 3, 5]
 [0, 3, 5, 3]
 [0, 5, 3, 3]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 [2, 3, 1, 5]
 [2, 3, 5, 1]
 [2, 5, 1, 3]
 [2, 5, 3, 1]
 [3, 0, 3, 5]
 [3, 0, 5, 3]
 [3, 1, 2, 5]
 [3, 1, 5, 2]
 [3, 2, 1, 5]
 [3, 2, 5, 1]
 [3, 3, 0, 5]
 [3, 3, 5, 0]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [3, 5, 2, 1]
 [3, 5, 3, 0]
 [5, 0, 3, 3]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 1, 3]
 [5, 2, 3, 1]
 [5, 3, 0, 3]
 [5, 3, 1, 2]
 [5, 3, 2, 1]
 [5, 3, 3, 0]

julia> p = p14 +p24+p13c+p23c+p34
Polyhedron in ambient dimension 4

julia> q = collect(boundary_lattice_points(p))
68-element Vector{PointVector{ZZRingElem}}:
 [0, 2, 4, 5]
 [0, 2, 5, 4]
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 3, 5, 3]
 [0, 4, 2, 5]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [0, 5, 2, 4]
 [0, 5, 3, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 3, 5, 2]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [1, 5, 3, 2]
 [2, 0, 4, 5]
 [2, 0, 5, 4]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 ⋮
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [3, 5, 0, 3]
 [3, 5, 1, 2]
 [4, 0, 2, 5]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 0, 5]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 0, 2, 4]
 [5, 0, 3, 3]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 1, 3, 2]
 [5, 2, 0, 4]
 [5, 2, 1, 3]
 [5, 2, 2, 2]
 [5, 3, 0, 3]
 [5, 3, 1, 2]

julia> v = collect(vertices(p))
24-element Vector{PointVector{QQFieldElem}}:
 [5, 3, 1, 2]
 [5, 3, 0, 3]
 [5, 1, 3, 2]
 [5, 0, 3, 3]
 [5, 2, 0, 4]
 [5, 0, 2, 4]
 [3, 5, 1, 2]
 [3, 5, 0, 3]
 [1, 5, 3, 2]
 [0, 5, 3, 3]
 [2, 5, 0, 4]
 [0, 5, 2, 4]
 [3, 1, 5, 2]
 [3, 0, 5, 3]
 [1, 3, 5, 2]
 [0, 3, 5, 3]
 [2, 0, 5, 4]
 [0, 2, 5, 4]
 [4, 2, 0, 5]
 [4, 0, 2, 5]
 [2, 4, 0, 5]
 [0, 4, 2, 5]
 [2, 0, 4, 5]
 [0, 2, 4, 5]

julia> q_not_v = filter(i->!(in(i,v)), q)
44-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 3, 5]
 [0, 3, 4, 4]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 4, 1, 5]
 [1, 4, 4, 2]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 [2, 2, 2, 5]
 [2, 2, 5, 2]
 [2, 3, 1, 5]
 [2, 3, 4, 2]
 [2, 4, 3, 2]
 [2, 5, 1, 3]
 [2, 5, 2, 2]
 [3, 0, 3, 5]
 [3, 0, 4, 4]
 [3, 1, 2, 5]
 [3, 2, 1, 5]
 [3, 2, 4, 2]
 [3, 3, 0, 5]
 [3, 3, 3, 2]
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 1, 5]
 [4, 1, 4, 2]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 2, 1, 3]
 [5, 2, 2, 2]

julia> filter(i->!in(5,i), q_not_v)
19-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 4, 4]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [1, 4, 4, 2]
 [2, 3, 4, 2]
 [2, 4, 3, 2]
 [3, 0, 4, 4]
 [3, 2, 4, 2]
 [3, 3, 3, 2]
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 4, 2]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]

julia> filter(i->in(5,i), q_not_v)
25-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 3, 5]
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 4, 1, 5]
 [1, 5, 1, 4]
 [1, 5, 2, 3]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 [2, 2, 2, 5]
 [2, 2, 5, 2]
 [2, 3, 1, 5]
 [2, 5, 1, 3]
 [2, 5, 2, 2]
 [3, 0, 3, 5]
 [3, 1, 2, 5]
 [3, 2, 1, 5]
 [3, 3, 0, 5]
 [4, 1, 1, 5]
 [5, 1, 1, 4]
 [5, 1, 2, 3]
 [5, 2, 1, 3]
 [5, 2, 2, 2]

julia> filter(i->in(5,i)&in(4,i), q_not_v)
6-element Vector{PointVector{ZZRingElem}}:
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 4, 1, 5]
 [1, 5, 1, 4]
 [4, 1, 1, 5]
 [5, 1, 1, 4]

julia> filter(i->in(5,i)&in(3,i), q_not_v)
15-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 3, 5]
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 5, 2, 3]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 [2, 3, 1, 5]
 [2, 5, 1, 3]
 [3, 0, 3, 5]
 [3, 1, 2, 5]
 [3, 2, 1, 5]
 [3, 3, 0, 5]
 [5, 1, 2, 3]
 [5, 2, 1, 3]

julia> filter(i->in(5,i)&in(3,i)&in(2,i), q_not_v)
12-element Vector{PointVector{ZZRingElem}}:
 [1, 2, 3, 5]
 [1, 2, 5, 3]
 [1, 3, 2, 5]
 [1, 5, 2, 3]
 [2, 1, 3, 5]
 [2, 1, 5, 3]
 [2, 3, 1, 5]
 [2, 5, 1, 3]
 [3, 1, 2, 5]
 [3, 2, 1, 5]
 [5, 1, 2, 3]
 [5, 2, 1, 3]

julia> filter(i->in(5,i)&in(3,i)&!in(2,i), q_not_v)
3-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 3, 5]
 [3, 0, 3, 5]
 [3, 3, 0, 5]

julia> filter(i->in(5,i)&!in(3,i), q_not_v)
10-element Vector{PointVector{ZZRingElem}}:
 [1, 1, 4, 5]
 [1, 1, 5, 4]
 [1, 4, 1, 5]
 [1, 5, 1, 4]
 [2, 2, 2, 5]
 [2, 2, 5, 2]
 [2, 5, 2, 2]
 [4, 1, 1, 5]
 [5, 1, 1, 4]
 [5, 2, 2, 2]

julia> filter(i->!in(5,i), q_not_v)
19-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 4, 4]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [1, 4, 4, 2]
 [2, 3, 4, 2]
 [2, 4, 3, 2]
 [3, 0, 4, 4]
 [3, 2, 4, 2]
 [3, 3, 3, 2]
 [3, 4, 0, 4]
 [3, 4, 2, 2]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 1, 4, 2]
 [4, 2, 3, 2]
 [4, 3, 0, 4]
 [4, 3, 2, 2]
 [4, 4, 0, 3]
 [4, 4, 1, 2]

julia> filter(i->!in(5,i)&!in(2,i), q_not_v)
9-element Vector{PointVector{ZZRingElem}}:
 [0, 3, 4, 4]
 [0, 4, 3, 4]
 [0, 4, 4, 3]
 [3, 0, 4, 4]
 [3, 4, 0, 4]
 [4, 0, 3, 4]
 [4, 0, 4, 3]
 [4, 3, 0, 4]
 [4, 4, 0, 3]

julia> filter(i->!in(5,i)&in(2,i), q_not_v)
10-element Vector{PointVector{ZZRingElem}}:
 [1, 4, 4, 2]
 [2, 3, 4, 2]
 [2, 4, 3, 2]
 [3, 2, 4, 2]
 [3, 3, 3, 2]
 [3, 4, 2, 2]
 [4, 1, 4, 2]
 [4, 2, 3, 2]
 [4, 3, 2, 2]
 [4, 4, 1, 2]

julia> filter(i->!in(5,i)&in(2,i)&in(3,i), q_not_v)
7-element Vector{PointVector{ZZRingElem}}:
 [2, 3, 4, 2]
 [2, 4, 3, 2]
 [3, 2, 4, 2]
 [3, 3, 3, 2]
 [3, 4, 2, 2]
 [4, 2, 3, 2]
 [4, 3, 2, 2]

julia> filter(i->!in(5,i)&in(2,i)&!in(3,i), q_not_v)
3-element Vector{PointVector{ZZRingElem}}:
 [1, 4, 4, 2]
 [4, 1, 4, 2]
 [4, 4, 1, 2]

julia> p = p13l +p23l+p13c+p23c+p33l
Polyhedron in ambient dimension 4

julia> q = collect(boundary_lattice_points(p))
12-element Vector{PointVector{ZZRingElem}}:
 [1, 3, 5, 2]
 [1, 4, 4, 2]
 [1, 5, 3, 2]
 [2, 2, 5, 2]
 [2, 5, 2, 2]
 [3, 1, 5, 2]
 [3, 5, 1, 2]
 [4, 1, 4, 2]
 [4, 4, 1, 2]
 [5, 1, 3, 2]
 [5, 2, 2, 2]
 [5, 3, 1, 2]

julia> v = collect(vertices(p))
6-element Vector{PointVector{QQFieldElem}}:
 [5, 3, 1, 2]
 [5, 1, 3, 2]
 [3, 5, 1, 2]
 [1, 5, 3, 2]
 [3, 1, 5, 2]
 [1, 3, 5, 2]

julia> q_not_v = filter(i->!(in(i,v)), q)
6-element Vector{PointVector{ZZRingElem}}:
 [1, 4, 4, 2]
 [2, 2, 5, 2]
 [2, 5, 2, 2]
 [4, 1, 4, 2]
 [4, 4, 1, 2]
 [5, 2, 2, 2]

julia> p = p13l +p23l+p24+p34+p33l
Polyhedron in ambient dimension 4

julia> q = collect(boundary_lattice_points(p))
43-element Vector{PointVector{ZZRingElem}}:
 [1, 3, 5, 2]
 [1, 4, 4, 2]
 [1, 4, 5, 1]
 [1, 5, 3, 2]
 [1, 5, 4, 1]
 [2, 2, 5, 2]
 [2, 3, 4, 2]
 [2, 3, 5, 1]
 [2, 4, 3, 2]
 [2, 4, 5, 0]
 [2, 5, 2, 2]
 [2, 5, 3, 1]
 [2, 5, 4, 0]
 [3, 1, 5, 2]
 [3, 2, 4, 2]
 [3, 2, 5, 1]
 [3, 3, 3, 2]
 [3, 3, 5, 0]
 [3, 4, 2, 2]
 [3, 4, 4, 0]
 [3, 5, 1, 2]
 [3, 5, 2, 1]
 [3, 5, 3, 0]
 [4, 1, 4, 2]
 [4, 1, 5, 1]
 [4, 2, 3, 2]
 [4, 2, 5, 0]
 [4, 3, 2, 2]
 [4, 3, 4, 0]
 [4, 4, 1, 2]
 [4, 4, 3, 0]
 [4, 5, 1, 1]
 [4, 5, 2, 0]
 [5, 1, 3, 2]
 [5, 1, 4, 1]
 [5, 2, 2, 2]
 [5, 2, 3, 1]
 [5, 2, 4, 0]
 [5, 3, 1, 2]
 [5, 3, 2, 1]
 [5, 3, 3, 0]
 [5, 4, 1, 1]
 [5, 4, 2, 0]

julia> v = collect(vertices(p))
18-element Vector{PointVector{QQFieldElem}}:
 [5, 4, 2, 0]
 [5, 4, 1, 1]
 [5, 3, 1, 2]
 [5, 2, 4, 0]
 [5, 1, 4, 1]
 [5, 1, 3, 2]
 [4, 5, 2, 0]
 [4, 5, 1, 1]
 [3, 5, 1, 2]
 [2, 5, 4, 0]
 [1, 5, 4, 1]
 [1, 5, 3, 2]
 [4, 2, 5, 0]
 [4, 1, 5, 1]
 [3, 1, 5, 2]
 [2, 4, 5, 0]
 [1, 4, 5, 1]
 [1, 3, 5, 2]

julia> q_not_v = filter(i->!(in(i,v)), q)
25-element Vector{PointVector{ZZRingElem}}:
 [1, 4, 4, 2]
 [2, 2, 5, 2]
 [2, 3, 4, 2]
 [2, 3, 5, 1]
 [2, 4, 3, 2]
 [2, 5, 2, 2]
 [2, 5, 3, 1]
 [3, 2, 4, 2]
 [3, 2, 5, 1]
 [3, 3, 3, 2]
 [3, 3, 5, 0]
 [3, 4, 2, 2]
 [3, 4, 4, 0]
 [3, 5, 2, 1]
 [3, 5, 3, 0]
 [4, 1, 4, 2]
 [4, 2, 3, 2]
 [4, 3, 2, 2]
 [4, 3, 4, 0]
 [4, 4, 1, 2]
 [4, 4, 3, 0]
 [5, 2, 2, 2]
 [5, 2, 3, 1]
 [5, 3, 2, 1]
 [5, 3, 3, 0]

julia> filter(i->!in(5,i), q_not_v)
13-element Vector{PointVector{ZZRingElem}}:
 [1, 4, 4, 2]
 [2, 3, 4, 2]
 [2, 4, 3, 2]
 [3, 2, 4, 2]
 [3, 3, 3, 2]
 [3, 4, 2, 2]
 [3, 4, 4, 0]
 [4, 1, 4, 2]
 [4, 2, 3, 2]
 [4, 3, 2, 2]
 [4, 3, 4, 0]
 [4, 4, 1, 2]
 [4, 4, 3, 0]

julia> filter(i->in(5,i), q_not_v)
12-element Vector{PointVector{ZZRingElem}}:
 [2, 2, 5, 2]
 [2, 3, 5, 1]
 [2, 5, 2, 2]
 [2, 5, 3, 1]
 [3, 2, 5, 1]
 [3, 3, 5, 0]
 [3, 5, 2, 1]
 [3, 5, 3, 0]
 [5, 2, 2, 2]
 [5, 2, 3, 1]
 [5, 3, 2, 1]
 [5, 3, 3, 0]

julia> filter(i->!in(5,i), q_not_v)
13-element Vector{PointVector{ZZRingElem}}:
 [1, 4, 4, 2]
 [2, 3, 4, 2]
 [2, 4, 3, 2]
 [3, 2, 4, 2]
 [3, 3, 3, 2]
 [3, 4, 2, 2]
 [3, 4, 4, 0]
 [4, 1, 4, 2]
 [4, 2, 3, 2]
 [4, 3, 2, 2]
 [4, 3, 4, 0]
 [4, 4, 1, 2]
 [4, 4, 3, 0]

