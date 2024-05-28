### Utils
function cols(c)
    return collect(c[:,i] for i in 1:size(c)[2])
end

function normalize_groundset(M::Matroid)
    nb = Vector{Vector{Int64}}([])
    for b in bases(M)
        push!(nb, [findfirst(i->i==j,M.groundset) for j in b])
    end
    return matroid_from_bases(nb,length(M.groundset))
end

function matroid_from_polytope(P::Polyhedron)
    verts = Vector{Vector{QQFieldElem}}(collect(vertices(P)))
    bases = Vector{Vector{Int64}}([])
    for v in verts
        v = Int.(v)
        push!(bases,findall(i -> i == 1,v))
    end
    return matroid_from_bases(bases, ambient_dim(P))
end


function add_loop(M::Matroid)
    return matroid_from_bases(bases(M),M.groundset[length(M.groundset)]+1)
end

function add_coloop(m)
    bases_new = Vector{Vector{Int64}}([])
    newelem = length(m.groundset)+1
    for b in bases(m)
        b_new = push!(b,newelem)
        push!(bases_new, b_new)
    end
    return matroid_from_bases(bases_new,newelem)
end

function matroid_polytope(M::Matroid)
    b = bases(M)
    v = Vector{Vector{Int64}}([])
    for i in b
        push!(v,Vector{Int64}(map(j->(j in i),M.groundset)))
    end
    return convex_hull(transpose(reshape(collect(Iterators.flatten(v)), (length(v[1]),length(v)))))
end


function is_flag_matroid(M::Matroid,N::Matroid)
    return issubset(sort.(flats(M)),sort.(flats(N)))
end

## requires f to be a map from [n] to [m]
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


function is_morphism_of_matroids(M::Matroid, N::Matroid, f::Dict)
    return is_flag_matroid(induced_matroid(N,f),M)
end

function is_strong_map(M::Matroid, N::Matroid, f::Dict)
    loopmap = Dict((length(M.groundset)+1) => (length(N.groundset)+1))
    pointed_f = merge(f,loopmap)
    return is_flag_matroid(induced_matroid(add_loop(N),pointed_f),add_loop(M))
end

function find_all_morphisms(M::Matroid,N::Matroid)
    all_maps = all_maps_m_to_n(length(M.groundset),length(N.groundset),false)
    return filter((f) -> is_morphism_of_matroids(M,N,f) , all_maps)
end

function all_maps_m_to_n(m::Int,n::Int,sm::Bool)
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

function find_all_strong_maps(M::Matroid,N::Matroid)
    all_maps = all_maps_m_to_n(length(M.groundset),length(N.groundset),true)
    return filter((f) -> is_strong_map(M,N,f) , all_maps)
end

function filter_by_rank(r::Int,l::Vector{Dict{Int64,Int64}},N::Matroid)
    return filter((f)-> (rank(induced_matroid(add_loop(N),f))==r),l)
end



## Other utils:
## reshape(collect(Iterators.flatten(V)), (length(V[1]),length(V))) turns vector(vector) into matrix

#rks = [rank(M,ele) for ele in applied_map]
#println(applied_map)
#println(rks)
#return
#max_elems = findall(rks .== maximum(rks))
#highest_rank = unique!(applied_map[max_elems])
#new_bases = bases(uniform_matroid(highest_rank,new_groundset))


#lgths = [length(subs) for subs in highest_rank]
#smallest_lgths = findall(lgths .== minimum(lgths))
#new_bases = Vector{Vector{Int64}}([])
#images_bases = highest_rank[smallest_lgths]
#for im_b in images_bases
#    append!(new_bases,preimages_map(images_bases))
#end
#return matroid_from_bases(new_bases,new_groundset)


#ind = [is_morphism_of_matroids(M,N,f) for f in all_maps]
