
function quiver_pluecker_relations(r::Int, s::Int, A::Matrix{Int})
    n,m = size(A)
    Iset = subsets(collect(1:m),r-1)
    Jset = subsets(collect(1:n),s+1)
    m_c_r =subsets(collect(1:m),r)
    n_c_s =subsets(collect(1:n),s)
    R,x,y = polynomial_ring(QQ, "x"=>m_c_r, "y"=>n_c_s)
    pluecker_list = Vector{MPolyRingElem}([])
    for I in Iset
        Icomp = filter!(k->!(in(k,I)), collect(1:m))
        for J in Jset
            f = R(0)
            for j in Icomp
            sj = sign_qpl(j,I,J)
                for i in J
                    I_c_j = sort([I;j]);
                    J_m_i = setdiff(J,[i]); 
                    f = f + sj*A[i,j]*x[findfirst(k -> k == I_c_j, m_c_r)]*y[findfirst(l -> l == J_m_i, n_c_s)]
                end
            end
            push!(pluecker_list,f)   
        end
    end
    return unique!(pluecker_list)
end

function sign_qpl(j::Int,I::Vector{Int},J::Vector{Int})
    return (-1)^(length(findall(k -> k>j, J))+length(findall(i -> i>j, I)))
end


function normalize_polynomial(f::MPolyRingElem)
    coeffs_f = collect(coefficients(f))
    f = f/gcd(coeffs_f)
    if coeffs_f[findfirst(i->i != 0,coeffs_f)]<0
        f = (-1)*f
    end
    return f
end

function remove_redundant_polys(p1::Vector{MPolyRingElem})
    unique!(filter!(i->i != 0, p1))    
    return unique!([normalize_polynomial(i) for i in p1])
end

#### Example loop quiver
n,m,r,s = (5,5,2,2);
Iset = subsets(collect(1:m),r-1);
Jset = subsets(collect(1:n),s+1);
m_c_r =subsets(collect(1:m),r);
n_c_s =subsets(collect(1:n),s);
R,x,y = polynomial_ring(QQ, "x"=>m_c_r, "y"=>n_c_s)
M = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 1 0] ##contraction matrix
B = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0] ## projection matrix
p1 = quiver_pluecker_relations(r,s,M)
p1 = [evaluate(i,vcat(x,x)) for i in p1]
p1 = remove_redundant_polys(p1)
I = ideal(p1)
dim(I)
primary_decomposition(I)
