

S2,(b2,b34,b4,b56,b6,b7)=PolynomialRing(QQ,["b2","b34","b4", "b56", "b6", "b7"])
R,(x, y, z)=PolynomialRing(S2,["x", "y", "z"])



#####
#TESTS 
#
#mon_k = (b4*b7*b6*b56^2,2)
#mon_all = ((b4*b7*b6*b56^2,2),(b4^2*b6*b56^2,2)) 
#cpoly = b2^2*b4*b7*b6
#cpoly2 = b2^2*b56^2*b6*b7
#filtered_monomial = 1
#unfiltered_monomial = 2
######

function cross_confirmation_faithful(monomial_and_k::Tuple{fmpq_mpoly,Int64},conf_vector)
    v = exponent_vector(monomial_and_k[1],1)
    out = Vector{Tuple{fmpq_mpoly,Any}}([])
        with_b34 = cross_confirmation_experimental(monomial_and_k,conf_vector[3],conf_vector[1],2)
        with_b2 = cross_confirmation_experimental(monomial_and_k,conf_vector[4],conf_vector[2],1)
    if with_b2[2] == 0
        push!(out,with_b2)
    else
        push!(out,(with_b2[1]*((gens(parent(with_b2[1]))[2])^(with_b2[2])),0))
    end
    if with_b34[2] == 0
        push!(out,with_b34)
    else
        push!(out,(with_b34[1]*((gens(parent(with_b34[1]))[1])^(with_b34[2])),0))
    end
    return out
end


##initial_terms with cross-confirmation, works for S b2,b34,b4,b56,b6,b7
function b4_dominant_terms(monomial_list::Vector{Tuple{fmpq_mpoly,Any}},conf_vector)
    exponents_of_monomial_list = [exponent_vector(i[1],1) for i in monomial_list]
    unique_b6b7_exponents = unique(map((j)->j[5:6],exponents_of_monomial_list))
    highest_exponents_b6_b7 = Vector{Tuple{fmpq_mpoly,Any}}([])
    out = Vector{Tuple{fmpq_mpoly,Any}}([])
    for un_exp in unique_b6b7_exponents
        terms_with_correct_monomial =  filter(i->(exponent_vector(i[1],1)[5:6]==un_exp),monomial_list)
        highest_b4 = maximum(i->(exponent_vector(i[1],1)[3]),terms_with_correct_monomial)
        potential_added_monomial_component = filter(i->(exponent_vector(i[1],1)[3]==highest_b4),terms_with_correct_monomial)
        lowest_b56 = minimum(i->(exponent_vector(i[1],1)[4]),potential_added_monomial_component)
        added_monomial_component = filter(i->(exponent_vector(i[1],1)[4]==lowest_b56),potential_added_monomial_component)[1]
        v = exponent_vector(added_monomial_component[1],1)
        k_old = added_monomial_component[2]
        k_new = k_old +v[1]+v[2]
        added_monomial = div(added_monomial_component[1],(gens(parent(added_monomial_component[1]))[1]^v[1])*(gens(parent(added_monomial_component[1]))[2]^v[2]))
        append!(out,cross_confirmation_faithful((added_monomial,k_new),conf_vector))
    end 
    return out
end


#function cone_identification_main(all_cones)
#    ineq = [inequalities_all_cones(all_cones)]
#    println(ineq)
#    out = Vector{Any}([])
#    for i in 1:length(ineq)
#        push!(out,(ineq[i][1], identify_all_potential_cones(ineq[i][2])))
#    end
#    return out
#end



function identify_all_potential_cones(entry)
    out = Vector{Any}([])
    for i in 1:length(entry[2])
        push!(out, (entry[2][i][1][1],inequalities_to_cone(entry[2][i][2])))
    end
    return out
end



function inequalities_to_cone(potential_inequalities)
    out = Vector{Any}([])
    for j in 1:length(potential_inequalities)
        if potential_inequalities[j] == "0 < 0"
            potential_inequalities[j] = ""
        elseif potential_inequalities[j] == "b56 < b4"
            push!(out, "1 or 3")
        elseif potential_inequalities[j] == "b4 < b56"
            push!(out, "2 or 4")
        elseif potential_inequalities[j] == "b34 < b2"
            push!(out, "a or c")
        elseif potential_inequalities[j] == "b2 < b36"
            push!(out, "b or d")
        elseif potential_inequalities[j] == "b4 + b7 < 2*b6"
            push!(out, "3 or 4")
        elseif potential_inequalities[j] == "b56 + b7 < 2*b6"
            push!(out, "3 or 4")
        elseif potential_inequalities[j] == "2*b6 < b4 + b7"
            push!(out, "1 or 2")
        elseif potential_inequalities[j] == "2*b6 < b56 + b7"
            push!(out, "1 or 2")
        elseif potential_inequalities[j] == "2*b4 < b2 + b6"
            push!(out, "a or b")
        elseif potential_inequalities[j] == "2*b4 < b34 + b6"
            push!(out, "a or b")
        elseif potential_inequalities[j] == "b2 + b6 < 2*b6"
            push!(out, "c or d")
        elseif potential_inequalities[j] == "b34 + b6 < 2*b6"
            push!(out, "c or d")
        end
    end
    return potential_inequalities, out
end


function inequalities_all_cones(output_of_initials)
    out = Vector{Any}([])
    for i in 1:length(output_of_initials)
        push!(out,(output_of_initials[i][1],identify_all_cones(output_of_initials[i][2])))
    end
    return out
end

function identify_all_cones(all_monomials)
    out = Vector{Any}([])
    for i in all_monomials
        push!(out, (i,identify_cone(i,all_monomials)))
    end
    return out
end



function identify_cone(monom,all_monomials)
    gen_s = gens(parent(monom[1]))
    v_monomial = exponent_vector(monom[1],1)
    potential_inequalities = Vector{String}([])
    for i in all_monomials 
        new_inequality = exponent_vector(i[1],1) - v_monomial 
        s = split_pos_neg(new_inequality)
        push!(potential_inequalities,exp_to_inequality(s[1],s[2],gen_s))
    end
    potential_inequalities = unique(potential_inequalities)
    return potential_inequalities
end

function exp_to_inequality(v1,v2,g)
    h = gcd(vcat(v1,v2))
    if h > 0
        v1 = Vector{Int64}((1/h)*v1)
        v2 = Vector{Int64}((1/h)*v2)
    end
    e1 = sum(g .* v1)
    e2 = sum(g.*v2)
    return("$e1 < $e2")
end

function exp_to_mon(v,g)
    return (g[1]^v[1])*(g[2]^v[2])*(g[3]^v[3])*(g[4]^v[4])*(g[5]^v[5])*(g[6]^v[6])
end


function split_pos_neg(v)
    v1 = Vector{Int64}(zeros(length(v)))
    v2 = Vector{Int64}(zeros(length(v)))
    for i in 1:length(v)
        if v[i] <0
            v1[i] = -v[i]
        else
            v2[i]= v[i]
        end
    end
    return [v1,v2]
end
        




function cross_confirmation_experimental(monomial_and_k::Tuple{fmpq_mpoly,Int64},cpoly::fmpq_mpoly,cpoly2::fmpq_mpoly,filtered_monomial::Int64)
    v = exponent_vector(monomial_and_k[1],1)
    k = monomial_and_k[2]
    all_monomials = collect(monomials(cpoly))
    expected_b65 = v[4]
    v[4]=0
    monomials_with_correct_exponents = filter(i->(exponent_vector(i,1)[3:6]==v[3:6]),all_monomials)
    if length(monomials_with_correct_exponents)>0
        highest_fm =min(maximum(i->(exponent_vector(i,1)[filtered_monomial]),monomials_with_correct_exponents),k)
        added_crossconv = filter(i->(exponent_vector(i,1)[filtered_monomial]==highest_fm),monomials_with_correct_exponents)
        if length(added_crossconv)==0
            highest_fm = highest_fm-2
            added_crossconv = filter(i->(exponent_vector(i,1)[filtered_monomial]==highest_fm),monomials_with_correct_exponents)
        end
        new_monomial = added_crossconv[1]*(gens(parent(cpoly))[4])^expected_b65
        k = k-highest_fm
        v2 = exponent_vector(new_monomial,1)
        expected_b4 = v2[3]
        v2[3] = 0
        all_monomials_b56 = collect(monomials(cpoly2))
        monomials_with_correct_exponents_b56 = filter(i->((exponent_vector(i,1)[3:6]==v2[3:6])&&exponent_vector(i,1)[filtered_monomial]==v2[filtered_monomial]),all_monomials_b56)
        if  length(monomials_with_correct_exponents)>0
            return (new_monomial,k)
        end
        return (new_monomial,(k,"ERROR:double verification failed"))
    end
    return monomial_and_k
end



function initial_terms_high(pyz::AbstractAlgebra.Generic.MPoly{fmpq_mpoly},confb56b34::AbstractAlgebra.Generic.MPoly{fmpq_mpoly},confb56b2::AbstractAlgebra.Generic.MPoly{fmpq_mpoly},confb4b34::AbstractAlgebra.Generic.MPoly{fmpq_mpoly},confb4b2::AbstractAlgebra.Generic.MPoly{fmpq_mpoly})
    monomial_list = collect(monomials(pyz))
    out =  Vector{Tuple{ AbstractAlgebra.Generic.MPoly{fmpq_mpoly}, Vector{Tuple{fmpq_mpoly,Any}}}}([])
    for i in 1:length(monomial_list)
        all_initial_forms = initial_terms_faithful_internal(coeff(pyz,monomial_list[i]), sum(exponent_vector(monomial_list[i],1)))
        processed_initial_forms = Vector{Tuple{fmpq_mpoly,Any}}([])
        conf_vector = [coeff(confb56b34,monomial_list[i]),coeff(confb56b2,monomial_list[i]),coeff(confb4b34,monomial_list[i]),coeff(confb4b2,monomial_list[i])]
        for j in all_initial_forms
            if j[2]>1
                append!(processed_initial_forms,cross_confirmation_faithful(j,conf_vector))
            else
                push!(processed_initial_forms,j)
            end
        end
        b4_dominant_forms = b4_dominant_terms(processed_initial_forms,conf_vector)
        append!(processed_initial_forms,b4_dominant_forms)
        push!(out,(monomial_list[i],unique(processed_initial_forms)))
    end
    return out
end






function initial_terms_faithful_internal(pyz::fmpq_mpoly,monomial_degree::Int64)
    
    expected_degree = max(1,7*(7-monomial_degree))
    out = []
    out2 = Vector{Tuple{fmpq_mpoly,Int64}}([])
    high_b7 = maximum(i->exponent_vector(i,1)[6], monomials(pyz)) 
    all_monomials = collect(monomials(pyz))
    contains_highb7 = filter(i->(exponent_vector(i,1)[6] == high_b7),all_monomials)
    high_b6= maximum(i->exponent_vector(i,1)[5], contains_highb7) 
    contains_highb6_and_highb7 = filter(i->(exponent_vector(i,1)[5] ==  high_b6),contains_highb7)
    high_b56 = maximum(i->exponent_vector(i,1)[4], contains_highb6_and_highb7) 
    high_b4 = maximum(i->exponent_vector(i,1)[3], contains_highb6_and_highb7) 
    v =[0,0,high_b4,high_b56,high_b6,high_b7]
    k = expected_degree - sum(v)
    append!(out,filter(i->(exponent_vector(i,1)[4] ==  high_b56),contains_highb6_and_highb7))
    append!(out,filter(i->(exponent_vector(i,1)[3] ==  high_b4),contains_highb6_and_highb7)) 

    initial_terms_internal(all_monomials,v,expected_degree,out,k)
    out = filter_bad_terms(out)
    for i in out
        k = expected_degree - sum(exponent_vector(i,1))
        push!(out2,(i,k))
    end
    return out2
end

function filter_bad_terms(out)
    out = unique(out)
    for i in out
        v = exponent_vector(i,1)
        filter!(j->(!(prod(exponent_vector(j,1).<=v)))||(j==i),out)
    end
    return out    
end


function initial_terms_internal(all_monomials,v,expected_degree,out,k)
    v1 = v+[0,0,-2,-2,4,-2]
    if (prod(v1[4:6] .>= 0) || (v1[3]>= 0 &&prod(v1[5:6].>=0)))&&v1[5]<=expected_degree 
        append!(out,filter(i->((exponent_vector(i,1)[4:6] == v1[4:6])||((exponent_vector(i,1)[5:6] == v1[5:6])&& (exponent_vector(i,1)[3] == v1[3])) ),all_monomials))
        initial_terms_internal(all_monomials,v1,expected_degree,out,k)
    end
   v2 = v+[0,0,0,4,-2,0]
   if (prod(v2[4:6] .>= 0) || (v2[3]>= 0 &&prod(v2[5:6].>=0)))&&v2[4]<=expected_degree 
        k = k-2
        if k>=0
            append!(out,filter(i->((exponent_vector(i,1)[4:6] == v2[4:6])||((exponent_vector(i,1)[5:6] == v2[5:6])&& (exponent_vector(i,1)[3] == v2[3])) ),all_monomials))
            initial_terms_internal(all_monomials,v2,expected_degree,out,k)
        end
    end
    v3 = v+[0,0,4,0,-2,0]
    if (prod(v3[4:6] .>= 0) || (v3[3]>= 0 &&prod(v3[5:6].>=0)))&&v3[3]<=expected_degree  
        k = k-2
        if k>=0
            append!(out,filter(i->((exponent_vector(i,1)[4:6] == v3[4:6])||((exponent_vector(i,1)[5:6] == v3[5:6])&& (exponent_vector(i,1)[3] == v3[3])) ),all_monomials))
            initial_terms_internal(all_monomials,v3,expected_degree,out,k)
        end
    end
    return out
end

