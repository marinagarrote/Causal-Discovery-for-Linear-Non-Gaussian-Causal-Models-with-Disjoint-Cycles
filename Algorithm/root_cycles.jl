


function root_cycles(moments, top_order, Λ, vertices_ind, method_det, eps)
    if length(moments.V) < 3
        return(moments, top_order, Λ, vertices_ind)
    end

    C₂ = find_root_cycles(moments, method_det, eps)
    if length(C₂) == 0 return(moments, top_order, Λ, vertices_ind) end
    C₂ = join_incompatible_cycles(C₂)

    for c in 1:length(C₂)
        C₂[c] = sortNodesCycle(moments, C₂[c])
        Λ = update_Λ_cycle(moments, Λ, C₂[c], vertices_ind)
    end

    top_order = push!(top_order, C₂)

    setC₂ = unique(reduce(vcat, C₂))
    mmts_res = residuals(moments, setC₂, setdiff(moments.V, setC₂))

    return(mmts_res, top_order, Λ, vertices_ind)
    #return(root_cycles(mmts_res, top_order, Λ, vertices_ind))
end


function join_incompatible_cycles(C)
    if length(C) == 1
        return C
    end
    
    newC = C
    idx=[]
    for i in 1:length(C)
        for j in (i+1):length(C)
            if length(intersect(newC[i], newC[j])) != 0
                newC[j] = unique(vcat(newC[i], newC[j]))
                idx = append!(idx, i)
                continue
            end
        end
    end
    return newC[1:end .∉ [idx]]
end

function find_root_cycles(moments, method_det, eps)
    C₂ = candidate_root_cycles(moments, method_det, eps)
    if length(C₂) > 1 
        C₂ = unique(identifyRootCycles(moments, C₂, eps)) 
    end
    return C₂
end

function candidate_root_cycles(moments, method_det, eps)
    D₂ = determinants₂(moments, method_det)
    D₃ = determinants₃(moments, method_det)
    
    n = size(D₂)[1]
    d₂₀ = findall(r -> all(abs(Float64.(r)).>eps), D₂)
    d₃₀ = findall(r -> all(abs(Float64.(r)).<eps), D₃)
    
    d₂₀ = filter(r -> r[1]!=r[2], d₂₀)
    d₃₀ = filter(r -> r[1]!=r[2], d₃₀)

    R = intersect(d₂₀, d₃₀)

    g = SimpleGraph(n,0)
    while length(R) > 1
        i = 1
        if CartesianIndex(R[i][2], R[i][1]) in R
            if !has_edge(g,R[i][2], R[i][1])
                add_edge!(g, R[i][1], R[i][2])
            end
            remove!(R,CartesianIndex(R[i][2], R[i][1]))
        end
        remove!(R, R[i])
    end
        
    maxCliques = maximal_cliques(g)
    iV = filter(r -> length(r) >= 3, maxCliques)
    
    return([moments.V[iV[s]] for s in 1:length(iV)])
end

function identifyRootCycles(moments, C₂, eps)
    childCycles = []
    for i in 1:length(C₂)
        for j in (i+1):length(C₂)
            cc = childCycle(moments, C₂[i], C₂[j], eps)
            if length(cc) != 0 childCycles = push!(childCycles, cc) end
        end
    end
    return(setdiff(C₂, childCycles))
end

function childCycle(moments, D₁, D₂, eps)
    m₁ = moments
    m₂ = moments
    m₁ = regress(m₁, D₁, D₂);
    m₂ = regress(m₂, D₂, D₁)
   
    iD₁ = getIndices(moments, D₁)
    iD₂ = getIndices(moments, D₂)

    l₁ = length(findall(r -> all(abs(r) .< eps), m₁.T[iD₁, iD₁, iD₂]))
    l₂ = length(findall(r -> all(abs(r) .< eps), m₂.T[iD₂, iD₂, iD₁]))
    d₁ = length(D₁)^2*length(D₂)
    d₂ = length(D₂)^2*length(D₁)

    if (l₁ == d₁) & (l₂ != d₂)
        return(D₂)
    elseif (l₂ == d₂) & (l₁ != d₁)
        return(D₁)
    else
        return()
    end
end

function sortNodesCycle(moments, C)
    C = sort(C)
    #println(C)
    v = C[1]
    cycle = [v]
   
    while length(C) > 1
        x = setdiff(vcat(cycle,C), v)
        dets = []
        for i in setdiff(C, v)
            y = setdiff(vcat(cycle,C), i)
            ix = getIndices(moments, x)
            iy = getIndices(moments, y)
            sing_vals = svd(moments.S[ix,iy]).S
            d = minimum(sing_vals)/maximum(sing_vals)
            # d = det(moments.S[ix,iy])
            dets = append!(dets, d)
        end
        vi = setdiff(C, v)[argmax(dets)]
        push!(cycle, vi)
        C = setdiff(C,v)
        v = vi
        # print(d)
        #     if abs(d) > eps  
        #         push!(cycle, i) 
        #         C = setdiff(C,v)
        #         v = i
        #         break
        #     end
    end

    ## PROVISIONAL
    sc = sum([cycle[i-1] < cycle[i] for i in 2:length(cycle)])
    revcycle = cycle[end:-1:1]
    src = sum([revcycle[i-1] < revcycle[i] for i in 2:length(cycle)])

    if sc > src
        return cycle
    else
        return revcycle
    end
end