

function root_nodes(moments, top_order, Λ, method_det, eps)
    if(length(moments.V) < 1) 
        return(moments, top_order, Λ)
    end

    C₁ = find_root_nodes(moments, method_det, eps)
    if length(C₁) == 0
        return(moments, top_order, Λ)
    end

    top_order = push!(top_order, C₁)
    mmts_res = residuals(moments, C₁, setdiff(moments.V, C₁))
    
    return(root_nodes(mmts_res, top_order, Λ, method_det, eps))
end

function find_root_nodes(moments, method_det, eps)
    
    D₂ = determinants₂(moments, method_det)
    # D₃ = determinants₃(moments, "svd")
    
    n = size(D₂)[1]
    d₂₀ = findall(r -> all(abs(Float64.(r)).<eps), D₂)
    # d₃₀ = findall(r -> all(abs(Float64.(r)).<eps), D₃)
    
    R = d₂₀; # intersect(d₂₀, d₃₀); 
    R = filter(r -> r[1]!=r[2], R)

    c = zeros(Bool, n)
    for i in 1:n
        rᵢ = filter(r -> r[1] == i, R)
        if length(rᵢ) == (n-1) c[i] = true end
    end

    C₁ = moments.V[findall(r -> r, c)]
    return(C₁)
end

