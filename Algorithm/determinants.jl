

function dᵤᵥ²(moments,u,v, method="det")
    S = moments.S
    T = moments.T
    M = [[S[u,u], S[u,v]] [T[u,u,u], T[u,u,v]]]
    if method == "det"
        return det(M)
    elseif method == "svd"
        sing_vals = svd(M).S
        return(minimum(sing_vals))
    elseif method == "cond_number"
        sing_vals = svd(M).S
        return(minimum(sing_vals)/maximum(sing_vals))
    end
end

function dᵤᵥ³(moments,u,v, method="det")
    S = moments.S
    T = moments.T
    M = [[S[u,u], S[u,v], S[v,v]] [T[u,u,u], T[u,u,v], T[u,v,v]] [T[u,u,v], T[u,v,v], T[v,v,v]]]
    if method == "det"
        return det(M)
    elseif method == "svd"
        sing_vals = svd(M).S
        return(minimum(sing_vals))
    elseif method == "cond_number"
        sing_vals = svd(M).S
        return(minimum(sing_vals)/maximum(sing_vals))
    end
end

function determinants₂(moments, method="svd")
    n = size(moments.S)[1]
    d₂ = zeros(n,n)
    [d₂[u,v] = dᵤᵥ²(moments,u,v, method) for u in 1:n for v in 1:n];
    return(d₂)
end

function determinants₃(moments, method="svd")
    n = size(moments.S)[1]
    d₃ = zeros(n,n)
    [d₃[u,v] = dᵤᵥ³(moments,u,v, method) for u in 1:n for v in 1:n];
    return(d₃)
end
