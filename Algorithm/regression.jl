

function regressCoefficients(moments, C, D)
    iC = getIndices(moments, C)
    iD = getIndices(moments, D)
    #R_DC = round.(moments.S[iD,iC]*inv(moments.S[iC,iC]), digits=dgts)
    R_DC = moments.S[iD,iC]*inv(moments.S[iC,iC])
    return(R_DC)
end

function regressS(moments, M)
    return(M*(moments.S)*transpose(M))
end

function regressT(moments, M)
    Mt = transpose(M)
    n = size(M)[1]
    T₁₂₋₃ = reshape(moments.T,n*n,n)
    T₂ = transpose(kron(Mt,Mt))* T₁₂₋₃ * Mt
    return(reshape(T₂,n,n,n))
end

function regress(moments, C, D)
    if length(C) == 0
        return(moments)
    end
    if length(D) == 0
        D =  setdiff(moments.V, C)
    end

    n = length(moments.V)
    R = regressCoefficients(moments, C, D)

    M = diagm(ones(n))
    iC = getIndices(moments, C)
    iD = getIndices(moments, D)
    M[iD,iC] = round.(-R, digits = 10)

    S = regressS(moments, M)
    T = regressT(moments, M)

    return(Moments(moments.V, S, T, indices(moments.V)))
end

function residuals(moments, C, D)
    iD = getIndices(moments, D)
    mom_res = regress(moments, C, D)

    S = mom_res.S
    T = mom_res.T
    return(Moments(D, S[iD,iD], T[iD,iD,iD], indices(D)))
end

