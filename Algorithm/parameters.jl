
function update_Λ_cycle(moments, Λ, C, vertices_ind)
    l = length(C)

    for i in 1:length(C)   
        x = mod(i,l)
        y = mod(i+1,l)
        z = mod(i+2,l)

        u = C[if x == 0 l else x end]
        v = C[if y == 0 l else y end]
        w = C[if z == 0 l else z end]

        uu = getIndices(moments,[u])[1]; vv = getIndices(moments,[v])[1]; ww = getIndices(moments,[w])[1]
        M = [[moments.S[uu,uu]    moments.S[uu,vv]    moments.S[uu,ww]     moments.S[vv,ww]]
             [moments.T[uu,uu,uu] moments.T[uu,uu,vv] moments.T[uu,uu,ww]  moments.T[uu,vv,ww]]
             [moments.T[uu,uu,vv] moments.T[uu,vv,vv] moments.T[uu,vv,ww]  moments.T[vv,vv,ww]]]

        iu = getIndices(vertices_ind,[u])[1]; iv = getIndices(vertices_ind,[v])[1]
        Λ[iu, iv] = det(M[:,[2,3,4]])/det(M[:,[1,3,4]])
    end
    return(Λ)
end

function update_Λ_sources(moments, top_order, Λ, dgts)

    for i in 2:length(top_order)
        #println(i)
        D = unique(reduce(vcat, top_order[i]))
        
        C = [reduce(vcat, top_order[j]) for j in 1:(i-1)]
        C = unique(reduce(vcat, C))

        R = regressCoefficients(moments, C, D)
        iC = getIndices(moments, C)
        iD = getIndices(moments, D)

        Λ[iC, iD] = round.(transpose(R)*(diagm(ones(length(iD))) - Λ[iD,iD]), digits=dgts)
    end
    return(Λ)
end