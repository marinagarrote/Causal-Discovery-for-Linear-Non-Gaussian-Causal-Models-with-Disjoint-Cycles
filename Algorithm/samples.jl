using Random, Distributions

function random_ϵ(n, size)
    ϵ = rand(TriangularDist(-1, 1, 1), size, n)
    #ϵ = rand(Exponential(1), size, n)
    # for i in 1:size
    #     ϵᵢ = rand(TriangularDist(-1, 1, 1), n, m)
    #     # ϵᵢ = rand(Exponential(1), n)
    #     ϵ = push!(ϵ, ϵᵢ)
    # end
    # return(mapreduce(permutedims, vcat, ϵ))
    return(ϵ)
end

function gaussian_squared_ϵ(n, size)
    # d = Normal()
    # ϵ = rand(d, size, n)

    d = Normal(1, 2)
    #ϵ = [rand(d, size) for i in 1:n]
    #ϵ =  mapreduce(permutedims, vcat, ϵ)
    ϵ = transpose(rand(d, n, size))

    for i in 1:size 
        for j in 1:n
            ϵ[i,j] = ϵ[i,j]^2 #sign(ϵ[i,j])*ϵ[i,j]^2 + 2
        end
    end
    return(ϵ)
end


function beta_ϵ(n, size)
    # d = Normal()
    # ϵ = rand(d, size, n)

    d = Beta(1,3)
    ϵ = transpose(rand(d, n, size))
    return(ϵ)
end


function samples(Λ, ϵ)
    n = size(Λ)[1]

    I = diagm(ones(n))
    A = transpose(inv(I-Λ))

    X = A*transpose(ϵ)
    return(transpose(X))
end

function tensorS_v2(X)
    n = size(X)[2]
    m = size(X)[1]
    S = zeros(n,n);
    for i in 1:n
        X₁ = X[:,i]; x₁ = mean(X₁)
        for j in i:n
            X₂ = X[:,j]; x₂ = mean(X₂)
            s = sum((X₁.-x₁).*(X₂.-x₂))
            S[i,j] = s/m
            S[j,i] = s/m
        end
    end

    return(S)
end

function tensorS(X)
    n = size(X)[2]
    m = size(X)[1]
    S = zeros(n,n);
    for i in 1:n
        for j in 1:n
            X₁ = X[:,i]; x₁ = mean(X₁)
            X₂ = X[:,j]; x₂ = mean(X₂)
            s = 0
            for l in 1:m
                s = s + (X₁[l]-x₁)*(X₂[l]-x₂)
            end
            S[i,j] = s/m
        end
    end

    return(S)
end

function tensorT(X)
    n = size(X)[2]
    m = size(X)[1]
    T = zeros(n,n,n);
    for i in 1:n
        Xᵢ = X[:,i]; mᵢ = mean(Xᵢ);
        for j in i:n
            Xⱼ = X[:,j]; mⱼ = mean(Xⱼ)
            for k in j:n
                Xₖ = X[:,k]; mₖ = mean(Xₖ)
                #print(i,j,k,"\n")
                #s = sum((Xᵢ .- mᵢ).*(Xⱼ .- mⱼ).*(Xₖ .- mₖ))
                s = 0
                for l in 1:m
                    s = s + (X[l,i] - mᵢ)*(X[l,j] - mⱼ)*(X[l,k] - mₖ)
                end
                T[i,j,k] = s/m
                T[i,k,j] = s/m
                T[j,i,k] = s/m
                T[j,k,i] = s/m
                T[k,i,j] = s/m
                T[k,j,i] = s/m
            end
        end
    end

    return(T)
end

function moments_ϵ(ϵ, ord)
    n = size(ϵ)[2]
    m = size(ϵ)[1]
    Ω = zeros(n);
    for i in 1:n  
        ϵᵢ = ϵ[:,i]; xᵢ = mean(ϵᵢ)
        s = 0
        for l in 1:m
            s = s + (ϵᵢ[l]-xᵢ)^ord
        end
        Ω[i] = s/m
    end

    return(Ω)
end
