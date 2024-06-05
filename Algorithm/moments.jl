using LinearAlgebra

#eps = 1e-1
#digits = 1

#-------------
#  Objects
#-------------
mutable struct Graph
    E::Vector{Vector{Int64}} 
    V::Vector{Int64}
    indices::Dict{Int64, Int64}
end

struct Moments
    V::Vector{Int64} 
    S::Matrix{Float64}
    T::Array{Float64, 3}
    indices::Dict{Int64, Int64}
end

function Graph(Edges, Vertices)
    Graph(Edges, Vertices, indices(Vertices))
end

function Moments(G, Λ, Ω₂, Ω₃)
    S = momentS(Λ, Ω₂)
    T = momentT(Λ, Ω₃)

    return(Moments(G.V, S, T, indices(G.V)))
end

#-----------------------
# Auxiliary functions
#-----------------------

function remove!(a, item)
    deleteat!(a, findall(x -> x == item,a))
end

function indices(V)
    return(Dict(V .=> collect(1:length(V))))
end

function getIndices(moment::Moments, C)
    return(getindex.(Ref(moment.indices), C))
end

function getIndices(indices::Dict{Int64, Int64}, C)
    return(getindex.(Ref(indices), C))
end

#-------------
#  Data
#-------------

function randomΛ(Edges, Vertices)
    ind = indices(Vertices)
    n = length(Vertices)
    Λ = zeros(n,n)
    for e in Edges
        Λ[ind[e[1]],ind[e[2]]] = rand()+rand()-rand()
    end
    return(Λ)
end

function randomΩ(n)
    return(diagm(rand(n)*10))
end

function momentS(Λ, Ω₂)
    n = size(Λ)[1]
    I = diagm(ones(n))
    A = inv(I-Λ)
    S = transpose(A)*Ω₂*A
end

function momentT(Λ, Ω₃)
    n = size(Λ)[1]
    I = diagm(ones(n))
    A = inv(I-Λ)

    T = zeros(n,n,n)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                T[i,j,k] = 0
                for s in 1:n
                    T[i,j,k] = T[i,j,k] + Ω₃[s,s]*A[s,i]*A[s,j]*A[s,k]
                end
            end
        end
    end

    return(T)
end

function randomParameters(G)
    n = length(G.V)
    Λ = randomΛ(G.E, G.V);
    Ω₂ = randomΩ(n);
    Ω₃ = randomΩ(n);
    return(Λ, Ω₂, Ω₃)
end

function Moments(G, Λ, Ω₂, Ω₃)
    S = momentS(Λ, Ω₂)
    T = momentT(Λ, Ω₃)

    return(Moments(G.V, S, T, indices(G.V)))
end
