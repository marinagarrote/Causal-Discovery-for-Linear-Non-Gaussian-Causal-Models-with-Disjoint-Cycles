include("Algorithm/algorithm.jl")
include("Algorithm/samples.jl")
using IterTools
using Graphs
using GraphPlot
using Gadfly
using Plots
using DataFrames
using CSV

function create_graph(n)
    found = false
    Λ = zeros(n,n)
    while !found
        Λ = zeros(n,n)
        n1 = rand((3:5))
        n2 = rand((3:(n-n1)))
        g = SimpleDiGraph(n) 
    
        for i in 1:(n1-1)
            print(i," ", i+1,"\n")
            add_edge!(g, i, i+1);
        end
        print(n1," ", 1,"\n")
        add_edge!(g, n1, 1);
    
        for i in 1:(n2-1)
            print(i+n1," ", i+n1+1,"\n")
            add_edge!(g, i+n1, i+n1+1);
        end
        add_edge!(g, n2+n1, 1+n1);
    
        edgs = collect(edges(g))
    
    
       
        for edge in edgs
            Λ[edge.src, edge.dst] = 1
        end
    
        singlenodes = vcat([1+n1+n2:n]...)
    
        idx_c1 = vec(collect(product(1:n1, n1+1:n)))
        idx_c2 = vec(collect(product(setdiff(1:n, n1+1:n1+n2), n1+1:n1+n2)))
        idx_singnodes = vec(collect(product(1:n1+n2, singlenodes)))
    
        idx = unique(vcat(idx_c1, idx_c2, idx_singnodes))
        rand_idx = rand(1:length(idx), rand(3:10, 1)[1])
    
        for i in vcat(rand_idx)
            Λ[idx[i][1], idx[i][2]] = 1
        end
    
    
        for i in 1:n
            for j in 1:n
                if i == j
                    Λ[i,j] = 0
                elseif Λ[i,j] != 0 && Λ[j,i] != 0
                    Λ[j,i] = 0
                end
            end
        end
    
    
        g₀ = SimpleDiGraph(Λ)
        cycles = vcat(simplecycles(g₀)...)
        un_cycles = unique(cycles)
    
        if length(cycles) == length(un_cycles)
            found = true
        end
    end
    return(Λ)
end


#number of vertices
n = 10
Λ = create_graph(n)