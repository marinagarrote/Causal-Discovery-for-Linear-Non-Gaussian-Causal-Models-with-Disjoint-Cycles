using LinearAlgebra, Graphs

include("moments.jl")
include("regression.jl")
include("determinants.jl")
include("root_cycles.jl")
include("root_nodes.jl")
include("parameters.jl")

# eps = 1e-7
# dgts = 9

function topological_order(moments, top_order, Λ, vertices_ind, method_det, eps, DAG=false)
    if(length(moments.V) < 1) 
        # source_edges
        return(top_order, Λ)
    end
    
    top_order_depth = length(top_order)
    (moments, top_order, Λ) = root_nodes(moments, top_order, Λ, method_det, eps);
    if !DAG
        (moments, top_order, Λ, vertices_ind) = root_cycles(moments, top_order, Λ, vertices_ind, method_det, eps);
    end

    if(top_order_depth ==  length(top_order)) 
        top_order = push!(top_order, moments.V)
        return(top_order, Λ)
    end

    return(topological_order(moments, top_order, Λ, vertices_ind, method_det, eps))
end

function main_algorithm(moments, top_order, Λ, vertices_ind, method_det, eps; DAG=false)
    (top_order, Λ) = topological_order(moments, top_order, Λ, vertices_ind, method_det, eps, DAG)
    Λ = update_Λ_sources(moments, top_order, Λ, 16)
    return(top_order, Λ)
end


function Λ_from_components(moments, top_order, Λ, digits)
    
    mom_original = moments
    for component in top_order
        for elem in component
            if length(elem) > 1
                Λ = update_Λ_cycle(moments, Λ, elem, vertices_ind)
            end
        end

        c = unique(reduce(vcat, component))
        moments = residuals(moments, c, setdiff(moments.V, c))
    end

    Λ = update_Λ_sources(mom_original, top_order, Λ, digits)
    return Λ
end