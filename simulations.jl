include("Algorithm/algorithm.jl")
include("Algorithm/samples.jl")
using IterTools
using Graphs
using GraphPlot
using Gadfly
using Plots
using DataFrames
using CSV

Λ₁ = [0 1 0 0 0
      0 0 1 0 1
      0 0 0 1 0
      0 1 0 0 0
      0 0 0 0 0
]
Λ₂ = [0 1 0 0 0 1 0 0
      0 0 1 0 0 0 0 0
      0 0 0 1 0 0 0 0
      1 0 0 0 0 0 0 1
      0 0 0 0 0 1 0 0
      0 0 0 0 0 0 1 0
      0 0 0 0 1 0 0 1
      0 0 0 0 0 0 0 0
]
Λ₃ = [0 1 0 0 0 0 1 0 0 0
      0 0 1 0 0 0 0 1 1 0
      0 0 0 1 0 0 0 0 0 0
      0 0 0 0 1 1 0 0 0 0
      1 0 0 0 0 0 0 0 1 0
      0 0 0 0 0 0 1 0 0 0
      0 0 0 0 0 0 0 1 0 0
      0 0 0 0 0 1 0 0 1 0
      0 0 0 0 0 0 0 0 0 0
      0 0 0 0 0 0 1 0 0 0
]

G1 = SimpleDiGraph(Λ₁); #gplot(G1,nodelabel=1:size(Λ₁)[1])
G2 = SimpleDiGraph(Λ₂); #gplot(G2,nodelabel=1:size(Λ₂)[1])
G3 = SimpleDiGraph(Λ₃); #gplot(G3,nodelabel=1:size(Λ₃)[1])

Gs = [G1, G2, G3]

sampleSizes = [100, 1000, 10000, 50000, 100000, 500000, 1000000, 5000000]

epsilons_model = [1e-9, 5e-9, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6]
#epsilons_sample = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]
epsilons_sample = [0.00001:0.0004:0.01;]

for g in 1:3

    g₀ = Gs[g]
    n = length(collect(vertices(g₀)))
    pg = gplot(g₀,nodelabel=1:n)
    #draw(PNG("plots/graph_"*string(g)*".png", 9inch, 6inch), pg)
    
    Λ =  collect(adjacency_matrix(g₀))

    Vertices = vcat([1:n]...)
    Edges = [[edge.src, edge.dst] for edge in collect(edges(g₀))]
    G = Graph(Edges, Vertices);

    P = length(Edges)
    N = n^2-n - P

    ########################
    df_model = DataFrame(it=Int64[], lambda=Matrix[], eps=Float64[], sample_size=Int64[], 
                P=Float64[], N=Float64[], TP=Float64[], FP=Float64[], TN=Int64[], FN=Int64[],
                accuracy=Float64[], F1=Float64[], time = Float64[], matrix_error=Float64[], matrix_error_TP=Union{Missing, Float64}[])
    df_sample = DataFrame(it=Int64[], lambda=Matrix[], eps=Float64[], sample_size=Int64[], 
                P=Float64[], N=Float64[], TP=Float64[], FP=Float64[], TN=Int64[], FN=Int64[],
                accuracy=Float64[], F1=Float64[], time = Float64[], matrix_error=Float64[], matrix_error_TP=Union{Missing, Float64}[])
    ########################
   
    for ii in 1:100
        Λ₀ = Λ.*rand(Uniform(0.1,1),n,n).*(-1).^rand(1:10,n,n)
        
        #Λ₀ = (Λ₀)./sum(Λ₀)*10
        print("\n", ii, " ")
        for ss in sampleSizes
            print(ss, " ")
            # println("new ss: ", ss,  "\n")
            ϵ = beta_ϵ(n, ss);
            X = samples(Λ₀, ϵ);

            Ω₂ = diagm(moments_ϵ(ϵ, 2));
            Ω₃ = diagm(moments_ϵ(ϵ, 3));

            S = cov(X);
            T = tensorT(X);

            moments_sample = Moments(Vertices, S, T, indices(Vertices));
            moments_model = Moments(G, Λ₀, Ω₂, Ω₃);

            # for e in 1:length(epsilons_model)
            #     println(ii, " ", ss," ", e)
            #     ##### MODEL #####
            #     # println("model \n")
            #     dgts = 7
            #     top_order = []; Λ_m = zeros(n,n);
            #     vertices_ind = moments_model.indices;

            #     (top_m, Λ_m) = main_algorithm(moments_model, top_order, Λ_m, vertices_ind, "cond_number", epsilons_model[e]);

            #     Λ_m = round.(Λ_m, digits=dgts)
            #     edges_model = [[edge.src, edge.dst] for edge in collect(edges(SimpleDiGraph(Λ_m)))]

            #     TP_model = length(intersect(Edges, edges_model))
            #     FP_model = length(setdiff(edges_model, Edges))
            #     TN_model = N - FP_model
            #     #TN2_model = n^2-n - length(edges_model) - FN_model
            #     FN_model = length(setdiff(Edges, edges_model))
            #     accuracy_model = (TP_model + TN_model)/(TP_model+FP_model+TN_model+FN_model)
            #     F1_model = 2*TP_model/(2*TP_model+FP_model+FN_model)

            #     data_model = [ii, Λ_m, epsilons_model[e], ss, P, N, TP_model, FP_model, TN_model, FN_model, accuracy_model, F1_model]
            #     push!(df_model, data_model)
            # end

            for e in 1:length(epsilons_sample)
                # println(ii, " ", ss," ", e, "\n")
               
                #code
                ##### SAMPLES #####
                #println("sample \n")
                dgts = 2
                top_order = []; Λ_s = zeros(n,n);
                vertices_ind = moments_sample.indices;
                t = @elapsed (top_s, Λ_s) = main_algorithm(moments_sample, top_order, Λ_s, vertices_ind, "cond_number", epsilons_sample[e]);
                Λ_s = round.(Λ_s, digits=dgts)
                err = sqrt(sum(Λ_s .- Λ₀).^2)
               

                edges_sample = [[edge.src, edge.dst] for edge in collect(edges(SimpleDiGraph(Λ_s)))]
                correctEdges = intersect(Edges, edges_sample)
                if length(correctEdges) > 0
                    err_TP = sqrt(sum([(Λ_s[edge[1], edge[2]]- Λ₀[edge[1], edge[2]])^2 for edge in correctEdges]))
                else
                    err_TP = missing
                end

                TP_sample = length(intersect(Edges, edges_sample))
                FP_sample = length(setdiff(edges_sample, Edges))
                TN_sample = N - FP_sample
                FN_sample = length(setdiff(Edges, edges_sample))
                accuracy_sample = (TP_sample + TN_sample)/(TP_sample+FP_sample+TN_sample+FN_sample)
                F1_sample = 2*TP_sample/(2*TP_sample+FP_sample+FN_sample)

                data_sample = [ii, Λ_m, epsilons_sample[e], ss, P, N, TP_sample, FP_sample, TN_sample, FN_sample, accuracy_sample, F1_sample, t, err, err_TP]
                push!(df_sample, data_sample);
            end
        end
    end

   ### df
   transform!(df_sample, [:TP,:P] => (./) => :TPR);
   transform!(df_sample, [:FP,:N] => (./) => :FPR);

   transform!(df_sample, [:TP,:TN] => (.+) => :TPTN);
   transform!(df_sample, [:FP,:FN] => (.+) => :FPFN);

   #CSV.write("data/df_graph_"*string(g)*".csv", df);

   ### Plots
    df_agg_s = combine(DataFrames.groupby(df_sample, [:eps, :sample_size]), [:accuracy, :F1] .=> maximum, [:accuracy, :F1] .=> mean);
 
    p1 = Gadfly.plot(df_agg_s, x=:sample_size, y=:accuracy_maximum, color=:eps, Geom.line, Scale.x_discrete());
    p2 = Gadfly.plot(df_agg_s, x=:eps, y=:accuracy_maximum, color=:sample_size, Geom.line, Scale.color_discrete());
    p3 = Gadfly.plot(df_agg_s, x=:sample_size, y=:F1_maximum, color=:eps, Geom.line, Scale.x_discrete());
    p4 = Gadfly.plot(df_agg_s, x=:eps, y=:F1_maximum, color=:sample_size, Geom.line, Scale.color_discrete());
    p_max = gridstack([p1 p2; p3 p4]);
    draw(PNG("plotsTest/results_max_graph"*string(g)*".png", 11inch, 6inch), p_max);

    p1 = Gadfly.plot(df_agg_s, x=:sample_size, y=:accuracy_mean, color=:eps, Geom.line, Scale.x_discrete());
    p2 = Gadfly.plot(df_agg_s, x=:eps, y=:accuracy_mean, color=:sample_size, Geom.line, Scale.color_discrete());
    p3 = Gadfly.plot(df_agg_s, x=:sample_size, y=:F1_mean, color=:eps, Geom.line, Scale.x_discrete());
    p4 = Gadfly.plot(df_agg_s, x=:eps, y=:F1_mean, color=:sample_size, Geom.line, Scale.color_discrete());
    p_mean = gridstack([p1 p2; p3 p4]);
    draw(PNG("plotsTest/results_mean_graph"*string(g)*".png", 11inch, 6inch), p_mean);


    df_agg_s = combine(DataFrames.groupby(df_sample, [:eps, :sample_size]), [:TPR, :FPR] .=> maximum, [:TPR, :FPR, :TPTN, :FPFN] .=> mean);
    df_ROC = filter(:sample_size => ==(50000), df_agg_s);
    df_ROC = combine(DataFrames.groupby(df_ROC, :FPR_mean), :TPR_mean => mean);

    # # df_ROC = filter(:eps => ==(0), df_sample);
    # #ROC =  Gadfly.plot(df_ROC, x=:FPR, y=:TPR, Geom.line);
    # ROC =  Gadfly.plot(df_ROC, x=:FPR_mean, y=:TPR_mean_mean, Geom.line, Coord.cartesian(xmin=0, xmax=1, ymin=0, ymax=1));
    # draw(PNG("plotsTest/ROC_graph_"*string(g)*".png", 8inch, 8inch), ROC);


    # #df_agg_s = combine(DataFrames.groupby(df_sample, [:eps, :sample_size]), [:TPTN, :FPFN] .=> mean);
    # df_ROC2 = filter(:sample_size => ==(50000), df_agg_s);
    # ROC2 =  Gadfly.plot(df_ROC2, x=:FPFN_mean, y=:TPTN_mean, Geom.line);
    # draw(PNG("plotsTest/ROC2_graph_"*string(g)*".png", 8inch, 8inch), ROC2);


    p1 = Gadfly.plot(df_agg_s, x=:eps, y=:TPR_maximum, color=:sample_size, Geom.line, Scale.color_discrete());
    p2 = Gadfly.plot(df_agg_s_mean, x=:eps, y=:TPR_mean, color=:sample_size, Geom.line, Scale.color_discrete());
    #p3 = Gadfly.plot(df_agg_s, x=:sample_size, y=:rP_maximum, color=:eps, Geom.line, Scale.x_discrete())
    p = gridstack([p1 p2]);
    draw(PNG("plotsTest/trueEdges_graph_"*string(g)*".png", 11inch, 6inch), p);
    
    df_time = combine(DataFrames.groupby(df_sample, [:sample_size]), [:time] .=> mean)
    println(df_time)
    CSV.write("data/df_time_"*string(g)*".csv", df_time);

    # df_error = combine(DataFrames.groupby(df_sample, [:eps, :sample_size]), :matrix_error => mean, :matrix_error_TP => matrix_error_TP -> mean(skipmissing(matrix_error_TP)))
    # df_error = filter(:matrix_error_mean => <=(sort(df_error[:, :matrix_error_mean], rev=true)[7]), df_error)
    # df_error = filter(:matrix_error_TP_function => <=(sort(df_error[:, :matrix_error_TP_function], rev=true)[7]), df_error)
    # p1 = Gadfly.plot(df_error, x=:eps, y=:matrix_error_mean, color=:sample_size, Geom.line, Scale.color_discrete());
    # p2 = Gadfly.plot(df_error, x=:eps, y=:matrix_error_TP_function, color=:sample_size, Geom.line, Scale.color_discrete());
    # p = gridstack([p1 p2])
    # draw(PNG("plotsTest/matrices_graph_"*string(g)*".png", 11inch, 6inch), p);
   
end
