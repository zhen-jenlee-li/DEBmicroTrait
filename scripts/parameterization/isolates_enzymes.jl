using DEBmicroTrait
using CSV, DataFrames, Statistics
using HypothesisTests
using JLD

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")
########################################

########################################
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
#α               = DEBmicroTrait.constrain_enzyme_allocation(V_cell, Min_gen_time, Gram_stain, df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)
zh              =  df_isolates.z_hydrolases./df_isolates.Genome_size*1e6
α_X             =  1e-2*(df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)./maximum(df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)
########################################

########################################
# statistics
df_isolates_out = DataFrame()
df_isolates_out.response = df_isolates.Rhizosphere_response
df_isolates_out.zh = df_isolates.z_hydrolases./df_isolates.Genome_size*1e6
df_isolates_out.α_X = α_X

df_p_α = filter(x->(x.response.=="positive") , df_isolates_out)
df_n_α = filter(x->(x.response.=="negative") , df_isolates_out)
df_u_α = filter(x->(x.response.=="undefined") , df_isolates_out)
zh_kw   = KruskalWallisTest(df_p_α.zh, df_n_α.zh)
a_kw   = KruskalWallisTest(df_p_α.α_X, df_n_α.α_X)
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_enzymes.jld", "zh", zh, "alpha", α_X)
########################################
