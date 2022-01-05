using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using HypothesisTests

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates2traits.csv", DataFrame, missingstring="N/A")
########################################

########################################
Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
########################################

########################################
gmax            = log(2)./Min_gen_time
V_p             = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_r             = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
k_E             = DEBmicroTrait.translation_power(V_p, V_r, Min_gen_time)
########################################

########################################
y_EV            = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
########################################

########################################
# I/O
df_out              = DataFrame()
df_out.k_E          = k_E
df_out.y_EV         = y_EV
df_out.response     = df_isolates.Rhizosphere_response

save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_protein_synthesis.jld", "kE", k_E, "yEV", y_EV, "mingt", Min_gen_time)
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_protein_synthesis.csv", df_out)

########################################

########################################
# statistics
df_p       = filter(x->(x.response.=="positive"), df_out)
df_n       = filter(x->(x.response.=="negative"), df_out)
df_u       = filter(x->(x.response.=="undefined"), df_out)

k_E_kw     = KruskalWallisTest(df_p.k_E, df_n.k_E)
k_E_p      = median(df_p.k_E)
k_E_n      = median(df_n.k_E)

y_EV_kw     = KruskalWallisTest(df_p.y_EV, df_n.y_EV)
y_EV_p      = median(df_p.y_EV)
y_EV_n      = median(df_n.y_EV)
########################################
