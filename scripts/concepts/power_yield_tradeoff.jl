using DEBmicroTrait
using Roots
using CSV, DataFrames
using GLM
using Statistics
using JLD

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudation_properties.csv", DataFrame, missingstring="N/A")


Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time    = df_isolates.Min_gen_time
gmax            = log(2)./Min_gen_time
Gram_stain      = convert(Array{String}, df_isolates.gram_stain)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
y_EM            = ones(size(V_cell,1))
α               = (df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)./maximum(df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)
z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)
ρ_ps            = zeros(size(V_cell,1))
y_DEs           = zeros(size(V_cell,1), 83)
FCR             = zeros(size(V_cell,1), 83)
Yrmax           = zeros(size(V_cell,1), 83)
yield           = zeros(size(V_cell,1), 100000)
rate            = zeros(size(V_cell,1), 100000)
N_C             = DEBmicroTrait.extract_composition(df_metabolites.Formula[29])[1]


for j in 1:83
    for i in 1:size(V_cell,1)
        find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
        ρ_p = Roots.find_zero(find_ρ, 1.0)
        closure = genome_distr[:,i]./sum(genome_distr[:,i])
        if df_metabolites.Ontology[j] == "Sugars"
            ρ_ps[i] = ρ_p[1].*closure[1]
        elseif df_metabolites.Ontology[j] == "Organic acids"
            ρ_ps[i] = ρ_p[1].*closure[2]
        elseif df_metabolites.Ontology[j] == "Amino acids"
            ρ_ps[i] = ρ_p[1].*closure[3]
        elseif df_metabolites.Ontology[j] == "Fatty acids"
            ρ_ps[i] = ρ_p[1].*closure[4]
        elseif df_metabolites.Ontology[j] == "Nucleotides"
            ρ_ps[i] = ρ_p[1].*closure[5]
        else
            ρ_ps[i] = ρ_p[1].*closure[6]
        end
        y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
        y_DEs[i,j] = y_DE[1]
        N_C             = DEBmicroTrait.extract_composition(df_metabolites.Formula[j])[1]
        yield_tmp, rate_tmp = DEBmicroTrait.rate_yield_trade_off(ρ_ps[i]*ones(1,1), [α[i]],[V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE, N_C, [y_EM[i]])
        for h in 1:size(rate_tmp,1)
            yield[i,h] = yield_tmp[h]
            rate[i,h] = rate_tmp[h]
        end

        rmax, rid = findmax(rate[i,:])
        Ymax, Yid = findmax(yield[i,:])
        FCR[i,j] = Ymax - yield[i,rid]
        Yrmax[i,j] = yield[i,rid]
    end
end


FCR = zeros(size(V_cell,1))
FCR_50 = zeros(size(V_cell,1))
for i in 1:39
    rmax, rid = findmax(rate[i,:])
    Ymax, Yid = findmax(yield[i,:])
    FCR[i] = Ymax - yield[i,rid]
    FCR_50[i] = yield[i,rid]
end

Ymaxmed  = 78.5
median(FCR_50*100)



df_out = DataFrame()
df_out.Vcell = V_cell
df_out.FCR = FCR
df_out.response = df_isolates.Rhizosphere_response

linearRegressor = lm(@formula(log(FCR) ~ log(Vcell)), df_out)
r2(linearRegressor)
summary(linearRegressor)


CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_power_yield.csv", df_out)

# Bionumber  	116954
L_DNA = collect(range(0.1, 13.5, length=100)).*1e6
V_cell = DEBmicroTrait.genome_size_to_cell_volume(L_DNA)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(L_DNA)
min_gentime = DEBmicroTrait.gmax_regression(rrn_copies)
gmax = log(2)./min_gentime
Gram_stain = repeat(["+"], 1000)

V_P = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_R = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
k_E = DEBmicroTrait.translation_power(V_P, V_R, min_gentime)
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
k_M = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell, min_gentime, Gram_stain)

minimum(k_E)
maximum(k_E)

minimum(y_EV)
maximum(y_EV)

minimum(k_M)
maximum(k_M)
median(k_M)

save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/power_yield_tradeoff.jld", "gmax", gmax, "kE", k_E, "yEV", y_EV, "kM", k_M, "rrn", rrn_copies, "Ldna", L_DNA)
