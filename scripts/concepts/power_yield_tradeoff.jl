using DEBmicroTrait
using Roots
using CSV, DataFrames
using GLM

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files/exudation_properties.csv", DataFrame, missingstring="N/A")


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
y_DEs           = zeros(size(V_cell,1))
FCR             = zeros(size(V_cell,1))
yield           = zeros(size(V_cell,1), 100000)
rate            = zeros(size(V_cell,1), 100000)

for i in 1:size(V_cell,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[29])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    closure = genome_distr[:,i]./sum(genome_distr[:,i])
    ρ_ps[i] = ρ_p[1].*closure[1]
    y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[29])
    y_DEs[i] = y_DE[1]
    yield_tmp, rate_tmp = DEBmicroTrait.rate_yield_trade_off(ρ_ps[i]*ones(1,1), [α[i]],[V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE, N_C, [y_EM[i]])
    for h in 1:size(rate_tmp,1)
        yield[i,h] = yield_tmp[h]
        rate[i,h] = rate_tmp[h]
    end
end

FCR = zeros(size(V_cell,1))
for i in 1:39
    rmax, rid = findmax(rate[i,:])
    Ymax, Yid = findmax(yield[i,:])
    FCR[i] = Ymax - yield[i,rid]
end


df_out = DataFrame()
df_out.Vcell = V_cell
df_out.FCR = FCR
df_out.response = df_isolates.Rhizosphere_response

linearRegressor = lm(@formula(log(FCR) ~ log(Vcell)), df_out)
r2(linearRegressor)


CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_power_yield.csv", df_out)
