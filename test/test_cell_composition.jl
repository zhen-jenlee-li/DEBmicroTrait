using DEBmicroTrait
using CSV, DataFrames
using HypothesisTests, Statistics

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudation_properties.csv", DataFrame, missingstring="N/A")

Min_gen_time    = df_isolates.Min_gen_time
gmax = log(2)./Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)


df_isolates.SAV = DEBmicroTrait.sav_from_growth_rate(gmax)
df_isolates.SAV = DEBmicroTrait.surface_area_volume_ratio(V_cell)
kcat = 180*60^2
df_isolates.sigma_p = @. gmax/(df_isolates.SAV*1e-12*kcat)

z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

ρ_ps            = zeros(size(df_metabolites.Name,1), size(V_cell,1))


for j in 1:size(df_metabolites.Name,1)
    if df_metabolites.Ontology[j] == "Sugars"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p = df_isolates.sigma_p[i]
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[1]
        end
    elseif df_metabolites.Ontology[j] == "Organic acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = df_isolates.sigma_p[i]
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[2]
        end
    elseif df_metabolites.Ontology[j] == "Amino acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = df_isolates.sigma_p[i]
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[3]
        end
    elseif df_metabolites.Ontology[j] == "Fatty acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = df_isolates.sigma_p[i]
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[4]
        end
    elseif df_metabolites.Ontology[j] == "Nucleotides"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = df_isolates.sigma_p[i]
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[5]
        end
    else
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = df_isolates.sigma_p[i]
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[6]
        end
    end
end


ρ_ps[ρ_ps.==0.0] .= 1e-12

N_C = zeros(size(df_metabolites.Name,1))
for i in 1:size(df_metabolites.Name,1)
    elementstring = convert(String, df_metabolites.Formula[i])
    N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
end

N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C

D_S               = DEBmicroTrait.aqueous_diffusivity(df_metabolites.Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, D_S)

a_s               = Vmax./K_D

df_assimilation = DataFrame()
df_assimilation.transporter_density = vec(ρ_ps)
df_assimilation.Vmax = vec(Vmax)
df_assimilation.KD = vec(K_D)
df_assimilation.affinity = vec(Vmax./K_D)

df_assimilation.ontology = repeat(df_metabolites.Ontology, size(V_cell,1))
response = Array{String}(undef, (size(df_metabolites.Name,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     response[:,i] .= df_isolates.Rhizosphere_response[i]
 end

df_assimilation.response = vec(response)
genome_size = Array{Int}(undef, (size(df_metabolites.Name,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     genome_size[:,i] .= df_isolates.Genome_size[i]
end
df_assimilation.genome_size = vec(genome_size)

cell_size = Array{Float64}(undef, (size(df_metabolites.Name,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     cell_size[:,i] .= V_cell[i]
end
df_assimilation.cell_size = vec(cell_size)
min_gt = Array{Float64}(undef, (size(df_metabolites.Name,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     min_gt[:,i] .= df_isolates.Min_gen_time[i]
end
df_assimilation.mingt = vec(min_gt)

SAV = Array{Float64}(undef, (size(df_metabolites.Name,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     SAV[:,i] .= df_isolates.SAV[i]
end
df_assimilation.sav = vec(SAV)

df_assimilation.yield = vec(abs.(y_DEs))

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_assimilation_sav.csv", df_assimilation)
