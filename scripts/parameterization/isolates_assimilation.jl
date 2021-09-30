using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using JLD

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/exudation_properties.csv", DataFrame, missingstring="N/A")
########################################

########################################
# isolate traits
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)

z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

y_EM            = ones(size(V_cell,1))
########################################

########################################
# calc transporter density
ρ_ps            = zeros(size(df_metabolites.Name,1), size(V_cell,1))
y_DEs           = zeros(size(df_metabolites.Name,1), size(V_cell,1))

for j in 1:size(df_metabolites.Name,1)
    if df_metabolites.Ontology[j] == "Sugars"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[1]
            y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i] = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Organic acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[2]
            y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i] = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Amino acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[3]
            y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i] = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Fatty acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[4]
            y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i] = y_DE[1]
        end
    elseif df_metabolites.Ontology[j] == "Nucleotides"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[5]
            y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i] = y_DE[1]
        end
    else
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]],[rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[6]
            y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_metabolites.Formula[j])
            y_DEs[j,i] = y_DE[1]
        end
    end
end

ρ_ps[ρ_ps.==0.0] .= 1e-12

N_C = zeros(size(df_metabolites.Name,1))
for i in 1:size(df_metabolites.Name,1)
    elementstring = df_metabolites.Formula[i]
    N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
end

N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C

D_S               = DEBmicroTrait.aqueous_diffusivity(df_metabolites.Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, D_S)

a_s               = Vmax./K_D
# I/O
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
     SAV[:,i] .= DEBmicroTrait.surface_area_volume_ratio([V_cell[i]])
end
df_assimilation.sav = vec(SAV)

df_assimilation.yield = vec(abs.(y_DEs))

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_assimilation.csv", df_assimilation)
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_assimilation.jld", "rho", ρ_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DEs, "NC", N_C)


# statistics
using HypothesisTests
#sugars
df_p_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="positive") , df_assimilation)
df_n_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="negative") , df_assimilation)
df_u_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="undefined") , df_assimilation)
ρ_sugars    = KruskalWallisTest(df_p_sugars.transporter_density, df_n_sugars.transporter_density)
ρ_sugars_p  = median(df_p_sugars.transporter_density)
ρ_sugars_n  = median(df_n_sugars.transporter_density)
ρ_sugars_u  = median(df_u_sugars.transporter_density)
Vmax_sugars    = KruskalWallisTest(df_p_sugars.Vmax, df_n_sugars.Vmax)
Vmax_sugars_p  = median(df_p_sugars.Vmax)
Vmax_sugars_n  = median(df_n_sugars.Vmax)
Vmax_sugars_u  = median(df_u_sugars.Vmax)
K_sugars    = KruskalWallisTest(df_p_sugars.KD, df_n_sugars.KD)
K_sugars_p  = median(df_p_sugars.KD)
K_sugars_n  = median(df_n_sugars.KD)
K_sugars_u  = median(df_u_sugars.KD)
a_sugars    = KruskalWallisTest(df_p_sugars.Vmax./df_p_sugars.KD, df_n_sugars.Vmax./df_n_sugars.KD)
a_sugars_p  = median(df_p_sugars.Vmax./df_p_sugars.KD)
a_sugars_n  = median(df_n_sugars.Vmax./df_n_sugars.KD)

#organic acids
df_p_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="positive") , df_assimilation)
df_n_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="negative") , df_assimilation)
df_u_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="undefined") , df_assimilation)
ρ_organics    = KruskalWallisTest(df_p_organics.transporter_density, df_n_organics.transporter_density)
ρ_organics_p  = median(df_p_organics.transporter_density)
ρ_organics_n  = median(df_n_organics.transporter_density)
ρ_organics_u  = median(df_u_organics.transporter_density)
Vmax_organics    = KruskalWallisTest(df_p_organics.Vmax, df_n_organics.Vmax)
Vmax_organics_p  = median(df_p_organics.Vmax)
Vmax_organics_n  = median(df_n_organics.Vmax)
Vmax_organics_u  = median(df_u_organics.Vmax)
K_organics    = KruskalWallisTest(df_p_organics.KD, df_n_organics.KD)
K_organics_p  = median(df_p_organics.KD)
K_organics_n  = median(df_n_organics.KD)
K_organics_u  = median(df_u_organics.KD)
a_organics    = KruskalWallisTest(df_p_organics.Vmax./df_p_organics.KD, df_n_organics.Vmax./df_n_organics.KD)
a_organics_p  = median(df_p_organics.Vmax./df_p_organics.KD)
a_organics_n  = median(df_n_organics.Vmax./df_n_organics.KD)

#amino acids
df_p_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="positive") , df_assimilation)
df_n_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="negative") , df_assimilation)
df_u_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="undefined") , df_assimilation)
ρ_aminos    = KruskalWallisTest(df_p_aminos.transporter_density, df_n_aminos.transporter_density,  df_u_aminos.transporter_density)
ρ_aminos_p  = median(df_p_aminos.transporter_density)
ρ_aminos_n  = median(df_n_aminos.transporter_density)
ρ_aminos_u  = median(df_u_aminos.transporter_density)
Vmax_aminos    = KruskalWallisTest(df_p_aminos.Vmax, df_n_aminos.Vmax)
Vmax_aminos_p  = median(df_p_aminos.Vmax)
Vmax_aminos_n  = median(df_n_aminos.Vmax)
Vmax_aminos_u  = median(df_u_aminos.Vmax)
K_aminos    = KruskalWallisTest(df_p_aminos.KD, df_n_aminos.KD)
K_aminos_p  = median(df_p_aminos.KD)
K_aminos_n  = median(df_n_aminos.KD)
K_aminos_u  = median(df_u_aminos.KD)
a_aminos    = KruskalWallisTest(df_p_aminos.Vmax./df_p_aminos.KD, df_n_aminos.Vmax./df_n_aminos.KD)
a_aminos_p  = median(df_p_aminos.Vmax./df_p_aminos.KD)
a_aminos_n  = median(df_n_aminos.Vmax./df_n_aminos.KD)

#fatty acids
df_p_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="positive") , df_assimilation)
df_n_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="negative") , df_assimilation)
df_u_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="undefined") , df_assimilation)
ρ_fattys    = KruskalWallisTest(df_p_fattys.transporter_density, df_n_fattys.transporter_density, df_u_fattys.transporter_density)
ρ_fattys_p  = median(df_p_fattys.transporter_density)
ρ_fattys_n  = median(df_n_fattys.transporter_density)
ρ_fattys_u  = median(df_u_fattys.transporter_density)
Vmax_fattys    = KruskalWallisTest(df_p_fattys.Vmax, df_n_fattys.Vmax)
Vmax_fattys_p  = median(df_p_fattys.Vmax)
Vmax_fattys_n  = median(df_n_fattys.Vmax)
Vmax_fattys_u  = median(df_u_fattys.Vmax)
K_fattys    = KruskalWallisTest(df_p_fattys.KD, df_n_fattys.KD)
K_fattys_p  = median(df_p_fattys.KD)
K_fattys_n  = median(df_n_fattys.KD)
K_fattys_u  = median(df_u_fattys.KD)
a_fattys    = KruskalWallisTest(df_p_fattys.Vmax./df_p_fattys.KD, df_n_fattys.Vmax./df_n_fattys.KD)
a_fattys_p  = median(df_p_fattys.Vmax./df_p_fattys.KD)
a_fattys_n  = median(df_n_fattys.Vmax./df_n_fattys.KD)

# Nucleotides
df_p_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="positive") , df_assimilation)
df_n_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="negative") , df_assimilation)
df_u_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="undefined") , df_assimilation)
ρ_nucleos    = KruskalWallisTest(df_p_nucleos.transporter_density, df_n_nucleos.transporter_density, df_u_nucleos.transporter_density)
ρ_nucleos_p  = median(df_p_nucleos.transporter_density)
ρ_nucleos_n  = median(df_n_nucleos.transporter_density)
ρ_nucleos_u  = median(df_u_nucleos.transporter_density)
Vmax_nucleos    = KruskalWallisTest(df_p_nucleos.Vmax, df_n_nucleos.Vmax)
Vmax_nucleos_p  = median(df_p_nucleos.Vmax)
Vmax_nucleos_n  = median(df_n_nucleos.Vmax)
Vmax_nucleos_u  = median(df_u_nucleos.Vmax)
K_nucleos    = KruskalWallisTest(df_p_nucleos.KD, df_n_nucleos.KD)
K_nucleos_p  = median(df_p_nucleos.KD)
K_nucleos_n  = median(df_n_nucleos.KD)
K_nucleos_u  = median(df_u_nucleos.KD)
a_nucleos    = KruskalWallisTest(df_p_nucleos.Vmax./df_p_nucleos.KD, df_n_nucleos.Vmax./df_n_nucleos.KD)
a_nucleos_p  = median(df_p_nucleos.Vmax./df_p_nucleos.KD)
a_nucleos_n  = median(df_n_nucleos.Vmax./df_n_nucleos.KD)

#Auxins
df_p_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="positive") , df_assimilation)
df_n_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="negative") , df_assimilation)
df_u_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="undefined") , df_assimilation)
ρ_auxins    = KruskalWallisTest(df_p_auxins.transporter_density, df_n_auxins.transporter_density)
ρ_auxins_p  = median(df_p_auxins.transporter_density)
ρ_auxins_n  = median(df_n_auxins.transporter_density)
ρ_auxins_u  = median(df_u_auxins.transporter_density)
Vmax_auxins    = KruskalWallisTest(df_p_auxins.Vmax, df_n_auxins.Vmax)
Vmax_auxins_p  = median(df_p_auxins.Vmax)
Vmax_auxins_n  = median(df_n_auxins.Vmax)
Vmax_auxins_u  = median(df_u_auxins.Vmax)
K_auxins    = KruskalWallisTest(df_p_auxins.KD, df_n_auxins.KD)
K_auxins_p  = median(df_p_auxins.KD)
K_auxins_n  = median(df_n_auxins.KD)
K_auxins_u  = median(df_u_auxins.KD)
a_auxins    = KruskalWallisTest(df_p_auxins.Vmax./df_p_auxins.KD, df_n_auxins.Vmax./df_n_auxins.KD)
a_auxins_p  = median(df_p_auxins.Vmax./df_p_auxins.KD)
a_auxins_n  = median(df_n_auxins.Vmax./df_n_auxins.KD)

# yields
df_y_sugars = filter(x->(x.ontology.=="Sugars"), df_assimilation)
y_sugars    = @. 1/df_y_sugars.yield
y_sugar_med = median(y_sugars)
df_y_organics = filter(x->(x.ontology.=="Organic acids"), df_assimilation)
y_organics     =  @. 1/df_y_organics.yield
y_organics_med = median(y_organics)
df_y_aminos = filter(x->(x.ontology.=="Amino acids"), df_assimilation)
y_aminos     =  @. 1/df_y_aminos.yield
y_aminos_med = median(y_aminos)
df_y_fattys = filter(x->(x.ontology.=="Fatty acids"), df_assimilation)
y_fattys     =  @. 1/df_y_fattys.yield
y_fattys_med = median(y_fattys)
df_y_nucleos = filter(x->(x.ontology.=="Nucleotides"), df_assimilation)
y_nucleos     =  @. 1/df_y_nucleos.yield
y_nucloes_med = median(y_nucleos)
df_y_auxins = filter(x->(x.ontology.=="Auxins"), df_assimilation)
y_auxins     =  @. 1/df_y_auxins.yield
y_auxins_med = median(y_auxins)

kw_sugars_organics = KruskalWallisTest(y_sugars, y_organics)
kw_sugars_aminos = KruskalWallisTest(y_sugars, y_aminos)
kw_sugars_fattys = KruskalWallisTest(y_sugars, y_fattys)
kw_sugars_auxins = KruskalWallisTest(y_sugars, y_auxins)
kw_sugars_nucleos = KruskalWallisTest(y_sugars, y_nucleos)

kw_organics_aminos =  KruskalWallisTest(y_organics, y_aminos)
kw_organics_fattys =  KruskalWallisTest(y_organics, y_fattys)
kw_organics_auxins =  KruskalWallisTest(y_organics, y_auxins)
kw_organics_nucleos =  KruskalWallisTest(y_organics, y_nucleos)

kw_aminos_fattys =  KruskalWallisTest(y_aminos, y_fattys)
kw_aminos_auxins =  KruskalWallisTest(y_aminos, y_auxins)
kw_aminos_nucleos =  KruskalWallisTest(y_aminos, y_nucleos)

kw_fattys_auxins =  KruskalWallisTest(y_auxins, y_fattys)
kw_fattys_nucleos =  KruskalWallisTest(y_nucleos, y_fattys)

kw_nucleos_auxins = KruskalWallisTest(y_nucleos, y_auxins)

################################################################################
# Aromatic organic acids

df_aromatics = filter(x->(x.Name.=="nicotinic acid" || x.Name.=="shikimic acid" ||
                          x.Name.=="shikimic acid" || x.Name.=="salicylic acid" ||
                          x.Name.=="thymidine" || x.Name.=="threonic acid" ||
                          x.Name.=="indole-3-acetic acid" || x.Name.=="deoxyguanosine" ||
                          x.Name.=="cytidine" || x.Name.=="choline" ||
                          x.Name.=="3-Dehydroshikimic acid" || x.Name.=="trans-cinnamic acid" ||
                          x.Name.=="aminobutyric acid" || x.Name.=="L-pyroglutamic acid") , df_metabolites)

########################################
# isolate traits
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)

z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

y_EM            = ones(size(V_cell,1))
########################################

########################################
# calc transporter density
ρ_ps            = zeros(size(df_aromatics.Name,1), size(V_cell,1))
y_DEs           = zeros(size(df_aromatics.Name,1), size(V_cell,1))

for j in 1:size(df_aromatics.Name,1)
  if df_aromatics.Ontology[j] == "Sugars"
      for i in 1:size(V_cell,1)
          #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
          find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          ρ_p = Roots.find_zero(find_ρ, 1.0)
          closure = genome_distr[:,i]./sum(genome_distr[:,i])
          ρ_ps[j,i] = ρ_p[1].*closure[1]
          y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          y_DEs[j,i] = y_DE[1]
      end
  elseif df_aromatics.Ontology[j] == "Organic acids"
      for i in 1:size(V_cell,1)
          #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
          find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          ρ_p = Roots.find_zero(find_ρ, 1.0)
          closure = genome_distr[:,i]./sum(genome_distr[:,i])
          ρ_ps[j,i] = ρ_p[1].*closure[2]
          y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          y_DEs[j,i] = y_DE[1]
      end
  elseif df_aromatics.Ontology[j] == "Amino acids"
      for i in 1:size(V_cell,1)
          #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
          find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          ρ_p = Roots.find_zero(find_ρ, 1.0)
          closure = genome_distr[:,i]./sum(genome_distr[:,i])
          ρ_ps[j,i] = ρ_p[1].*closure[3]
          y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          y_DEs[j,i] = y_DE[1]
      end
  elseif df_aromatics.Ontology[j] == "Fatty acids"
      for i in 1:size(V_cell,1)
          #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
          find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          ρ_p = Roots.find_zero(find_ρ, 1.0)
          closure = genome_distr[:,i]./sum(genome_distr[:,i])
          ρ_ps[j,i] = ρ_p[1].*closure[4]
          y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          y_DEs[j,i] = y_DE[1]
      end
  elseif df_aromatics.Ontology[j] == "Nucleotides"
      for i in 1:size(V_cell,1)
          #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
          find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          ρ_p = Roots.find_zero(find_ρ, 1.0)
          closure = genome_distr[:,i]./sum(genome_distr[:,i])
          ρ_ps[j,i] = ρ_p[1].*closure[5]
          y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          y_DEs[j,i] = y_DE[1]
      end
  else
      for i in 1:size(V_cell,1)
          #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]],[rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
          find_ρ(x) = DEBmicroTrait.constrain_transporter_density_cost(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          ρ_p = Roots.find_zero(find_ρ, 1.0)
          closure = genome_distr[:,i]./sum(genome_distr[:,i])
          ρ_ps[j,i] = ρ_p[1].*closure[6]
          y_DE = DEBmicroTrait.yield_transporter_density_cost(ρ_ps[j,i], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], [y_EM[i]], df_aromatics.Formula[j])
          y_DEs[j,i] = y_DE[1]
      end
  end
end

ρ_ps[ρ_ps.==0.0] .= 1e-12

N_C = zeros(size(df_aromatics.Name,1))
for i in 1:size(df_aromatics.Name,1)
  elementstring = df_aromatics.Formula[i]
  N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
end

N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C

D_S               = DEBmicroTrait.aqueous_diffusivity(df_aromatics.Molecular_weight)
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, D_S)

df_aromatics_out = DataFrame()
df_aromatics_out.Vmax = vec(Vmax)
df_aromatics_out.KD = vec(K_D)
df_aromatics_out.affinity = vec(Vmax./K_D)

df_aromatics_out.ontology = repeat(df_aromatics.Ontology, size(V_cell,1))
response = Array{String}(undef, (size(df_aromatics.Name,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     response[:,i] .= df_isolates.Rhizosphere_response[i]
 end

df_aromatics_out.response = vec(response)

df_aromatics_out.yield = vec(abs.(y_DEs))
df_aromatics_out.name = repeat(df_aromatics.Name, size(V_cell,1))

df_aromatics_out.transporter_density = vec(ρ_ps)

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_assimilation_aromatics.csv", df_aromatics_out)
