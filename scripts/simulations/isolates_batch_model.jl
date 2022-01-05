using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations
using HypothesisTests

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudation_properties.csv", DataFrame, missingstring="N/A")
id_isolates     = 1:39
id_isolates_p   = id_isolates[df_isolates.Rhizosphere_response .== "positive"]
id_isolates_n   = id_isolates[df_isolates.Rhizosphere_response .== "negative"]
id_monomers     = 1:83
id_monomers_sugars = id_monomers[df_metabolites.Ontology .== "Sugars"]
id_monomers_auxins = id_monomers[df_metabolites.Ontology .== "Auxins"]
id_monomers_fattys = id_monomers[df_metabolites.Ontology .== "Fatty acids"]
id_monomers_nucleos = id_monomers[df_metabolites.Ontology .== "Nucleotides"]
id_monomers_aminos = id_monomers[df_metabolites.Ontology .== "Amino acids"]
id_monomers_organics = id_monomers[df_metabolites.Ontology .== "Organic acids"]


assimilation      = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_assimilation.jld")
enzymes           = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_enzymes.jld")
maintenance       = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_maintenance.jld")
protein_synthesis = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_protein_synthesis.jld")
turnover          = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_turnover.jld")
initb             = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_init.jld")

n_transporters    = @. (df_isolates.z_sugars + df_isolates.z_auxins + df_isolates.z_fatty_acids + df_isolates.z_nucleotides + df_isolates.z_amino_acids + df_isolates.z_organic_acids)/df_isolates.Genome_size*1e6

df_pca            = DataFrame()
df_pca.genomesize = df_isolates.Genome_size
df_pca.mingentime = df_isolates.Min_gen_time
df_pca.rrncopies  = df_isolates.rRNA_genes
df_pca.zporters   = round.(n_transporters)
df_pca.zsugars    = round.(df_isolates.z_sugars./(df_isolates.Genome_size/1e6))
df_pca.zauxins    = round.(df_isolates.z_auxins./(df_isolates.Genome_size/1e6))
df_pca.zfattys    = round.(df_isolates.z_fatty_acids./(df_isolates.Genome_size/1e6))
df_pca.znucleos   = round.(df_isolates.z_nucleotides./(df_isolates.Genome_size/1e6))
df_pca.zaminos    = round.(df_isolates.z_amino_acids./(df_isolates.Genome_size/1e6))
df_pca.zorganics  = round.(df_isolates.z_organic_acids./(df_isolates.Genome_size/1e6))
df_pca.zhydrolase = round.(df_isolates.z_hydrolases./(df_isolates.Genome_size/1e6))
V_cell            = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
df_pca.sav        = DEBmicroTrait.surface_area_volume_ratio(V_cell)./1e6

df_pca.response   = df_isolates.Rhizosphere_response
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_pca.csv", df_pca)

df_names = DataFrame()
df_names.Name = df_isolates.Abbreviation
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_name.csv", df_names)


condition(u,t,integrator) = u[1] - 1e-5
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)

BGE_tseries       = zeros(39, 83, 500)
r_tseries         = zeros(39, 83, 500)
BR_tseries        = zeros(39, 83, 500)
BP_tseries        = zeros(39, 83, 500)
x_tseries         = zeros(39, 83, 500)
t_tseries         = zeros(39, 83, 500)
D_tseries         = zeros(39, 83, 500)
E_tseries         = zeros(39, 83, 500)
V_tseries         = zeros(39, 83, 500)
X_tseries         = zeros(39, 83, 500)
CO2_tseries       = zeros(39, 83, 500)
N_cells_tseries   = zeros(39, 83, 500)
assimilation_tseries    = zeros(39, 83, 500)
uptake_tseries    = zeros(39, 83, 500)
maintenance_tseries    = zeros(39, 83, 500)

for i in 1:39
    for j in 1:83
        id_isolate = i
        id_monomer = j

        p                 = DEBmicroTrait.init_batch_model(id_isolate, id_monomer, assimilation, enzymes, maintenance, protein_synthesis, turnover)
        n_polymers        = p.setup_pars.n_polymers
        n_monomers        = p.setup_pars.n_monomers
        n_microbes        = p.setup_pars.n_microbes


        u0                                                                         = zeros(p.setup_pars.dim)
        u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*initb["Bio0"][id_isolate]
        u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*initb["Bio0"][id_isolate]
        u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.25

        tspan             = (0.0,1000.0)
        prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
        sol               = solve(prob, alg_hints=[:stiff], callback=cb)

        for k in 1:length(sol.t)
            t_tseries[i,j,k] = sol.t[k]
            D_tseries[i,j,k] =  sol[k][1]
            E_tseries[i,j,k] =  sol[k][2]
            V_tseries[i,j,k] =  sol[k][3]
            X_tseries[i,j,k] =  sol[k][4]
            CO2_tseries[i,j,k] =  sol[k][5]
            Bio = sol[k][2].+sol[k][3]
            N_cells_tseries[i,j,k]  = @. Bio[1]*1e-6*12.011/(initb["rhoB"]*initb["Md"])[1]
        end
        du   = zeros(p.setup_pars.dim)
        BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
        BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
        BGE  = @. BP/(BP + BR)
        for k in 1:length(sol.t)
            BGE_tseries[i,j,k] = BGE[k]
        end
        for k in 1:length(sol.t)
            BP_tseries[i,j,k] = BP[k]
        end
        for k in 1:length(sol.t)
            BR_tseries[i,j,k] = BR[k]
        end
        r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [sol[i][2]], [sol[i][3]])[1] for i in 1:size(sol.t,1)]
        for k in 1:length(sol.t)
            r_tseries[i,j,k] = r[k]
        end

        for k in 1:length(sol.t)
            x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r[k], p.metabolism_pars, [sol[k][2]], [sol[k][3]])
            x_tseries[i,j,k] = x[1]
        end

        for k in 1:length(sol.t)
            x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r[k], p.metabolism_pars, [sol[k][2]], [sol[k][3]])
            maintenance_tseries[i,j,k] = rM_CO2[1]./(rG_CO2[1]+rM_CO2[1]+rX_CO2[1])
        end

        for k in 1:length(sol.t)
            J_DE  = DEBmicroTrait.assimilation!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
            assimilation_tseries[i,j,k] = J_DE[1]
        end

        for k in 1:length(sol.t)
            J_D  = DEBmicroTrait.uptake!(zeros(1), p.assimilation_pars, [sol[k][1]], [sol[k][3]])
            uptake_tseries[i,j,k] = J_D[1]
        end
    end
end

save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_batch_model.jld", "tseries", t_tseries, "Dseries", D_tseries,
        "Eseries", E_tseries, "Vseries", V_tseries, "Xseries", X_tseries, "CO2series", CO2_tseries, "rseries", r_tseries, "BGEseries", BGE_tseries,
        "xseries", x_tseries, "Nseries", N_cells_tseries)


BGE_tseries[BGE_tseries.>1.0].=0
BGE_tseries[BGE_tseries.<=0.0].=0
BGE_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            BGE_median[i,j]  = median(filter(!iszero, BGE_tseries[i,j,:]))
        catch
            BGE_median[i,j]  = NaN
        end
    end
end

BP_norm_tseries = BP_tseries./N_cells_tseries
BP_norm_tseries[BP_norm_tseries.<0.0] .=NaN

BP_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            BP_median[i,j]  = median(filter(!isnan, BP_norm_tseries[i,j,:]))
        catch
            BP_median[i,j]  = NaN
        end
    end
end



BR_norm_tseries = BR_tseries./N_cells_tseries
BR_norm_tseries[BR_norm_tseries.<0.0] .=NaN
BR_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            BR_median[i,j]  = median(filter(!isnan, BR_norm_tseries[i,j,:]))
        catch
            BR_median[i,j]  = NaN
        end
    end
end



r_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            r_median[i,j]  = median(filter(!iszero, r_tseries[i,j,:]))
        catch
            r_median[i,j]  = NaN
        end
    end
end

x_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            x_median[i,j]  = median(filter(!iszero, x_tseries[i,j,:]))
        catch
            x_median[i,j]  = NaN
        end
    end
end

jM_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            jM_median[i,j]  = median(filter(!iszero, maintenance_tseries[i,j,:]))
        catch
            jM_median[i,j]  = NaN
        end
    end
end



assimilation_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            assimilation_median[i,j]  = median(filter(!iszero, assimilation_tseries[i,j,:]))
        catch
            assimilation_median[i,j]  = NaN
        end
    end
end

uptake_median = zeros(39,83)
for i in 1:39
    for j in 1:83
        try
            uptake_median[i,j]  = median(filter(!iszero, uptake_tseries[i,j,:]))
        catch
            uptake_median[i,j]  = NaN
        end
    end
end

dreserve_median = zeros(39,83)
mE_tseries = E_tseries./V_tseries

for i in 1:39
    for j in 1:83
        try
            dreserve_median[i,j]  = median(filter(!isnan, mE_tseries[i,j,:]))
        catch
            dreserve_median[i,j]  = NaN
        end
    end
end



df_out_bge = DataFrame()
df_out_bge.BGE = vec(BGE_median)    # BGE
df_out_bge.BP = vec(BP_median)    # BP
df_out_bge.BR = vec(BR_median)
df_out_bge.rgrowth = vec(r_median)    # growth rate
df_out_bge.renzyme = vec(x_median)    # enzyme rate
df_out_bge.rmaint = vec(jM_median)    # enzyme rate
df_out_bge.rassim = vec(assimilation_median)
df_out_bge.ruptake = vec(uptake_median)
df_out_bge.rdreserve = vec(dreserve_median)    # enzyme rate
df_out_bge.isolate = repeat(df_isolates.Isolate, 83)  # species
df_out_bge.class = repeat(df_isolates.Class, 83)
df_out_bge.phylum = repeat(df_isolates.Phylum, 83)
df_out_bge.response = repeat(df_isolates.Rhizosphere_response, 83)
df_out_bge.genomesize = repeat(df_isolates.Genome_size, 83)
monomer = Array{String}(undef,39,83)
for i in 1:83
     monomer[:,i] .= df_metabolites.Name[i]
 end
df_out_bge.monomer = vec(monomer)
ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
df_out_bge.ontology = vec(ontology)

rhos = Array{Float64}(undef,39,83)
for i in 1:83
    rhos[:,i] .= assimilation["rho"][i,:]
end
df_out_bge.rho = vec(rhos)

Vmaxs = Array{Float64}(undef,39,83)
for i in 1:83
 Vmaxs[:,i] .= 180.0*60^2*assimilation["NSB"][i,:]
end

df_out_bge.Vmax = vec(Vmaxs)

Ks = Array{Float64}(undef,39,83)
for i in 1:83
    Ks[:,i] .= assimilation["KD"][i,:]
end

df_out_bge.KD = vec(Ks)

yields = Array{Float64}(undef,39,83)
for i in 1:83
    yields[:,i] .= assimilation["yDE"][i,:]
end
df_out_bge.yield = vec(yields)

df_out_bge.kE = repeat(protein_synthesis["kE"], 83)
df_out_bge.yEV = repeat(protein_synthesis["yEV"], 83)
df_out_bge.aX = repeat(enzymes["alpha"], 83)


CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_BGE_all.csv", df_out_bge)

df_out_levin = DataFrame(assimilation_median, :auto)
df_out_levin.response = df_isolates.Rhizosphere_response
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/first_manuscript/files_pub/isolates_batch_model_levin.csv", df_out_levin)


df_out_met = DataFrame()
df_out_met.BGE = vec(BGE_median)
name = Array{String}(undef,39,83)
for i in 1:83
     name[:,i] .= df_metabolites.Name[i]
 end
df_out_met.name = vec(name)

yieldass = Array{String}(undef,39,83)
df_out_met.yieldass = vec(assimilation["yDE"])

df_low = filter(x->(x.BGE.<0.25) , df_out_met)

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_batch_model_BGE_met.csv", df_low)


df_out = DataFrame()

df_out.response = repeat(df_isolates.Rhizosphere_response, 83)
df_out.genomesize = repeat(df_isolates.Genome_size./1e6, 83)
df_out.mingt = repeat(df_isolates.Min_gen_time, 83)
df_out.rrn = repeat(df_isolates.rRNA_genes, 83)
df_out.gramstain = repeat(df_isolates.gram_stain, 83)
df_out.gramstain[df_out.gramstain .== "(+)"] .= "1"
df_out.gramstain[df_out.gramstain .== "(-)"] .= "2"
df_out.kE = repeat(protein_synthesis["kE"], 83)
df_out.yEV = repeat(protein_synthesis["yEV"], 83)
df_out.aX = repeat(enzymes["alpha"], 83)
df_out.gV = repeat(turnover["gV0"], 83)

V_cell = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
SAV = DEBmicroTrait.surface_area_volume_ratio(V_cell)
df_out.sav = repeat(SAV./1e6, 83)

yields = Array{Float64}(undef,39,83)
for i in 1:83
    yields[:,i] .= assimilation["yDE"][i,:]
end
df_out.yield = vec(yields)

df_lambda = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudates_lambda.csv", DataFrame, missingstring="N/A")
lambdas = Array{Float64}(undef,39,83)
for i in 1:83
    lambdas[:,i] .= df_lambda.lambda[i]
end
lambdas
df_out.lambda = vec(lambdas)

rhos = Array{Float64}(undef,39,83)
for i in 1:83
    rhos[:,i] .= assimilation["rho"][i,:]
end
df_out.rho = vec(rhos)

Vmaxs = Array{Float64}(undef,39,83)
for i in 1:83
 Vmaxs[:,i] .= 180.0*60^2*assimilation["NSB"][i,:]
end

df_out.Vmax = vec(Vmaxs)

Ks = Array{Float64}(undef,39,83)
for i in 1:83
    Ks[:,i] .= assimilation["KD"][i,:]
end

df_out.KD = vec(Ks)
df_out.rgrowth = vec(r_median)
df_out.BGE = vec(BGE_median)

df_out_pn = filter(x->(x.response.=="negative" || x.response.=="positive") , df_out)
df_out_u = filter(x->(x.response.=="undefined") , df_out)
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_train_bge_growth.csv", df_out_pn[:,2:17])
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/isolates_batch_model_test_bge_growth.csv", df_out_u[:,2:17])


ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
df_out.ontology = vec(ontology)

df_metabolites.Ontology

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_batch_model_bge.csv", df_out)




metabolite = Array{String}(undef,39,83)
for i in 1:39
     metabolite[i,:] .= df_metabolites.Name[i]
 end
df_out.metabolite = vec(metabolite)

df_out = DataFrame()
df_out.BGE = vec(BGE_median)
df_out.rgrowth = vec(r_median)
df_out.response = repeat(df_isolates.Rhizosphere_response, 83)
ontology = Array{String}(undef,39,83)
for i in 1:83
     ontology[:,i] .= df_metabolites.Ontology[i]
 end
df_out.ontology = vec(ontology)
########################################
# statistics
using HypothesisTests
#sugars

df_p_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="positive") , df_out)
df_n_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="negative") , df_out)
df_u_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="undefined") , df_out)
BGE_sugars_pn  = KruskalWallisTest(df_p_sugars.BGE, df_n_sugars.BGE)
BGE_sugars_p = median(filter(!isnan, df_p_sugars.BGE))
BGE_sugars_n = median(filter(!isnan, df_n_sugars.BGE))
#organics
df_p_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="positive") , df_out)
df_n_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="negative") , df_out)
df_u_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="undefined") , df_out)
BGE_organics_pn  = KruskalWallisTest(df_p_organics.BGE, df_n_organics.BGE)
BGE_organics_p = median(filter(!isnan, df_p_organics.BGE))
BGE_organics_n = median(filter(!isnan, df_n_organics.BGE))
#aminos
df_p_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="positive") , df_out)
df_n_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="negative") , df_out)
df_u_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="undefined") , df_out)
BGE_aminos_pn  = KruskalWallisTest(df_p_aminos.BGE, df_n_aminos.BGE)
BGE_aminos_p = median(filter(!isnan, df_p_aminos.BGE))
BGE_aminos_n = median(filter(!isnan, df_n_aminos.BGE))
# fattys
df_p_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="positive") , df_out)
df_n_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="negative") , df_out)
df_u_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="undefined") , df_out)
BGE_fattys_pn  = KruskalWallisTest(df_p_fattys.BGE, df_n_fattys.BGE)
BGE_fattys_p = median(filter(!isnan, df_p_fattys.BGE))
BGE_fattys_n = median(filter(!isnan, df_n_fattys.BGE))
# nucleos
df_p_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="positive") , df_out)
df_n_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="negative") , df_out)
df_u_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="undefined") , df_out)
BGE_nucleos_pn  = KruskalWallisTest(df_p_nucleos.BGE, df_n_nucleos.BGE)
BGE_nucleos_p = median(filter(!isnan, df_p_nucleos.BGE))
BGE_nucleos_n = median(filter(!isnan, df_n_nucleos.BGE))
# auxins
df_p_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="positive") , df_out)
df_n_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="negative") , df_out)
df_u_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="undefined") , df_out)
BGE_auxins_pn  = KruskalWallisTest(df_p_auxins.BGE, df_n_auxins.BGE)
BGE_auxins_p = median(filter(!isnan, df_p_auxins.BGE))
BGE_auxins_n = median(filter(!isnan, df_n_auxins.BGE))

# growth rates

df_p_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="positive") , df_out)
df_n_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="negative") , df_out)
df_u_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="undefined") , df_out)
BGE_sugars_pn  = KruskalWallisTest(df_p_sugars.rgrowth, df_n_sugars.rgrowth)
BGE_sugars_p = median(filter(!isnan, df_p_sugars.rgrowth))
BGE_sugars_n = median(filter(!isnan, df_n_sugars.rgrowth))
#organics
df_p_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="positive") , df_out)
df_n_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="negative") , df_out)
df_u_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="undefined") , df_out)
BGE_organics_pn  = KruskalWallisTest(df_p_organics.rgrowth, df_n_organics.rgrowth)
BGE_organics_p = median(filter(!isnan, df_p_organics.rgrowth))
BGE_organics_n = median(filter(!isnan, df_n_organics.rgrowth))
#aminos
df_p_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="positive") , df_out)
df_n_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="negative") , df_out)
df_u_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="undefined") , df_out)
BGE_aminos_pn  = KruskalWallisTest(df_p_aminos.rgrowth, df_n_aminos.rgrowth)
BGE_aminos_p = median(filter(!isnan, df_p_aminos.rgrowth))
BGE_aminos_n = median(filter(!isnan, df_n_aminos.rgrowth))
# fattys
df_p_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="positive") , df_out)
df_n_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="negative") , df_out)
df_u_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="undefined") , df_out)
BGE_fattys_pn  = KruskalWallisTest(df_p_fattys.rgrowth, df_n_fattys.rgrowth)
BGE_fattys_p = median(filter(!isnan, df_p_fattys.BGE))
BGE_fattys_n = median(filter(!isnan, df_n_fattys.BGE))
# nucleos
df_p_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="positive") , df_out)
df_n_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="negative") , df_out)
df_u_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="undefined") , df_out)
BGE_nucleos_pn  = KruskalWallisTest(df_p_nucleos.rgrowth, df_n_nucleos.rgrowth)
BGE_nucleos_p = median(filter(!isnan, df_p_nucleos.BGE))
BGE_nucleos_n = median(filter(!isnan, df_n_nucleos.BGE))
# auxins
df_p_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="positive") , df_out)
df_n_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="negative") , df_out)
df_u_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="undefined") , df_out)
BGE_auxins_pn  = KruskalWallisTest(df_p_auxins.rgrowth, df_n_auxins.rgrowth)
BGE_auxins_p = median(filter(!isnan, df_p_auxins.BGE))
BGE_auxins_n = median(filter(!isnan, df_n_auxins.BGE))




df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/final_manuscript/files/isolates_batch_model_bge.csv", DataFrame)
df_low = filter(x->(x.BGE.<0.25) , df)
df_p_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="positive") , df_low)
df_n_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="negative") , df_low)
df_p = filter(x->(x.response.=="positive") , df_low)
df_n = filter(x->(x.response.=="negative") , df_low)
BGE_so  = KruskalWallisTest(df_p.BGE, df_n.BGE)
median(df_p_sugars.BGE)
median(df_p_organics.BGE)






# du = zeros(p.setup_pars.dim)
# mass_balance = zeros(size(sol.t,1))
#
# for i in 1:size(sol.t,1)
#     out = DEBmicroTrait.batch_model!(du, sol.u[i], p, sol.t[i])
#     mass_balance[i] = sum(out)
# end
