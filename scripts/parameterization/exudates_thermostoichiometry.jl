using DEBmicroTrait
using CSV, DataFrames


df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudation_properties.csv", DataFrame, missingstring="N/A")

N_metabolites      = 83
dGcox              = zeros(N_metabolites)
dGcat              = zeros(N_metabolites)
dGAn               = zeros(N_metabolites)
λ_base             = zeros(N_metabolites)
N_C                = zeros(N_metabolites)
eta                = zeros(N_metabolites)
chemFormBiom       = [1, 1.64, 0.2, 0.38, 0, 0, 0]

for i in 1:N_metabolites
    elementstring = convert(String,df_metabolites.Formula[i])
    N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
    out = DEBmicroTrait.get_lambda(elementstring, chemFormBiom)
    dGcox[i] = out[2][3]
    dGcat[i] = out[2][5]
    dGAn[i]  = out[2][8]
    λ_base[i]     = out[1][1]
    eta      = @. dGAn/(λ_base*dGcat)
end


df_out = DataFrame()
df_out.monomer = df_metabolites.Name
df_out.lambda = λ_base
df_out.eta = eta
df_out.ontology = df_metabolites.Ontology

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/exudates_eta.csv", df_out)
