using DEBmicroTrait
using CSV, DataFrames
using Roots
using Statistics
using HypothesisTests

df     = CSV.read("/Users/glmarschmann/Data/Zhen/IsogenieGenomes.ecosysguilds.csv", DataFrame)

################################################################################
###################### define functional guilds
# aerobic heterotroph
df_aerob_hetero                 = df[(df.heterotroph_withETC.==1) , :]
Genome_size_aerob_hetero        = convert(Array{Float64}, df_aerob_hetero.bp)
V_cell_aerob_hetero             = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_aerob_hetero)
Min_gen_time_aerob_hetero       = df_aerob_hetero.minimum_generation_time_hours
rrn_copies_aerob_hetero         = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_aerob_hetero)
Gram_stain_aerob_hetero         = repeat(["+"], size(V_cell_aerob_hetero,1))

# acetotrophic methanogens
df_am = df[(df.acetoclastic_methanogenesis.==1) , :]
Genome_size_am = convert(Array{Float64}, df_am.bp)
V_cell_am = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_am)
Min_gen_time_am = df_am.minimum_generation_time_hours
rrn_copies_am = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_am)
Gram_stain_am = repeat(["+"], size(V_cell_am,1))
SAV_am = DEBmicroTrait.surface_area_volume_ratio(V_cell_am./1e-18)

# acetotrophic methanogens methanol
df_amm = df_am[(df_am.methanol_methylCoM_or.==1) , :]
Genome_size_amm = convert(Array{Float64}, df_amm.bp)
V_cell_amm = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_amm)
Min_gen_time_amm = df_amm.minimum_generation_time_hours
rrn_copies_amm = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_amm)
Gram_stain_amm = repeat(["+"], size(V_cell_amm,1))
SAV_amm = DEBmicroTrait.surface_area_volume_ratio(V_cell_amm./1e-18)

# hydrogenotrophic methanogens
df_hmo = df[(df.hydrogenotrophic_methanogenesis.==1) , :]
Genome_size_hmo = convert(Array{Float64}, df_hmo.bp)
V_cell_hmo = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_hmo)
Min_gen_time_hmo = df_hmo.minimum_generation_time_hours
rrn_copies_hmo = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_hmo)
Gram_stain_hmo = repeat(["+"], size(V_cell_hmo,1))
SAV_hmo = DEBmicroTrait.surface_area_volume_ratio(V_cell_hmo./1e-18)

# methane oxidizers
df_mo = df[(df.methane_oxidation.==1) , :]
Genome_size_mo = convert(Array{Float64}, df_mo.bp)
V_cell_mo = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_mo)
Min_gen_time_mo = df_mo.minimum_generation_time_hours
rrn_copies_mo = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_mo)
Gram_stain_mo = repeat(["+"], size(V_cell_mo,1))
SAV_mo = DEBmicroTrait.surface_area_volume_ratio(V_cell_mo./1e-18)

# Fermenters
df_ferm = df[(df.heterotroph_withETC.==0) , :]
Genome_size_ferm = convert(Array{Float64}, df_ferm.bp)
V_cell_ferm = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_ferm)
Min_gen_time_ferm = df_ferm.minimum_generation_time_hours
rrn_copies_ferm = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_ferm)
Gram_stain_ferm = repeat(["+"], size(V_cell_ferm,1))
SAV_ferm = DEBmicroTrait.surface_area_volume_ratio(V_cell_ferm./1e-18)

# Nitrate reductase
df_nitrate_red = df[(df.denitrification.==1) , :]
Genome_size_nitrate_red = convert(Array{Float64}, df_nitrate_red.bp)
V_cell_nitrate_red = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitrate_red)
Min_gen_time_nitrate_red = df_nitrate_red.minimum_generation_time_hours
rrn_copies_nitrate_red = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitrate_red)
Gram_stain_nitrate_red= repeat(["+"], size(V_cell_nitrate_red,1))
SAV_nitrate_red = DEBmicroTrait.surface_area_volume_ratio(V_cell_nitrate_red./1e-18)

# Nitrite reductase
df_nitrite_red = df[(df.nitrite_nitricoxide_partial.==1) , :]
Genome_size_nitrite_red = convert(Array{Float64}, df_nitrite_red.bp)
V_cell_nitrite_red = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitrite_red)
Min_gen_time_nitrite_red = df_nitrite_red.minimum_generation_time_hours
rrn_copies_nitrite_red = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitrite_red)
Gram_stain_nitrite_red= repeat(["+"], size(V_cell_nitrite_red,1))
SAV_nitrite_red = DEBmicroTrait.surface_area_volume_ratio(V_cell_nitrite_red./1e-18)

# Nitric oxide reductase
df_nitoxide_red = df[(df.nitricoxide_nitrousoxide_partial.==1) , :]
Genome_size_nitoxide_red  = convert(Array{Float64}, df_nitoxide_red.bp)
V_cell_nitoxide_red  = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitoxide_red )
Min_gen_time_nitoxide_red  = df_nitoxide_red.minimum_generation_time_hours
rrn_copies_nitoxide_red  = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitoxide_red )
Gram_stain_nitoxide_red = repeat(["+"], size(V_cell_nitoxide_red ,1))
SAV_nitoxide_red  = DEBmicroTrait.surface_area_volume_ratio(V_cell_nitoxide_red./1e-18)

# Nitrous oxide reductase
df_nitrousoxide_red = df[(df.nitrousoxide_dinitrogen.==1) , :]
Genome_size_nitrousoxide_red  = convert(Array{Float64}, df_nitrousoxide_red.bp)
V_cell_nitrousoxide_red  = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitrousoxide_red)
Min_gen_time_nitrousoxide_red  = df_nitrousoxide_red.minimum_generation_time_hours
rrn_copies_nitrousoxide_red  = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitrousoxide_red)
Gram_stain_nitrousoxide_red = repeat(["+"], size(V_cell_nitrousoxide_red ,1))
SAV_nitrousoxide_red  = DEBmicroTrait.surface_area_volume_ratio(V_cell_nitrousoxide_red./1e-18)

# Complete denitrification
df_nitrousoxide_red = df[(df.denitrification.==1), :]
Genome_size_nitrousoxide_red  = convert(Array{Float64}, df_nitrousoxide_red.bp)
V_cell_nitrousoxide_red  = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_nitrousoxide_red)
Min_gen_time_nitrousoxide_red  = df_nitrousoxide_red.minimum_generation_time_hours
rrn_copies_nitrousoxide_red  = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_nitrousoxide_red)
Gram_stain_nitrousoxide_red = repeat(["+"], size(V_cell_nitrousoxide_red ,1))
SAV_nitrousoxide_red  = DEBmicroTrait.surface_area_volume_ratio(V_cell_nitrousoxide_red./1e-18)

# aerobic diazotroph
df_aerob_diazo = df[(df.aerobic_diazotroph.==1), :]
Genome_size_aerob_diazo  = convert(Array{Float64}, df_aerob_diazo.bp)
V_cell_aerob_diazo  = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_aerob_diazo)
Min_gen_time_aerob_diazo  = df_aerob_diazo.minimum_generation_time_hours
rrn_copies_aerob_diazo  = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_aerob_diazo)
Gram_stain_aerob_diazo = repeat(["+"], size(V_cell_aerob_diazo ,1))
SAV_aerob_diazo  = DEBmicroTrait.surface_area_volume_ratio(V_cell_aerob_diazo./1e-18)

# anaerobic diazotroph
df_anaerob_diazo = df[(df.anaerobic_diazotroph.==1), :]
Genome_size_anaerob_diazo  = convert(Array{Float64}, df_anaerob_diazo.bp)
V_cell_anaerob_diazo  = DEBmicroTrait.genome_size_to_cell_volume(Genome_size_anaerob_diazo)
Min_gen_time_anaerob_diazo  = df_anaerob_diazo.minimum_generation_time_hours
rrn_copies_anaerob_diazo  = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size_anaerob_diazo)
Gram_stain_anaerob_diazo = repeat(["+"], size(V_cell_anaerob_diazo ,1))
SAV_aerob_andiazo  = DEBmicroTrait.surface_area_volume_ratio(V_cell_anaerob_diazo./1e-18)
################################################################################

################################################################################
                        # Glucose oxidation- aerobic heterotroph
# transporter estimation
N_C  = [6]
y_DE = [2.61]./N_C
y_EM = ones(size(V_cell_aerob_hetero,1))
ρ_ps_aerob_hetero_glucose = zeros(size(V_cell_aerob_hetero,1))
r_ms = zeros(size(V_cell_aerob_hetero,1))
for i in 1:size(V_cell_aerob_hetero,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_aerob_hetero_glucose[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_aerob_hetero)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_aerob_hetero, Min_gen_time_aerob_hetero, Gram_stain_aerob_hetero)
N_SB_aerob_hetero_glucose      = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_aerob_hetero[i]], ρ_ps_aerob_hetero_glucose[i].*ones(1,1), [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]])[1] for i in 1:size(V_cell_aerob_hetero,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_aerob_hetero_glucose
r_resp_aerob_hetero_glucose    = r_growth .+ r_main .+r_assim
median(r_resp_aerob_hetero_glucose)
std(r_resp_aerob_hetero_glucose)
# half-saturation constant
Molecular_weight       = [180.16]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_aerob_hetero       = DEBmicroTrait.specific_reference_affinity(V_cell_aerob_hetero, reshape(ρ_ps_aerob_hetero_glucose*100, 714,1), D_S)
median(K_D_aerob_hetero.*N_C.*12.011)
std(K_D_aerob_hetero.*N_C*12.011)
K_D_aerob_hetero_glucose = K_D_aerob_hetero.*N_C.*12.011
                        # Acetate oxidation- aerobic heterotroph
# transporter estimation
N_C  = [2]
y_DE = [0.66]./N_C
y_EM = ones(size(V_cell_aerob_hetero,1))
ρ_ps_aerob_hetero_acetate = zeros(size(V_cell_aerob_hetero,1))
r_ms = zeros(size(V_cell_aerob_hetero,1))
for i in 1:size(V_cell_aerob_hetero,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_aerob_hetero_acetate[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_aerob_hetero[i]], [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]], [rrn_copies_aerob_hetero[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV                           = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_aerob_hetero)
r_growth                       = @. (1-1/y_EV)*r_ms
r_main                         = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_aerob_hetero, Min_gen_time_aerob_hetero, Gram_stain_aerob_hetero)
N_SB_aerob_hetero_acetate      = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_aerob_hetero[i]], ρ_ps_aerob_hetero_acetate[i].*ones(1,1), [Min_gen_time_aerob_hetero[i]], [Gram_stain_aerob_hetero[i]])[1] for i in 1:size(V_cell_aerob_hetero,1)]
k2p                            = 180.0*60^2
r_assim                        = @. y_DE*N_C*k2p*N_SB_aerob_hetero_acetate
r_resp_aerob_hetero_acetate    = r_growth .+ r_main .+r_assim
median(r_resp_aerob_hetero_acetate)
std(r_resp_aerob_hetero_acetate)
# half-saturation constant
Molecular_weight       = [59.04]
D_S                    = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_aerob_hetero       = DEBmicroTrait.specific_reference_affinity(V_cell_aerob_hetero, reshape(ρ_ps_aerob_hetero_acetate*100, 714,1), D_S)
median(K_D_aerob_hetero.*N_C.*12.011)
std(K_D_aerob_hetero.*N_C*12.011)
K_D_aerob_hetero_acetate = K_D_aerob_hetero.*N_C.*12.011
################################################################################
r_resp_aerob_hetero = vcat(r_resp_aerob_hetero_glucose, r_resp_aerob_hetero_acetate)
K_D_aerob_hetero = vcat(K_D_aerob_hetero_glucose, K_D_aerob_hetero_acetate)
################################################################################
                        # Glucose fermentation
N_C = [6]
y_DE = [1.06]./N_C
y_EM = ones(size(V_cell_ferm,1))
ρ_ps_ferm= zeros(size(V_cell_ferm,1))
r_ms = zeros(size(V_cell_ferm,1))
for i in 1:size(V_cell_ferm,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_ferm[i]], [Min_gen_time_ferm[i]], [Gram_stain_ferm[i]], [rrn_copies_ferm[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_ferm[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_ferm[i]], [Min_gen_time_ferm[i]], [Gram_stain_ferm[i]], [rrn_copies_ferm[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_ferm)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_ferm, Min_gen_time_ferm, Gram_stain_ferm)
N_SB_ferm = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_ferm[i]], ρ_ps_ferm[i]*ones(1,1), [Min_gen_time_ferm[i]], [Gram_stain_ferm[i]])[1] for i in 1:size(V_cell_ferm,1)]
k2p = 1800.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_ferm
r_resp_ferm = r_growth .+ r_main .+r_assim
median(r_resp_ferm)
std(r_resp_ferm)
# half saturation constant
Molecular_weight = [180.16]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_ferm = DEBmicroTrait.specific_reference_affinity(V_cell_ferm, reshape(ρ_ps_ferm*100, 815,1), D_S)
median(K_D_ferm*6*12.011)
std(K_D_ferm*6*12.011)

minimum(K_D_ferm*6*12.011)
maximum(K_D_ferm*6*12.011)
K_D_ferm = K_D_ferm*6*12.011
################################################################################

################################################################################
                        # Acetotrophic methanogens - acetate
N_C = [2]
y_DE = [1.49]./N_C
y_EM = ones(size(V_cell_am,1))
ρ_ps_am= zeros(size(V_cell_am,1))
r_ms = zeros(size(V_cell_am,1))
for i in 1:size(V_cell_am,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_am[i]], [Min_gen_time_am[i]], [Gram_stain_am[i]], [rrn_copies_am[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_am[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_am[i]], [Min_gen_time_am[i]], [Gram_stain_am[i]], [rrn_copies_am[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_am)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_am, Min_gen_time_am, Gram_stain_am)
N_SB_am = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_am[i]], ρ_ps_am[i]*ones(1,1), [Min_gen_time_am[i]], [Gram_stain_am[i]])[1] for i in 1:size(V_cell_am,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_am
r_resp_am = r_growth .+ r_main .+r_assim
median(r_resp_am)
std(r_resp_am)
mean(r_resp_am[r_resp_am .< 0.3])
median(r_growth)
median(r_assim)
# half saturation constants
Molecular_weight = [59.04]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_am = DEBmicroTrait.specific_reference_affinity(V_cell_am, reshape(ρ_ps_am*100, 65,1), D_S)
median(K_D_am.*N_C*12.011)
std(K_D_am.*N_C*12.011)
K_D_am = K_D_am.*N_C*12.011

################################################################################

################################################################################
                        # Acetotrophic methanogens - methanol
N_C = [1]
y_DE = [0.44]./N_C
y_EM = ones(size(V_cell_am,1))
ρ_ps_amm= zeros(size(V_cell_am,1))
r_msm = zeros(size(V_cell_am,1))
for i in 1:size(V_cell_am,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_am[i]], [Min_gen_time_am[i]], [Gram_stain_am[i]], [rrn_copies_am[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_amm[i] = ρ_p[1]
    r_msm[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_am[i]], [Min_gen_time_am[i]], [Gram_stain_am[i]], [rrn_copies_am[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_am)
r_growth = @. (1-1/y_EV)*r_msm
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_am, Min_gen_time_am, Gram_stain_am)
N_SB_amm = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_am[i]], ρ_ps_amm[i]*ones(1,1), [Min_gen_time_am[i]], [Gram_stain_am[i]])[1] for i in 1:size(V_cell_am,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_amm
r_resp_amm = r_growth .+ r_main .+r_assim
median(r_resp_amm)
std(r_resp_amm)
# half saturation constants
Molecular_weight = [32.04]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_amm = DEBmicroTrait.specific_reference_affinity(V_cell_am, reshape(ρ_ps_amm*100, 65,1), D_S)
median(K_D_amm.*N_C*12.011)
std(K_D_amm.*N_C*12.011)
K_D_amm = K_D_amm.*N_C*12.011
################################################################################
r_resp_aceto = vcat(r_resp_am, r_resp_amm)
K_D_aceto = vcat(K_D_am, K_D_amm)


################################################################################
                        # hydrogenotrophic methanogens
N_C = [1]
y_DE = [0.166]./N_C
y_EM = ones(size(V_cell_hmo,1))
ρ_ps_hmo= zeros(size(V_cell_hmo,1))
r_ms = zeros(size(V_cell_hmo,1))
for i in 1:size(V_cell_hmo,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_hmo[i]], [Min_gen_time_hmo[i]], [Gram_stain_hmo[i]], [rrn_copies_hmo[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_hmo[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_hmo[i]], [Min_gen_time_hmo[i]], [Gram_stain_hmo[i]], [rrn_copies_hmo[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_hmo)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_hmo, Min_gen_time_hmo, Gram_stain_hmo)
N_SB_hmo = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_hmo[i]], ρ_ps_hmo[i]*ones(1,1), [Min_gen_time_hmo[i]], [Gram_stain_hmo[i]])[1] for i in 1:size(V_cell_hmo,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_hmo
r_resp_hmo = r_growth .+ r_main .+r_assim
median(r_resp_hmo)
std(r_resp_hmo)
# half saturation constant
Molecular_weight = [2.016]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_hmo = DEBmicroTrait.specific_reference_affinity(V_cell_mo, reshape(ρ_ps_hmo, 69,1), D_S)
median(K_D_hmo.*N_C*12.011)
std(K_D_hmo.*N_C*12.011)
K_D_hmo = K_D_hmo.*N_C*12.011
################################################################################

################################################################################
                        # hydrogenotrophic methanogens methanol
N_C = [1]
y_DE = [0.44]./N_C
y_EM = ones(size(V_cell_hmo,1))
ρ_ps_hmom= zeros(size(V_cell_hmo,1))
r_msm = zeros(size(V_cell_hmo,1))
for i in 1:size(V_cell_hmo,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_hmo[i]], [Min_gen_time_hmo[i]], [Gram_stain_hmo[i]], [rrn_copies_hmo[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_hmom[i] = ρ_p[1]
    r_msm[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_hmo[i]], [Min_gen_time_hmo[i]], [Gram_stain_hmo[i]], [rrn_copies_hmo[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_hmo)
r_growth = @. (1-1/y_EV)*r_msm
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_hmo, Min_gen_time_hmo, Gram_stain_hmo)
N_SB_hmom = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_hmo[i]], ρ_ps_hmom[i]*ones(1,1), [Min_gen_time_hmo[i]], [Gram_stain_hmo[i]])[1] for i in 1:size(V_cell_hmo,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_hmom
r_resp_hmom = r_growth .+ r_main .+r_assim
median(r_resp_hmom)
std(r_resp_hmom)
# half saturation constant
Molecular_weight = [32.04]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_hmom = DEBmicroTrait.specific_reference_affinity(V_cell_mo, reshape(ρ_ps_hmom, 69,1), D_S)
median(K_D_hmom.*N_C*12.011)
std(K_D_hmom.*N_C*12.011)
K_D_hmom = K_D_hmom.*N_C*12.011
################################################################################
r_resp_hydro = vcat(r_resp_hmo, r_resp_hmom)
K_D_hydro = vcat(K_D_hmo, K_D_hmom)

################################################################################
                        # methanotrophs
N_C = [6]
y_DE = [5.83]./N_C
y_EM = ones(size(V_cell_mo,1))

ρ_ps_mo= zeros(size(V_cell_mo,1))
r_ms = zeros(size(V_cell_mo,1))
for i in 1:size(V_cell_mo,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_mo[i]], [Min_gen_time_mo[i]], [Gram_stain_mo[i]], [rrn_copies_mo[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_mo[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_mo[i]], [Min_gen_time_mo[i]], [Gram_stain_mo[i]], [rrn_copies_mo[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_mo)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_mo, Min_gen_time_mo, Gram_stain_mo)
N_SB_mo = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_mo[i]], ρ_ps_mo[i]*ones(1,1), [Min_gen_time_mo[i]], [Gram_stain_mo[i]])[1] for i in 1:size(V_cell_mo,1)]
k2p = 18.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_mo
r_resp_mo = r_growth .+ r_main .+r_assim
median(r_resp_mo)
std(r_resp_mo)
################################################################################################
# affinity
Molecular_weight = [16.043]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_mo = DEBmicroTrait.specific_reference_affinity(V_cell_mo, reshape(ρ_ps_mo./100, 15,1), D_S)
median(K_D_mo.*1*12.011)
std(K_D_mo*1*12.011)
K_D_mo = K_D_mo*1*12.011
################################################################################
                        # aerobic diazotroph
N_C = [6]
y_DE = [2.61]./N_C
y_EM = ones(size(V_cell_aerob_diazo,1))

ρ_ps_aerob_diazo= zeros(size(V_cell_aerob_diazo,1))
r_ms = zeros(size(V_cell_aerob_diazo,1))
for i in 1:size(V_cell_aerob_diazo,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_aerob_diazo[i]], [Min_gen_time_aerob_diazo[i]], [Gram_stain_aerob_diazo[i]], [rrn_copies_aerob_diazo[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_aerob_diazo[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_aerob_diazo[i]], [Min_gen_time_aerob_diazo[i]], [Gram_stain_aerob_diazo[i]], [rrn_copies_aerob_diazo[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_aerob_diazo)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_aerob_diazo, Min_gen_time_aerob_diazo, Gram_stain_aerob_diazo)
N_SB_aerob_diazo = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_aerob_diazo[i]], ρ_ps_aerob_diazo[i]*ones(1,1), [Min_gen_time_aerob_diazo[i]], [Gram_stain_aerob_diazo[i]])[1] for i in 1:size(V_cell_aerob_diazo,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_aerob_diazo
r_resp_aerob_diazo = r_growth .+ r_main .+r_assim
median(r_resp_aerob_diazo)
std(r_resp_aerob_diazo)
################################################################################################
# affinity
Molecular_weight = [180.16]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_aerob_diazo = DEBmicroTrait.specific_reference_affinity(V_cell_aerob_diazo, reshape(ρ_ps_aerob_diazo, 45,1), D_S)
median(K_D_aerob_diazo.*N_C*12.011)
std(K_D_aerob_diazo*N_C*12.011)
K_D_aerob_diazo=K_D_aerob_diazo*N_C*12.011
################################################################################
                        # anaerobic diazotroph
N_C = [6]
y_DE = [4.36]./N_C
y_EM = ones(size(V_cell_anaerob_diazo,1))

ρ_ps_anaerob_diazo= zeros(size(V_cell_anaerob_diazo,1))
r_ms = zeros(size(V_cell_anaerob_diazo,1))
for i in 1:size(V_cell_anaerob_diazo,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_anaerob_diazo[i]], [Min_gen_time_anaerob_diazo[i]], [Gram_stain_anaerob_diazo[i]], [rrn_copies_anaerob_diazo[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_anaerob_diazo[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_anaerob_diazo[i]], [Min_gen_time_anaerob_diazo[i]], [Gram_stain_anaerob_diazo[i]], [rrn_copies_anaerob_diazo[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_anaerob_diazo)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_anaerob_diazo, Min_gen_time_anaerob_diazo, Gram_stain_anaerob_diazo)
N_SB_anaerob_diazo = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_anaerob_diazo[i]], ρ_ps_anaerob_diazo[i]*ones(1,1), [Min_gen_time_anaerob_diazo[i]], [Gram_stain_anaerob_diazo[i]])[1] for i in 1:size(V_cell_anaerob_diazo,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_anaerob_diazo
r_resp_anaerob_diazo = r_growth .+ r_main .+r_assim
median(r_resp_anaerob_diazo)
std(r_resp_anaerob_diazo)
# half saturation constant
Molecular_weight = [180.16]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_anaerob_diazo = DEBmicroTrait.specific_reference_affinity(V_cell_anaerob_diazo, reshape(ρ_ps_anaerob_diazo, 159,1), D_S)
median(K_D_anaerob_diazo.*N_C*12.011)
std(K_D_anaerob_diazo*N_C*12.011)
K_D_anaerob_diazo=K_D_anaerob_diazo*N_C*12.011
################################################################################

################################################################################
                        # Glucose nitrite reductase
N_C = [6]
y_DE = [5.38]./N_C
y_EM = ones(size(V_cell_nitrite_red,1))

ρ_ps_nitrite_red= zeros(size(V_cell_nitrite_red,1))
r_ms = zeros(size(V_cell_nitrite_red,1))
for i in 1:size(V_cell_nitrite_red,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitrite_red[i]], [Min_gen_time_nitrite_red[i]], [Gram_stain_nitrite_red[i]], [rrn_copies_nitrite_red[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitrite_red[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitrite_red[i]], [Min_gen_time_nitrite_red[i]], [Gram_stain_nitrite_red[i]], [rrn_copies_nitrite_red[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitrite_red)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitrite_red, Min_gen_time_nitrite_red, Gram_stain_nitrite_red)
N_SB_nitrite_red = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitrite_red[i]], ρ_ps_nitrite_red[i]*ones(1,1), [Min_gen_time_nitrite_red[i]], [Gram_stain_nitrite_red[i]])[1] for i in 1:size(V_cell_nitrite_red,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_nitrite_red
r_resp_nitrite_red_glucose = r_growth .+ r_main .+r_assim
median(r_resp_nitrite_red)
std(r_resp_nitrite_red)
# half saturation constant
Molecular_weight = [180.16]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitrite_red = DEBmicroTrait.specific_reference_affinity(V_cell_nitrite_red, reshape(ρ_ps_nitrite_red, 63,1), D_S)
median(K_D_nitrite_red.*N_C*12.011)
std(K_D_nitrite_red*N_C*12.011)
K_D_nitrite_red_glucose = K_D_nitrite_red*N_C*12.011
################################################################################

################################################################################
                        # Glucose nitric oxide reductase
N_C = [6]
y_DE = [5.43]./N_C
y_EM = ones(size(V_cell_nitoxide_red,1))

ρ_ps_nitoxide_red= zeros(size(V_cell_nitoxide_red,1))
r_ms = zeros(size(V_cell_nitoxide_red,1))
for i in 1:size(V_cell_nitoxide_red,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitoxide_red[i]], [Min_gen_time_nitoxide_red[i]], [Gram_stain_nitoxide_red[i]], [rrn_copies_nitoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitoxide_red[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitoxide_red[i]], [Min_gen_time_nitoxide_red[i]], [Gram_stain_nitoxide_red[i]], [rrn_copies_nitoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitoxide_red)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitoxide_red, Min_gen_time_nitoxide_red, Gram_stain_nitoxide_red)
N_SB_nitoxide_red = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitoxide_red[i]], ρ_ps_nitoxide_red[i]*ones(1,1), [Min_gen_time_nitoxide_red[i]], [Gram_stain_nitoxide_red[i]])[1] for i in 1:size(V_cell_nitoxide_red,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_nitoxide_red
r_resp_nitoxide_red_glucose = r_growth .+ r_main .+r_assim
median(r_resp_nitoxide_red)
std(r_resp_nitoxide_red)
# half saturation constant
Molecular_weight = [180.16]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitoxide_red = DEBmicroTrait.specific_reference_affinity(V_cell_nitoxide_red, reshape(ρ_ps_nitoxide_red, 200,1), D_S)
median(K_D_nitoxide_red.*N_C*12.011)
std(K_D_nitoxide_red*N_C*12.011)
K_D_nitoxide_red_glucose = K_D_nitoxide_red*N_C*12.011
################################################################################

################################################################################
                        # Glucose nitrous oxide reductase
N_C = [6]
y_DE = [4.19]./N_C
y_EM = ones(size(V_cell_nitrousoxide_red,1))

ρ_ps_nitrousoxide_red= zeros(size(V_cell_nitrousoxide_red,1))
r_ms = zeros(size(V_cell_nitrousoxide_red,1))
for i in 1:size(V_cell_nitrousoxide_red,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitrousoxide_red[i]], [Min_gen_time_nitrousoxide_red[i]], [Gram_stain_nitrousoxide_red[i]], [rrn_copies_nitrousoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitrousoxide_red[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitrousoxide_red[i]], [Min_gen_time_nitrousoxide_red[i]], [Gram_stain_nitrousoxide_red[i]], [rrn_copies_nitrousoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitrousoxide_red)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitrousoxide_red, Min_gen_time_nitrousoxide_red, Gram_stain_nitrousoxide_red)
N_SB_nitrousoxide_red = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitrousoxide_red[i]], ρ_ps_nitrousoxide_red[i]*ones(1,1), [Min_gen_time_nitrousoxide_red[i]], [Gram_stain_nitrousoxide_red[i]])[1] for i in 1:size(V_cell_nitrousoxide_red,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_nitrousoxide_red
r_resp_nitrousoxide_red_glucose = r_growth .+ r_main .+r_assim
median(r_resp_nitrousoxide_red)
std(r_resp_nitrousoxide_red)
# half saturation constant
Molecular_weight = [180.16]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitrousoxide_red = DEBmicroTrait.specific_reference_affinity(V_cell_nitrousoxide_red, reshape(ρ_ps_nitrousoxide_red, 57,1), D_S)
median(K_D_nitrousoxide_red.*N_C*12.011)
std(K_D_nitrousoxide_red*N_C*12.011)
K_D_nitrousoxide_red_glucose = K_D_nitrousoxide_red*N_C*12.011
################################################################################

################################################################################
                        # Acetate nitrite reductase
N_C = [2]
y_DE = [0.38]./N_C
y_EM = ones(size(V_cell_nitrite_red,1))

ρ_ps_nitrite_red= zeros(size(V_cell_nitrite_red,1))
r_ms = zeros(size(V_cell_nitrite_red,1))
for i in 1:size(V_cell_nitrite_red,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitrite_red[i]], [Min_gen_time_nitrite_red[i]], [Gram_stain_nitrite_red[i]], [rrn_copies_nitrite_red[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitrite_red[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitrite_red[i]], [Min_gen_time_nitrite_red[i]], [Gram_stain_nitrite_red[i]], [rrn_copies_nitrite_red[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitrite_red)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitrite_red, Min_gen_time_nitrite_red, Gram_stain_nitrite_red)
N_SB_nitrite_red = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitrite_red[i]], ρ_ps_nitrite_red[i]*ones(1,1), [Min_gen_time_nitrite_red[i]], [Gram_stain_nitrite_red[i]])[1] for i in 1:size(V_cell_nitrite_red,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_nitrite_red
r_resp_nitrite_red_acetate = r_growth .+ r_main .+r_assim
median(r_resp_nitrite_red)
std(r_resp_nitrite_red)
# half saturation constant
Molecular_weight = [59.04]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitrite_red = DEBmicroTrait.specific_reference_affinity(V_cell_nitrite_red, reshape(ρ_ps_nitrite_red, 63,1), D_S)
median(K_D_nitrite_red.*N_C*12.011)
std(K_D_nitrite_red*N_C*12.011)
K_D_nitrite_red_acetate = K_D_nitrite_red*N_C*12.011
################################################################################

################################################################################
                        # Acetate nitric oxide reductase
N_C = [2]
y_DE = [0.88]./N_C
y_EM = ones(size(V_cell_nitoxide_red,1))

ρ_ps_nitoxide_red= zeros(size(V_cell_nitoxide_red,1))
r_ms = zeros(size(V_cell_nitoxide_red,1))
for i in 1:size(V_cell_nitoxide_red,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitoxide_red[i]], [Min_gen_time_nitoxide_red[i]], [Gram_stain_nitoxide_red[i]], [rrn_copies_nitoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitoxide_red[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitoxide_red[i]], [Min_gen_time_nitoxide_red[i]], [Gram_stain_nitoxide_red[i]], [rrn_copies_nitoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitoxide_red)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitoxide_red, Min_gen_time_nitoxide_red, Gram_stain_nitoxide_red)
N_SB_nitoxide_red = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitoxide_red[i]], ρ_ps_nitoxide_red[i]*ones(1,1), [Min_gen_time_nitoxide_red[i]], [Gram_stain_nitoxide_red[i]])[1] for i in 1:size(V_cell_nitoxide_red,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_nitoxide_red
r_resp_nitoxide_red_acetate = r_growth .+ r_main .+r_assim
median(r_resp_nitoxide_red)
std(r_resp_nitoxide_red)
# half saturation constant
Molecular_weight = [59.04]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitoxide_red = DEBmicroTrait.specific_reference_affinity(V_cell_nitoxide_red, reshape(ρ_ps_nitoxide_red, 200,1), D_S)
median(K_D_nitoxide_red.*N_C*12.011)
std(K_D_nitoxide_red*N_C*12.011)
K_D_nitoxide_red_acetate= K_D_nitoxide_red*N_C*12.011
################################################################################

################################################################################
                        # Acetate nitrous oxide reductase
N_C = [2]
y_DE = [1.14]./N_C
y_EM = ones(size(V_cell_nitrousoxide_red,1))

ρ_ps_nitrousoxide_red= zeros(size(V_cell_nitrousoxide_red,1))
r_ms = zeros(size(V_cell_nitrousoxide_red,1))
for i in 1:size(V_cell_nitrousoxide_red,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell_nitrousoxide_red[i]], [Min_gen_time_nitrousoxide_red[i]], [Gram_stain_nitrousoxide_red[i]], [rrn_copies_nitrousoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps_nitrousoxide_red[i] = ρ_p[1]
    r_ms[i] = DEBmicroTrait.check_growth_rate(ρ_p[1], [V_cell_nitrousoxide_red[i]], [Min_gen_time_nitrousoxide_red[i]], [Gram_stain_nitrousoxide_red[i]], [rrn_copies_nitrousoxide_red[i]], y_DE[1], N_C[1], [y_EM[i]])
end
# respiration
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies_nitrousoxide_red)
r_growth = @. (1-1/y_EV)*r_ms
r_main   = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell_nitrousoxide_red, Min_gen_time_nitrousoxide_red, Gram_stain_nitrousoxide_red)
N_SB_nitrousoxide_red = [DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell_nitrousoxide_red[i]], ρ_ps_nitrousoxide_red[i]*ones(1,1), [Min_gen_time_nitrousoxide_red[i]], [Gram_stain_nitrousoxide_red[i]])[1] for i in 1:size(V_cell_nitrousoxide_red,1)]
k2p = 180.0*60^2
r_assim = @. y_DE*N_C*k2p*N_SB_nitrousoxide_red
r_resp_nitrousoxide_red_acetate = r_growth .+ r_main .+r_assim
median(r_resp_nitrousoxide_red)
std(r_resp_nitrousoxide_red)
# half saturation constant
Molecular_weight = [59.04]
D_S = DEBmicroTrait.aqueous_diffusivity(Molecular_weight)
K_D_nitrousoxide_red = DEBmicroTrait.specific_reference_affinity(V_cell_nitrousoxide_red, reshape(ρ_ps_nitrousoxide_red, 57,1), D_S)
median(K_D_nitrousoxide_red.*N_C*12.011)
std(K_D_nitrousoxide_red*N_C*12.011)
K_D_nitrousoxide_red_acetate = K_D_nitrousoxide_red*N_C*12.011
################################################################################

df_out = DataFrame()

df_out.functionalgroup = repeat(["Acetotrophic methanogen"],size(df_am.genome,1))
df_out.rrespg =

df_out.genome = vcat(df_aerob_hetero.genome, df_aerob_hetero.genome, df_ferm.genome, df_am.genome, df_am.genome, df_hmo.genome, df_hmo.genome, df_mo.genome)
df_out.habitat = vcat(df_aerob_hetero.habitat_binnedfrom, df_aerob_hetero.habitat_binnedfrom, df_ferm.habitat_binnedfrom, df_am.habitat_binnedfrom,  df_am.habitat_binnedfrom, df_hmo.habitat_binnedfrom, df_hmo.habitat_binnedfrom, df_mo.habitat_binnedfrom)
df_out.rresp = vcat(r_resp_aerob_hetero, r_resp_ferm, r_resp_aceto, r_resp_hydro, r_resp_mo)
df_out.K = vcat(K_D_aerob_hetero[:,1], K_D_ferm[:,1], K_D_aceto[:,1], K_D_hydro[:,1], K_D_mo[:,1])


functional_groups = vcat(repeat(["Aerobic heterotroph"],size(df_aerob_hetero.genome,1)), repeat(["Aerobic heterotroph"],size(df_aerob_hetero.genome,1)), repeat(["Fermenter"],size(df_ferm.genome,1)), repeat(["Acetotrophic methanogen"],size(df_am.genome,1)), repeat(["Acetotrophic methanogen"],size(df_am.genome,1)),
                        repeat(["Hydrogenotrophic methanogen"],size(df_hmo.genome,1)), repeat(["Hydrogenotrophic methanogen"],size(df_hmo.genome,1)),repeat(["Methanotroph"],size(df_mo.genome,1)))

df_out.functionalgroup = functional_groups

CSV.write("/Users/glmarschmann/Data/Zhen/TraitPredictionsIsogenieGenomes.csv", df_out)
################################################################################

################################################################################
                              # Literature comparison
df_lit     = dropmissing(CSV.read("/Users/glmarschmann/Data/Zhen/LiteratureTraits.csv", DataFrame, missingstring="NaN"))

df_lit_aerob_hetero = filter(x->(x.functionalgroup.=="Aerobic heterotroph"), df_lit)
df_lit_ferm = filter(x->(x.functionalgroup.=="Fermenter"), df_lit)
df_lit_aceto = filter(x->(x.functionalgroup.=="Acetotrophic methanogen"), df_lit)
df_lit_hydro = filter(x->(x.functionalgroup.=="Hydrogenotrophic methanogen"), df_lit)
df_lit_meth = filter(x->(x.functionalgroup.=="Methanotroph"), df_lit)

df_pred_aerob_hetero = filter(x->(x.functionalgroup.=="Aerobic heterotroph"), df_out)
df_pred_ferm = filter(x->(x.functionalgroup.=="Fermenter"), df_out)
df_pred_aceto = filter(x->(x.functionalgroup.=="Acetotrophic methanogen"), df_out)
df_pred_hydro = filter(x->(x.functionalgroup.=="Hydrogenotrophic methanogen"), df_out)
df_pred_meth = filter(x->(x.functionalgroup.=="Methanotroph"), df_out)

HypothesisTests.KruskalWallisTest(df_pred_aerob_hetero.rresp, df_lit_aerob_hetero.rresp)
HypothesisTests.KruskalWallisTest(df_pred_aerob_hetero.K, df_lit_aerob_hetero.K)
HypothesisTests.KruskalWallisTest(df_pred_ferm.rresp, df_lit_ferm.rresp)
HypothesisTests.KruskalWallisTest(df_pred_ferm.K, df_lit_ferm.K)
HypothesisTests.KruskalWallisTest(df_pred_aceto.rresp, df_lit_aceto.rresp)
HypothesisTests.KruskalWallisTest(df_pred_aceto.K, df_lit_aceto.K)
HypothesisTests.KruskalWallisTest(df_pred_hydro.rresp, df_lit_hydro.rresp)
HypothesisTests.KruskalWallisTest(df_pred_hydro.K, df_lit_hydro.K)
HypothesisTests.KruskalWallisTest(df_pred_meth.rresp, df_lit_meth.rresp)
HypothesisTests.KruskalWallisTest(df_pred_meth.K, df_lit_meth.K)

(median(df_pred_aerob_hetero.rresp) - 0.125)/0.125
(median(df_pred_ferm.rresp) - 0.125)/0.125
(median(df_pred_aceto.rresp) - 0.125)/0.125
(median(df_pred_hydro.rresp) - 0.125)/0.125
(median(df_pred_meth.rresp) - 0.375)/0.375

(median(df_pred_aerob_hetero.K) - 12)/12
(median(df_pred_ferm.K) - 12)/12
(median(df_pred_aceto.K) - 12)/12
(median(df_pred_hydro.K) - 0.01)/0.01
(median(df_pred_meth.K) - 0.0012)/0.0012



df_fen_aerob_hetero = filter(x->(x.functionalgroup.=="Aerobic heterotroph" && x.habitat.=="fen"), df_out)
df_fen_ferm = filter(x->(x.functionalgroup.=="Fermenter"  && x.habitat.=="fen"), df_out)
df_fen_aceto = filter(x->(x.functionalgroup.=="Acetotrophic methanogen"  && x.habitat.=="fen"), df_out)
df_fen_hydro = filter(x->(x.functionalgroup.=="Hydrogenotrophic methanogen"  && x.habitat.=="fen"), df_out)
df_fen_meth = filter(x->(x.functionalgroup.=="Methanotroph"  && x.habitat.=="fen"), df_out)

df_bog_aerob_hetero = filter(x->(x.functionalgroup.=="Aerobic heterotroph" && x.habitat.=="bog"), df_out)
df_bog_ferm = filter(x->(x.functionalgroup.=="Fermenter"  && x.habitat.=="bog"), df_out)
df_bog_aceto = filter(x->(x.functionalgroup.=="Acetotrophic methanogen"  && x.habitat.=="bog"), df_out)
df_bog_hydro = filter(x->(x.functionalgroup.=="Hydrogenotrophic methanogen"  && x.habitat.=="bog"), df_out)
df_bog_meth = filter(x->(x.functionalgroup.=="Methanotroph"  && x.habitat.=="bog"), df_out)

df_palsa_aerob_hetero = filter(x->(x.functionalgroup.=="Aerobic heterotroph" && x.habitat.=="palsa"), df_out)
df_palsa_ferm = filter(x->(x.functionalgroup.=="Fermenter"  && x.habitat.=="palsa"), df_out)
df_palsa_aceto = filter(x->(x.functionalgroup.=="Acetotrophic methanogen"  && x.habitat.=="palsa"), df_out)
df_palsa_hydro = filter(x->(x.functionalgroup.=="Hydrogenotrophic methanogen"  && x.habitat.=="palsa"), df_out)
df_palsa_meth = filter(x->(x.functionalgroup.=="Methanotroph"  && x.habitat.=="palsa"), df_out)

HypothesisTests.KruskalWallisTest(df_fen_aerob_hetero.rresp, df_bog_aerob_hetero.rresp)
HypothesisTests.KruskalWallisTest(df_fen_aerob_hetero.rresp, df_palsa_aerob_hetero.rresp)
HypothesisTests.KruskalWallisTest(df_palsa_aerob_hetero.rresp, df_bog_aerob_hetero.rresp)

HypothesisTests.KruskalWallisTest(df_fen_ferm.rresp, df_bog_ferm.rresp)
HypothesisTests.KruskalWallisTest(df_fen_ferm.rresp, df_palsa_ferm.rresp)
HypothesisTests.KruskalWallisTest(df_palsa_ferm.rresp, df_bog_ferm.rresp)

HypothesisTests.KruskalWallisTest(df_fen_aceto.rresp, df_bog_aceto.rresp)
HypothesisTests.KruskalWallisTest(df_fen_aceto.rresp, df_palsa_aceto.rresp)
HypothesisTests.KruskalWallisTest(df_palsa_aceto.rresp, df_bog_aceto.rresp)

HypothesisTests.KruskalWallisTest(df_fen_hydro.rresp, df_bog_hydro.rresp)
HypothesisTests.KruskalWallisTest(df_fen_hydro.rresp, df_palsa_hydro.rresp)
HypothesisTests.KruskalWallisTest(df_palsa_hydro.rresp, df_bog_hydro.rresp)

HypothesisTests.KruskalWallisTest(df_fen_meth.rresp, df_bog_meth.rresp)
HypothesisTests.KruskalWallisTest(df_fen_meth.rresp, df_palsa_meth.rresp)
HypothesisTests.KruskalWallisTest(df_palsa_meth.rresp, df_bog_meth.rresp)


HypothesisTests.KruskalWallisTest(df_fen_aerob_hetero.K, df_bog_aerob_hetero.K)
HypothesisTests.KruskalWallisTest(df_fen_aerob_hetero.K, df_palsa_aerob_hetero.K)
HypothesisTests.KruskalWallisTest(df_palsa_aerob_hetero.K, df_bog_aerob_hetero.K)

HypothesisTests.KruskalWallisTest(df_fen_ferm.K, df_bog_ferm.K)
HypothesisTests.KruskalWallisTest(df_fen_ferm.K, df_palsa_ferm.K)
HypothesisTests.KruskalWallisTest(df_palsa_ferm.K, df_bog_ferm.K)

HypothesisTests.KruskalWallisTest(df_fen_aceto.K, df_bog_aceto.K)
HypothesisTests.KruskalWallisTest(df_fen_aceto.K, df_palsa_aceto.K)
HypothesisTests.KruskalWallisTest(df_palsa_aceto.K, df_bog_aceto.K)

HypothesisTests.KruskalWallisTest(df_fen_hydro.K, df_bog_hydro.K)
HypothesisTests.KruskalWallisTest(df_fen_hydro.K, df_palsa_hydro.K)
HypothesisTests.KruskalWallisTest(df_palsa_hydro.K, df_bog_hydro.K)

HypothesisTests.KruskalWallisTest(df_fen_meth.K, df_bog_meth.K)
HypothesisTests.KruskalWallisTest(df_fen_meth.K, df_palsa_meth.K)
HypothesisTests.KruskalWallisTest(df_palsa_meth.K, df_bog_meth.K)
