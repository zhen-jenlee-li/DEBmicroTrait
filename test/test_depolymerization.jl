using DEBmicroTrait, Test


n_polymers = 5
n_enzymes  = 6
P          = rand(n_polymers)
E          = rand(n_enzymes)
#
Χ_0        = [10,10,10,10,10]
ρ_P        = [1.35, 1.5, 1.4, 0.4, 1.0]
R_E        = [3.29e-9, 4e-9, 4e-9, 4.4e-9, 5e-9,3.9e-9]
V_E        = [5, 0.24, 4,0.19, 3,4]
MW_E       = [48e3, 65e3, 59.7e3, 70e3, 50e3, 52e3]
D_E        = DEBmicroTrait.aqueous_diffusivity(MW_E)
#
R_P        = 20e-3*ones(n_polymers)
#
α_kin_P    = DEBmicroTrait.binding_sites(ρ_P, R_E)
K_EP_0     = DEBmicroTrait.polymer_affinity(R_E, R_P, V_E, D_E)


#######################################################################
#                     test paper values
# enzyme characteristics
R_E        = [3.29e-9]
V_E        = [5.0]
MW_E       = [48e3]
#D_E        = DEBmicroTrait.aqueous_diffusivity(MW_E)
D_E        = [1e-10]
# Polymer characteristics
ρ_P        = [1.35]
R_P        = [1].*1e-3
#
K_EP_0     = DEBmicroTrait.polymer_affinity(R_E, R_P, V_E, D_E)
@test K_EP_0 ≈ [0.244] atol=0.001
#
k_max_RP   = DEBmicroTrait.max_hydrolyis_rate(R_E, R_P, V_E)
#
R_P        = [1, 10, 20, 30, 40, 50].*1e-3
K_EP_0     = DEBmicroTrait.polymer_affinity(R_E, R_P, V_E, D_E)


V_E        = [42.3, 9.2, 0.41, 34.1, 31.5, 113.5, 128.6, 0.18, 112.2, 0.092, 4.6, 32896, 165, 92.8, 152.13, 69, 0.02, 54, 7590, 16.0, 18.6,
              300.3, 0.39, 0.39]

k_max_RP   = DEBmicroTrait.max_hydrolyis_rate(R_E, R_P, V_E)

kmaxs      = zeros(24,50)
R_P        = collect(range(1,50, length=50))

for i in 1:24
    for j in 1:50
        kmaxs[i,j] = DEBmicroTrait.max_hydrolyis_rate(R_E, [R_P[j]], [V_E[i]])[1]
    end
end

using JLD
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/cellulase_depolymerization.jld", "kmaxs", kmaxs)

using Statistics

median(V_E)
