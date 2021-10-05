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

@test 1==1
