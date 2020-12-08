# units
const eV, T, K, ° =  u"eV", u"T", u"K", u"°"

# universal constants
const e = 1.602e-19u"C" 
const h = 4.136e-15u"eV*s"
const ħ = h/2π
const k_B = 8.617e-5u"eV/K"
const μ_B = 5.788e-5u"eV/T"
const m_e = 9.109e-31u"kg"

# graphene constants
const a = 2.46e-10u"m"               # lattice constant
λ(j) = ([3.16, 0.381]u"eV")[j+1]     # λ(0), λ(1) parametrize in-plane, dimer hoppings
v(j) = (√3a*λ(j)) / 2ħ               # v(0), v(1) fermi velocities
const m = λ(1) / (2v(0)^2)           # effective mass of quasiparticles
const γ = 0.9μ_B                     # Zeeman coupling
const γ̄  = 1.7μ_B                    # pseudo-Zeeman coupling
const A = (1e3a)^2                   # area of sample

# some helper functions
l(B) = sqrt(ħ/(e*abs(B)))            # Landau radius
β(T) = (k_B*T)^-1                    # thermodynamic β
