using Printf
using Unitful, UnitfulRecipes
using Plots
using Roots
using LinearAlgebra
using LaTeXStrings

# upscaled plot
upscale = 2.4
fntsm = Plots.font(pointsize=round(5.5*upscale))
fntlg = Plots.font(pointsize=round(6.0*upscale))
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm, 
        legendtitlefont=fntlg)
default(size=(400*upscale,300*upscale)) #Plot canvas size
default(linewidth=upscale)

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

l(B) = sqrt(ħ/(e*abs(B)))            # Landau radius
β(T) = (k_B*T)^-1                    # thermodynamic β

Ω(ξ, B, B_S) = v(0) * sqrt(2*e*abs(B + ξ*B_S)/ħ)
Δ(ξ, B, B_S) = (γ̄*B_S) / (ħ*Ω(ξ, B, B_S))
N(ξ, B, B_S) = floor((abs(B + ξ*B_S) * A) / (h/e))

function E(T, B)
    # From Equation B6
    return (-γ*B*tanh(β(T)*γ*B) 
            + (sum([N(ξ, B, B_S) * (ξ*γ̄*B_S*exp(-ξ*β(T)*γ̄*B_S) 
                                    + (ħ*Ω(ξ, B, B_S)*sqrt(1 + Δ(ξ, B, B_S)^2) 
                                       * exp(-β(T)*ħ*Ω(ξ, B, B_S)*sqrt(1 + Δ(ξ, B, B_S)^2)))) 
                    for ξ in [-1, 1]]) 
               / sum([N(ξ, B, B_S) * (exp(-ξ*β(T)*γ̄*B_S) 
                                      + exp(-β(T)*ħ*Ω(ξ, B, B_S)*sqrt(1 + Δ(ξ, B, B_S)^2))) 
                      for ξ in [-1, 1]])))
end

function S(T, B)
    # From Equation B7
    return (β(T)*E(T, B) 
            + log(2*cosh(β(T)*γ*B)) 
            + log(sum([N(ξ, B, B_S) * (exp(-ξ*β(T)*γ̄*B_S) 
                                       + exp(-β(T)*ħ*Ω(ξ, B, B_S)*sqrt(1 + Δ(ξ, B, B_S)^2)))
                      for ξ in [-1, 1]])))
end

function otto(r)
    B_2 = r^2 * B_1

    # Intermediate temperatures determined numerically. Equation 36.
    T_2 = find_zero(T -> S(T_C, B_1) - S(T, B_2), T_C)
    T_4 = find_zero(T -> S(T_H, B_2) - S(T, B_1), T_H)

    Q_H = E(T_H, B_2) - E(T_2, B_2)   # Equation 33
    Q_C = E(T_C, B_1) - E(T_4, B_1)   # Equation 34
    W = Q_C + Q_H                     # Equations 31 and 32
    η = 1 - abs(Q_C/Q_H)              # Equation 35

    @printf("r_B: %.3f. T_2: %.4f K. T_4: %.4f K. W: %4f meV. η: %4f.\n", 
            r, ustrip(K, T_2), ustrip(K, T_4), ustrip(u"meV", W), η)
end

B_1, B_S = 4.0T, 20.0T
T_C, T_H = 30.0K, 100.0K
