# using Interpolations
using LinearAlgebra
using LaTeXStrings
using Plots
using Printf
using Roots
using Unitful, UnitfulRecipes

include("PlotSettings.jl")
include("Constants.jl")

const dim = 100

const m_dot = 0.067m_e
const l_dot = 50u"nm"
const ω_d = ħ / (m_dot*l_dot^2)

ω_dot(B) = e*B / m_dot
Ω_dot(B) = sqrt(ω_d^2 + ω_dot(B)^2/4)

Ω(ξ, B, B_S) = v(0) * sqrt(2*e*abs(B + ξ*B_S)/ħ)
Δ(ξ, B, B_S) = (γ̄*B_S) / (ħ*Ω(ξ, B, B_S))
N(ξ, B, B_S) = floor((abs(B + ξ*B_S) * A) / (h/e))

plottable(sp) = transpose(hcat(sp...))

spectrum(E, B) = reshape([E(n, s, ξ, B) 
                          for n in 0:dim, 
                              s in [+1, -1], 
                              ξ in [+1, -1]],
                         :, 1)

# semiconductor quantum dot
E_dot(n, m, B) = ħ*Ω_dot(B)*(2n + abs(m) +1) - m*ħ*ω_dot(B)/2
sp_dot(B) = (rdim = floor(sqrt(dim)); 
             reshape([E_dot(n, m, B) for n in 0:rdim, m in -rdim:rdim], :, 1))

# monolayer graphene
E_mono(n, s, ξ, B) = (sqrt(2n)*ħ*v(0) / l(B))
sp_mono(B) = spectrum(E_mono, B)

# strained monolayer graphene
E_strm(n, s, ξ, B, B_S) = ((n == 0 ?  ξ*ħ*Ω(ξ,B,B_S)*Δ(ξ,B,B_S) : 
                                      ħ*Ω(ξ,B,B_S)*sqrt(n+Δ(ξ,B,B_S)^2)) - s*γ*B)
sp_strm(B, B_S) = spectrum((n, s, ξ, B) -> E_strm(n, s, ξ, B, B_S), B)

# bilayer graphene
E_bi(n, s, ξ, B) = ħ*(e*B/m_eff) * sqrt(n*(n-1))
sp_bi(B) = spectrum(E_bi, B)

# twisted bilayer graphene
h_m(k) = -ħ*v(0)*[0u"m^-1"  k[1] + im*k[2];
                  k[1] - im*k[2]  0u"m^-1"]

function sp_twist(θ, B, positive = true)
    Π = Bidiagonal(zeros(dim+1), sqrt.(1:dim), :U)
    U = Bidiagonal(zeros(8), [1., 0, 1, 0, 1, 0, 1], :U)
    L = U'

    H_c = ustrip.(eV, [O*eV     w*T_b        w*T_tr        w*T_tl;
                       w*T_b'   h_m(q_b(θ))  O*eV          O*eV;
                       w*T_tr'  O*eV         h_m(q_tr(θ))  O*eV;
                       w*T_tl'  O*eV         O*eV  h_m(q_tl(θ))])

    H_k = (√2*ħ*v(0)/l(B)) * (kron(U, Π') + kron(L, Π))

    H_b = ustrip.(eV, H_k) + kron(H_c, I(dim+1))

    if positive
        return eigvals(H_b)[4*(dim+1)+1:end]*eV
    else
        return eigvals(H_b)*eV
    end
end

# Plotting the spectrum
# B = (0.5:0.05:15)*T                 # magnetic field range
# θ = 1.1°                            # twist angle
# twi = plottable(sp_twist.(θ, B))    # computing the energy levels
# p = plot(B, twi, title = "θ = $(θ)°, N = $(dim)", ylims = (-0.2, 0.2),
#          legend = false, color = :red)
# savefig(p, "plot.png")
