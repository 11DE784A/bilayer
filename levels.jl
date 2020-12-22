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
# m_e = 9.109e-31u"kg"

# graphene constants
const a = 2.46e-10u"m"				# lattice constant
γ(j) = ([3.16, 0.381]u"eV")[j+1] 	# γ(0), γ(1) parametrize in-plane, dimer hoppings
v(j) = (√3a*γ(j)) / 2ħ 				# v(0), v(1) fermi velocities
const m = γ(1) / (2v(0)^2)			# effective mass of quasiparticles
const w = 0.11eV
const λ = 0.9μ_B
const λ_bar = 1.7μ_B

const dim = 100

# quantum dot constants
const m_e = 9.109e-31u"kg"
const m_dot = 0.067m_e
const l_dot = 50u"nm"
const ω_d = ħ / (m_dot*l_dot^2)

ω_dot(B) = e*B / m_dot
Ω_dot(B) = sqrt(ω_d^2 + ω_dot(B)^2/4)

l(B) = sqrt(ħ/(e*abs(B)))
ω(B) = e*B / m 				 # cyclotron frequency

rad(θ) = ustrip(u"rad", θ) 	 # input θ in degrees, output is 'dimensionless'
k(θ) = (8π/3a)*sin(θ/2)
α(θ) = w / (ħ*v(0)*k(θ))
v_F(θ) = (1 / (1 + 6α(θ)^2)) * v(0) # renormalized fermi velocity

β(T) = (k_B*T)^-1 			 # thermodynamic β

ΔK(θ) = (4π/3) * (rad(θ)/a)
α(θ, B) = (ΔK(θ)*l(B)) / sqrt(8)

plottable(sp) = transpose(hcat(sp...))

δ(i, j) = (i == j ? 1 : 0)
lg(x) = x == 0 ? 0 : log(x)

O = zeros(2, 2)

σ_1, σ_2, σ_3 = [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]
σ = [σ_1, σ_2]

d(σ, k) = sum(σ .* k)

q_1(θ) = k(θ)*[0, 1]
q_2(θ) = k(θ)*[√3/2, 1/2]
q_3(θ) = k(θ)*[-√3/2, 1/2]

T1, T2, T3 = I + σ_1, I - (1/2)σ_1 - (√3/2)σ_2, I - (1/2)σ_1 + (√3/2)σ_2

# Landau Levels

E_mono(B, n) = uconvert(eV, sqrt(2n)*ħ*v(0) / l(B))
sp_mono(B) = [E_mono(B, n) for n in 0:dim]

# semiconductor quantum dot
E_dot(n, m, B) = ħ*Ω_dot(B)*(2n + abs(m) +1) - m*ħ*ω_dot(B)/2
sp_dot(B) = (rdim = floor(sqrt(dim)); 
             reshape([E_dot(n, m, B) for n in 0:rdim, m in -rdim:rdim], :, 1))


sp_no(B) = uconvert(eV, ħ*(e*B/m))*[sqrt(n*(n-1)) for n in 0:dim]

function sp_twist(θ, B)
	e = [1.0, 0, 1, 0, 1, 0, 1]
	up = Bidiagonal(zeros(8), e, :U)
	low = Bidiagonal(zeros(8), e, :L)
	A = Bidiagonal(zeros(dim), [√n for n in 1:(dim-1)], :U)
	
	c = [O*eV   w*T1  w*T2  w*T3; 
		 w*T1'  -ħ*v(0)*d(σ, q_1(θ)) O*eV  O*eV; 
		 w*T2'  O*eV  -ħ*v(0)*d(σ, q_2(θ))  O*eV; 
		 w*T3'  O*eV  O*eV -ħ*v(0)*d(σ, q_3(θ))]
	
	H = ustrip.(eV, (√2*ħ*v(0)/l(B)))*(kron(A, low) + kron(A', up)) + kron(I(dim), ustrip.(eV, c))

    return eigvals(H)*eV
end

#  thermodynamics
P(spec, T) = (p = exp.(-β(T)*spec); p ./ sum(p))
Z(spec, T) = sum(exp.(-β(T)*spec))
U(spec, T) = sum(P(spec, T) .* spec)
S(spec, T) = (p = P(spec, T); -sum(p .* lg.(p)))

function otto(sp, ratio; q = "η")
	B_2 = ratio^2 * B_1
	# B_2 = ratio^2 * (B_S + B_1) - B_S

	sp_1, sp_2 = sp(B_1), sp(B_2)

	T_2, T_4 = 0K, 0K
	try
		T_2 = find_zero(T -> S(sp_1, T_1) - S(sp_2, T), (0.1K, 300K))
		T_4 = find_zero(T -> S(sp_2, T_3) - S(sp_1, T), (0.1K, 300K))
	catch err
		T_2 = find_zero(T -> S(sp_1, T_1) - S(sp_2, T), T_3)
		T_4 = find_zero(T -> S(sp_2, T_3) - S(sp_1, T), T_3)
	end
	T_C, T_H = min(T_1, T_2, T_3, T_4), max(T_1, T_2, T_3, T_4)

	Q_in = abs(U(sp_2, T_3) - U(sp_2, T_2))
	Q_out = abs(U(sp_1, T_1) - U(sp_1, T_4))
	Q_C, Q_H = min(Q_in, Q_out), max(Q_in, Q_out)

	@printf("r_B: %.3f. T_2: %.4f K. T_4: %.4f K. Work: %.4f meV. η: %4f. η_C: %4f. \n", 
			ratio, ustrip(T_2), ustrip(T_4), ustrip(u"meV", abs(Q_H-Q_C)), (1 - Q_C/Q_H), (1 - T_C/T_H))
	# @printf("r_B: %3f. Q_C: %4f meV. Q_H: %4f meV. \n", ratio, ustrip(u"meV", Q_C), ustrip(u"meV", Q_H))

	if q == "η"
		return 1 - Q_C/Q_H, 1 - T_C/T_H
	elseif q == "w"
		return abs(Q_H - Q_C)
	elseif q == "a"
		return 1 - Q_C/Q_H, 1 - T_C/T_H, abs(Q_H - Q_C)
	end
end

function plot_S(T, θ, r, q)
	sp_1, sp_2 = sp_twist(θ, B_1), sp_twist(θ, r^2*B_1)
	if q == "C"
		plot(T, t -> S(sp_1, T_1) - S(sp_2, t))
	elseif q == "H"
		plot(T, t -> S(sp_2, T_3) - S(sp_1, t))
	end
end

function plot_S(T, r, q)
	sp_1, sp_2 = sp_no(B_1), sp_no(r^2*B_1)
	if q == "C"
		plot(T, t -> S(sp_1, T_1) - S(sp_2, t))
	elseif q == "H"
		plot(T, t -> S(sp_2, T_3) - S(sp_1, t))
	end
end

B_1 = 4.0T
B_S = 20.0T
T_1, T_3 = 30.0K, 100.0K

r_B = 1.0:0.01:3.2
η_A = otto.(B -> sp_twist(5.0°, B), r_B)

# Tests. Reproducing plots in Python's thesis.
# B_test = [726.294, 7.26294, 3*0.0726294]*T
# 
# p1 = scatter(sp_twist(1.05°, B_test[1])[4dim:4dim+50], marker=:diamond, legend = false,
# 			 ylims = (-0.1, 4),
# 			 ylabel = "Energy (eV)",
# 			 annotations = (25, 3.5, Plots.text("θ = 1.05°, ω = 1.0", :center)))
# p2 = scatter(sp_twist(1.05°, B_test[2])[4dim:4dim+50], marker=:diamond, legend = false,
# 			 ylims = (-0.01, 0.5),
# 			 ylabel = "Energy (eV)",
# 			 annotations = (25, .45, Plots.text("θ = 1.05°, ω = 0.1", :center)))
# p3 = scatter(sp_twist(1.05°, B_test[3])[4dim:4dim+50], marker=:diamond, legend = false,
# 			 ylims = (-0.001, 0.14),
# 			 ylabel = "Energy (eV)",
# 			 xlabel = "Landau level",
# 			 annotations = (25, 0.12, Plots.text("θ = 1.05°, ω = 0.03", :center)))
# 
# q1 = scatter(sp_twist(5.0°, B_test[1])[4dim:4dim+50], marker=:rect, legend = false,
# 			 ylims = (-0.1, 4),
# 			 annotations = (25, 3.5, Plots.text("θ = 5°, ω = 1.0", :center)))
# q2 = scatter(sp_twist(5.0°, B_test[2])[4dim:4dim+50], marker=:rect, legend = false,
# 			 ylims = (-0.01, 0.5),
# 			 annotations = (25, .45, Plots.text("θ = 5°, ω = 0.1", :center)))
# q3 = scatter(sp_twist(5.0°, B_test[3])[4dim:4dim+50], marker=:rect, legend = false,
# 			 ylims = (-0.001, 0.14),
# 			 xlabel = "Landau level",
# 			 annotations = (25, 0.12, Plots.text("θ = 5°, ω = 0.03", :center)))
# 
# p = plot(p1, q1, p2, q2, p3, q3, layout = (3, 2), markersize = 5)

# Everything is as expected.
