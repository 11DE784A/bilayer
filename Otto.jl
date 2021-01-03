include("Spectra.jl")

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

B_1 = 5.0T
T_1, T_3 = 30.0K, 100.0K

r = 1.2
θ = 1.0°
otto_res = otto(B -> sp_twist(θ, B), r, q = "a")

