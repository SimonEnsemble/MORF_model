### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 86b87b78-fc7f-11ea-2e71-13a4099f7ca2
begin
	using PyPlot, Printf, SymPy, ForwardDiff, Roots
	
	# plotting business
	PyPlot.matplotlib.font_manager.fontManager.addfont("OpenSans-Regular.ttf")
	PyPlot.matplotlib.style.use("grandbudapest.mplstyle")
	
	function draw_axes(;color="0.6")
	    axvline(x=0, color=color, lw=2, zorder=1)
	    axhline(y=0, color=color, lw=2, zorder=1)
	end
	
	colorz = Dict("no gas"      => "k", 
		          "peace"       => "C2", 
		          "cooperation" => "C0", 
		          "competition" => "C3")
	marker = "X"
end

# ╔═╡ d30e022c-fc7f-11ea-2b57-775cc0121f5e
const kT_room = 2.479 # kJ/mol

# ╔═╡ d28bb768-fc7f-11ea-29c2-235136a39b07
md"
# the entropy of the wheel
"

# ╔═╡ 66d4b94a-fc85-11ea-0654-f72ec2b2e49f
begin
	# θ = ϵ♡ + δ - ϵ△
	w_n1(βθ) = exp(-βθ) / (1.0 + exp(-βθ))
	S(w) = -(w * log(w) + (1.0 - w) * log.(1.0 - w))
end

# ╔═╡ ed4675d6-fc7f-11ea-0afa-3de8c3953620
begin
	# θ = ϵ♡ + δ - ϵ△
	βθs = range(-6.5, 6.5, length=100)
	
	ht = 0.8          # height
	wd = βθs[end] * 2 # width
	function plot_θ(βθ::Float64, marker::String, a_color, label, dx::Float64, dy::Float64)
	    scatter([βθ], [S(w_n1(βθ))], marker=marker, zorder=4, 
			s=marker=="*" ? 1.35 * 150 : 150, color=a_color)
	    if label == "none"
	        return
	    end
	    if label != ""
	        annotate(label,
	                 xy=(βθ, S(w_n1(βθ))), xycoords="data",
	                 xytext=(βθ + dx * wd, S(w_n1(βθ)) + dy * ht), textcoords="data",
	                 arrowprops=Dict(:arrowstyle => "->", #linestyle="dashed",
	                                 :color => "k"
	                                ),
	                 zorder=1241
	                 )
	    end
	end
	
	figure(figsize=(7, 4))
	
	###
	#   fixed n = 0 dot
	###
	βδ = 3.0
# 	plot_θ(βδ, "o", colorz["no gas"], L"fixed $n=0$", -0.03, 0.12)
	scatter([βδ], [S(w_n1(βδ))], marker="o", zorder=4, 
			s=150, color="k", label=L"fixed $n=0$")
	
	###
	#   fixed n = 1 curve
	###
	plot(βθs, S.(w_n1.(βθs)), zorder=2, color="0.5", lw=3, label=L"fixed $n=1$")
	ylabel(L"wheel entropy, $S_w/k_B$")
	xlabel(L"$\beta(\epsilon_\heartsuit+\delta -\epsilon_\bigtriangleup)$")
	legend(loc="upper left")
	
	###
	#   draw boxes for peace, competition, and cooperation
	###
	dx = 0.4
	bbox_props = Dict{Symbol, Any}(:lw=>2, :alpha=>0.2)
	
	bbox_props[:fc] = colorz["cooperation"]
	bbox_props[:ec] = colorz["cooperation"]
	bbox_props[:boxstyle] = "rarrow,pad=0.2"
	text([βδ + dx], [S(w_n1(βδ))], "wheel-gas\ncooperation ", 
		ha="left", va="center", rotation=0,
	    size=14, bbox=bbox_props)
	
	bbox_props[:fc] = colorz["competition"]
	bbox_props[:ec] = colorz["competition"]
	bbox_props[:boxstyle] = "larrow,pad=0.2"
	text([βδ - dx], [S(w_n1(βδ))], "wheel-gas\ncompetition", 
		ha="right", va="center", rotation=0,
	    size=14, bbox=bbox_props)
	
	bbox_props[:fc] = colorz["peace"]
	bbox_props[:ec] = colorz["peace"]
	bbox_props[:boxstyle] = "Round,pad=0.2"
	text(βδ, 0.3, "wheel-gas" * "\n" * "peace", 
		 ha="center", bbox=bbox_props)

	###
	#  fixed n = 1 points
	###
	# peace
	plot_θ(βδ, marker, colorz["peace"], "", 0.0, 0.0)
	
	# with gas. wheel-gas peace
	βϵ♡ = -1.0
	βϵ△ = -3.0
	θ = βϵ♡+βδ-βϵ△
	plot_θ(θ , marker, colorz["cooperation"], "none", 0.1/sqrt(2), 0.1/sqrt(2))
	
	# with gas. wheel-gas competition. wheel wins.
	βϵ♡ = -2.5
	βϵ△ = -1.0
	θ = βϵ♡+βδ-βϵ△
	plot_θ(θ , marker, colorz["competition"], "wheel wins", 0.0, 0.1)
	
	# with gas. wheel-gas competition. detente
	βϵ♡ = -4.0
	βϵ△ = -1.0
	θ = βϵ♡+βδ-βϵ△
	plot_θ(θ , marker, colorz["competition"], "detente", 0.1/sqrt(2), 0.1/sqrt(2))
	
	# with gas. wheel-gas competition. gas wins.
	βϵ♡ = -5.5
	βϵ△ = -1.0
	θ = βϵ♡+βδ-βϵ△
	plot_θ(θ , marker, colorz["competition"], "gas wins", -0.2, 0.1)
	
	# with gas. wheel-gas competition. gas wins by far.
	βϵ♡ = -9.0
	βϵ△ = -1.0
	θ = βϵ♡+βδ-βϵ△
	plot_θ(θ , marker, colorz["competition"], "gas wins by far", -0.135, 0.125)
	
	# # legend(numpoints=1, bbox_to_anchor=(1, 0.8))
	# title(L"how does adsorbed gas affect $S_\mathrm{wheel}$?")
	# # draw_axes()
	xlim([-wd/2, wd/2])
	ylim([-0.01, ht])
	
	fill_between([-βδ, βδ], zeros(2), ht * ones(2), color="0.5", alpha=0.1, zorder=0)
	
	vlines(βδ, 0.0, S(w_n1(βδ)), linestyle="--", color="0.5")
	vlines(-βδ, 0.0, S(w_n1(-βδ)), linestyle="--", color="0.5")
	x_t = collect(-6:2:6)
	xticks(x_t)
	text(βδ, -0.05, L"$\beta\delta$", ha="center")
	text(-βδ, -0.05, L"$-\beta\delta$", ha="center")
	
	savefig("wheel_entropy.pdf", format="pdf", bbox_inches="tight")
	
	gcf()
end

# ╔═╡ 596e3acc-fc89-11ea-1fe3-cb9bf6d05345
md"
# materials space
"

# ╔═╡ 73c59f08-fc89-11ea-153a-7987634c790f
begin
	δ_ms = 1.0
	ϵ♡s = range(-5.0, 0.0, length=10)
	ϵΔs = range(-5.0, 0.0, length=10)
	
	fig, ax = subplots(figsize=(4, 4))
	plot(ϵ♡s, ϵΔs, color=colorz["peace"], lw=3, label="peace")
	plot(ϵ♡s, ϵ♡s .+ δ_ms, color="k", linestyle="--")
	
	xticks([])
	yticks([])
	
	fill_between(ϵ♡s, -5.0, ϵ♡s, color=colorz["cooperation"], label="cooperation")
	fill_between(ϵ♡s, 0.0, ϵ♡s, color=colorz["competition"], label="competition")
	xlim([-5.4, 0.1])
	ylim([-5.4, 0.1])
	
	props = Dict("boxstyle"=>"round", "facecolor"=>"white", "alpha"=>0.5)
	text(-3.8, -0.75, "competition", ha="center", va="center", bbox=props)
	text(-1.25, -3.8, "cooperation", ha="center", va="center", bbox=props)
	text(-4, -4, "peace", ha="center", va="center", color="k", bbox=props, rotation=45)
	
	text(-2.5, -2.0, "wheel wins", ha="center", va="center", rotation=45.0, color="k")
	text(-3.0, -1.5, "gas wins", ha="center", va="center", rotation=45.0, color="k")
	
	plot([-δ_ms, -δ_ms], [-0.15, 0.15], color="k", clip_on=false)
	text(-δ_ms, 0.4, L"$-\delta$", ha="center", va="center")
	text(0.0, 0.4, L"$0$", ha="center", va="center")
	text(0.4, 0.0, L"$0$", ha="center", va="center")
	text(0.0, -5.5, L"$\epsilon_\bigtriangleup$", ha="center", va="center", fontsize=18)
	text(-5.5, 0.0, L"$\epsilon_\heartsuit$", ha="center", va="center", fontsize=18)
	ax.annotate("", xy=(-5.2, 0.0), xytext=(0.2, 0.0),
	            arrowprops=Dict("facecolor"=>"k", "lw"=>1))
	ax.annotate("", xy=(0.0, -5.2), xytext=(0.0, 0.2),
	            arrowprops=Dict("facecolor"=>"k", "lw"=>1))
	
	ax.set_aspect("equal", adjustable="box")
	tight_layout()
	
	savefig("mat_space.pdf", format="pdf")
	
	gcf()
end

# ╔═╡ 4b8617e4-fc8a-11ea-0cd1-2ffba5c42303
md"
# Langmuir isotherm
"

# ╔═╡ 517c173e-fc8a-11ea-28e6-ab0c43713125
begin
	kp = range(0.0, 10.0, length=200)
	exp_n = kp ./ (1 .+ kp)
	
	figure(figsize=(6.4*.85, 4.8*.85))
	plot(kp, exp_n, linewidth=4, color="C4")
	ϵ=1e-1
	xlim([0-ϵ, maximum(kp)+ϵ])
	ϵ=1e-2
	ylim([0-ϵ, 1+ϵ])
	xlabel(L"$K\beta P$")
	ylabel(L"$\langle n \rangle$")
	draw_axes()
	text(6, 0.4, L"$\langle n \rangle=\dfrac{K\beta P}{1+K\beta P}$", fontsize=18, 
		bbox=Dict("boxstyle"=>"round", "facecolor"=>"wheat", "alpha"=>0.3))
	vlines(1, 0, 0.5, color="0.7", linestyle="--")
	hlines(0.5, 0, 1, color="0.7", linestyle="--")
	# title("Langmuir adsorption")
	tight_layout()
	savefig("langmuir.pdf", format="pdf")
	
	gcf()
end

# ╔═╡ c5d406d2-fc8a-11ea-378d-7388e9ad53d4
md"
# RMS-MOF analysis
"

# ╔═╡ dfdda3b0-fc8a-11ea-3777-c7fc194bd218
begin
	struct Material
	    δ::Float64
	    ϵ♡::Float64
	    ϵ△::Float64
	end
	
	function Base.show(io::IO, material::Material)
	    @printf("\tδ = %.3f\n\tϵ♡ = %.3f\n\tϵ△ = %.3f\n", 
			    material.δ, material.ϵ♡, material.ϵ△)
	end
end

# ╔═╡ fa8835e0-fc8a-11ea-00cf-25aa4e5c3ffb
begin
	# kT::Real as opposed to kT::Float64 for ForwardDiff.jl to work
	
	function K(material::Material, kT::Real) 
	    p♡_n0 = 1 / (1 + exp(-material.δ / kT))
	    return p♡_n0 * exp(-material.ϵ△ / kT) + (1 - p♡_n0) * exp(-material.ϵ♡ / kT)
	end
	
	function n(material::Material, kT::Real, P::Float64)
	    KβP = K(material, kT) * P / kT
	    return KβP / (1.0 + KβP)
	end
	
	function w(material::Material, kT::Real, P::Float64)
		pΔ_n0 = exp(-material.δ / kT) / (1 + exp(-material.δ / kT))
		K♡βP = exp(-material.ϵ♡ / kT) * P / kT
		KβP = K(material, kT) * P / kT
	    return pΔ_n0  * (1 + K♡βP) / (1 + KβP)
	end
	
	function ∂E_∂n(material::Material, kT::Float64)
	    β = 1 / (kT)
	    w_n0 = exp(-β * material.δ) / (1 + exp(-β * material.δ))
	    w_n1 = exp(-β * (material.ϵ♡ + material.δ)) / (exp(-β * material.ϵ△) + exp(-β * (material.ϵ♡ + material.δ)))
	    return material.δ * (w_n1 - w_n0) + material.ϵ♡ * w_n1 + material.ϵ△ * (1 - w_n1)
	end
	
	function ∂E□_∂n□(ϵ□::Float64)
	    return ϵ□
	end
	
	@assert ∂E□_∂n□(-4.0) ≈ ∂E_∂n(Material(0.0, -4.0, -4.0), rand())
	@assert ∂E□_∂n□(-4.0) ≈ ∂E_∂n(Material(2.0, -4.0, -4.0), rand())
	@assert ∂E□_∂n□(-2.0) ≈ ∂E_∂n(Material(2.0, -4.0, 0.0), rand()) # weird case
end

# ╔═╡ 409daf48-fc8c-11ea-0811-d3bbaa86a15b
md"
use `ForwardDiff.jl` to autodiff
"

# ╔═╡ 5b5c269a-fc8c-11ea-3eee-c3767d47d678
begin
	function ∂n_∂kT(material::Material, kT::Float64, P::Float64)
		# view n as a function of kT only.
		n_of_kT(x) = n(material, x[1], P) # x = [kT]
		∂n_∂kT_of_T = x -> ForwardDiff.gradient(n_of_kT, x) # x plays role of kT
		return ∂n_∂kT_of_T([kT])[1]
	end
	
	function ∂K_∂kT(material::Material, kT::Float64)
		# view K as a function of kT only.
		K_of_kT(x) = K(material, x[1]) # x = [kT]
		∂K_∂kT_of_kT = x -> ForwardDiff.gradient(K_of_kT, x) # x plays role of kT
		return ∂K_∂kT_of_kT([kT])[1]
	end
end

# ╔═╡ 8f1e2f58-fc8c-11ea-1fde-9d1d72b5d631
md"
## check these with symbolic differentiation in `SymPy`
"

# ╔═╡ 9a054b92-fc8c-11ea-0aaf-f90ed61e9586
begin
	# define symbol versions of the variables
	_kT = Sym("kT")
	_ϵ♡ = Sym("ϵ♡")
	_ϵ△ = Sym("ϵ△")
	_δ = Sym("δ")
	_P = Sym("P")
end

# ╔═╡ 323ced64-fc8d-11ea-0ffe-7d5a015476ed
K_symbolic = 1 / (1+exp(-_δ/_kT)) * (exp(-_ϵ△/_kT) + exp(-(_ϵ♡+_δ)/_kT))

# ╔═╡ 2a47bede-fc8d-11ea-2234-517b57146be1
n_symbolic = K_symbolic * _P / _kT / (1 + K_symbolic * _P / _kT)

# ╔═╡ dc88b1ec-fc8d-11ea-0de4-f51f94533fea
w_symbolic = exp(-_δ / _kT) / (1 + exp(-_δ / _kT)) * (1 + exp(-_ϵ♡ / _kT) * _P / _kT) / (1 + K_symbolic * _P / _kT)

# ╔═╡ 45f2e688-fc8d-11ea-169a-fd11486da805
begin
	# too complicated to parse anyway
	∂n_∂kT_symbolic = diff(n_symbolic, "kT")
	∂K_∂kT_symbolic = diff(K_symbolic, "kT")
end

# ╔═╡ 6b68afd0-fc8d-11ea-30f2-cd181b6d5fcf
md"
to evaluate these symbolic expressions at numerical values
"

# ╔═╡ 5719a1d6-fc8d-11ea-000b-ddcb985f72f0
begin
	# if P is not involved
	function eval_sym(sym::Sym, material::Material, kT::Float64)
	    params = Dict("kT" => kT, "ϵ♡" => material.ϵ♡, "ϵ△" => material.ϵ△, "δ" => material.δ)
	    return convert(Float64, sym.subs(params))
	end
	
	# if P is involved
	function eval_sym(sym::Sym, material::Material, kT::Float64, P::Float64)
	    params = Dict("kT" => kT, "ϵ♡" => material.ϵ♡, "ϵ△" => material.ϵ△, "δ" => material.δ, "P" => P)
	    return convert(Float64, sym.subs(params))
	end
	
	function test_symbolic_expressions()
		material = Material(rand(), -rand(), -rand())
		kT = rand()
		P = rand()

		# K
		@assert isapprox(K(material, kT), eval_sym(K_symbolic, material, kT))
		# dK/d(kT)
		@assert isapprox(∂K_∂kT(material, kT), eval_sym(∂K_∂kT_symbolic, material, kT))
		# n (P=0.2)
		@assert isapprox(n(material, kT, P), eval_sym(n_symbolic, material, kT, P))
		# dn/d(kT)
		@assert isapprox(∂n_∂kT(material, kT, P), eval_sym(∂n_∂kT_symbolic, material, kT, P))
		
		# w
		@assert isapprox(w(material, kT, P), eval_sym(w_symbolic, material, kT, P))
		"tests pass"
	end
	
	test_symbolic_expressions()
end

# ╔═╡ 551fba04-fc8e-11ea-2058-2bda590970f8
md"
## cognate Langmuir material comparison
"

# ╔═╡ 4658aca4-fc8e-11ea-04b5-b78a558415a3
begin
	"""
	What is the energy of adsorption ϵ□ for a Langmuir material 
	that exhibits the same K as the RMS-MOF at temperature kT?
	K□ = e^(-ϵ□/kT)
	set K = K□ => ϵ□ = - kT * log(K)
	"""
	get_ϵ□(material::Material, kT::Real) = -kT * log(K(material, kT))
	
	function ∂K□_∂kT(ϵ□::Float64, kT::Float64)
	    K□ = exp(-ϵ□ / kT)
	    return K□ * ϵ□ / kT ^ 2
	end
	
	function n□(ϵ□::Float64, kT::Real, P::Float64)
	    K□ = exp(-ϵ□ / kT)
	    KβP = K□ * P / kT
	    return KβP / (1 + KβP)
	end
	
	# ! IMPORTANT! keep ϵ□ fixed here.
	#    i.e. don't re-set ϵ□...
	function ∂n□_∂kT(ϵ□::Float64, kT::Real, P::Float64)
	    # view n as a function of kT only.
	    n□_of_kT(x) = n□(ϵ□, x[1], P) # x = [kT]
	    ∂n□_∂kT_of_T = x -> ForwardDiff.gradient(n□_of_kT, x) # x plays role of kT
	    return ∂n□_∂kT_of_T([kT])[1]
	end
	
	function test_□_stuff()
		ϵ□_test = randn()
		RMSMOF□□ = Material(randn(), ϵ□_test, ϵ□_test)
		RMSMOFrand = Material(rand(), randn(), randn())
		kT = rand()
		P = rand()
		
		@assert exp(-get_ϵ□(RMSMOFrand, kT_room) / kT_room) ≈ K(RMSMOFrand, kT_room)
		@assert get_ϵ□(RMSMOF□□, kT_room) ≈ ϵ□_test
		@assert n□(ϵ□_test, kT, P) ≈ n(RMSMOF□□, kT, P)
		"tests pass"
	end
	
	test_□_stuff()
end

# ╔═╡ 20c6b734-fc8f-11ea-3244-55e8ff942305
begin
	δ = 3.0
	kT = kT_room
	nb_pts = 50
	ϵ_range = range(-10.0, stop=0.0, length=nb_pts)
	
	ϵ□s = zeros(nb_pts, nb_pts)
	∂K_∂T_diffs = zeros(nb_pts, nb_pts)
	∂E_∂N_diffs = zeros(nb_pts, nb_pts)
	for (i, ϵ♡) in enumerate(ϵ_range)
	    for (j, ϵ△) in enumerate(ϵ_range)
	        material = Material(δ, ϵ♡, ϵ△)
	        
	        # find epsilon that gives same K for Langmuir material
	        ϵ□ = get_ϵ□(material, kT)
	        K□ = exp(-ϵ□ / kT)
	        @assert K□ ≈ K(material, kT)
	        
	        ϵ□s[j, i] = ϵ□
	        ∂E_∂N_diffs[j, i] = ∂E_∂n(material, kT) - ∂E□_∂n□(ϵ□)
	        ∂K_∂T_diffs[j, i] = (∂K_∂kT(material, kT) - ∂K□_∂kT(ϵ□, kT)) / K□
			
	        @assert ∂K_∂kT(material, kT) <= 0.0
	        @assert ∂K□_∂kT(ϵ□, kT) <= 0.0
	#         z[j, i] = ∂K_∂kT(material, kT) - ∂K□_∂kT(ϵ□, kT)
	        # j, i... this is not a bug.
	        #  try the following two to see:
	        #    Z[j, i] = βϵ♡
	        #    Z[j, i] = βϵ△
	    end
	end
	
	function decorate_fig(cbar_label::LaTeXString)
	    gca().set_aspect("equal", "box")
	
	    cbar = colorbar(label=cbar_label)
	    xlabel(L"$\epsilon_\heartsuit$ [kJ/mol]")
	    ylabel(L"$\epsilon_\bigtriangleup$ [kJ/mol]")
	    ylim(ymax=0.0)
	    tight_layout()
	end
end

# ╔═╡ 7f6dbbdc-fc91-11ea-1212-63b9919b72f5
begin
	figure()
	pcolor(ϵ_range, ϵ_range, ϵ□s, vmax=0.0, cmap=plt.cm.turbo)
	decorate_fig(L"$\epsilon_□$ [kJ/mol]")
	text(-9.5, -1.1, 
		L"$k_BT_0=$" * @sprintf("%.2f kJ/mol", kT) * "\n" * 
		L"$\delta=$" * @sprintf("%.2f kJ/mol", δ), ha="left", va="center")
	tight_layout()
	savefig("epsilon_square.pdf", format="pdf")
	gcf()
end

# ╔═╡ 392ad090-fc93-11ea-109b-5772e9d1ead3
begin
	props□ = Dict("boxstyle"=>"round", "facecolor"=>"white", "alpha"=>0.85)

	# dK/dT
	figure()
	pcolor(ϵ_range, ϵ_range, ∂K_∂T_diffs, cmap=plt.cm.PiYG, 
		vmax=maximum(abs.(∂K_∂T_diffs)), vmin=-maximum(abs.(∂K_∂T_diffs)))
	plot(ϵ_range, ϵ_range .+ 2 * δ, linestyle="--", color="0.5")
	plot(ϵ_range, ϵ_range, linestyle="--", color="0.5")
	text(-5, -5, L"$\epsilon_\bigtriangleup=\epsilon_\heartsuit$", 
		ha="center", va="center", rotation=45, bbox=props□)
	text(-8, -8+2*δ, L"$\epsilon_\bigtriangleup=\epsilon_\heartsuit+2\delta$",
		ha="center", va="center", rotation=45, bbox=props□)
	text(-5, -9, L"$k_BT_0=$" * @sprintf("%.2f kJ/mol", kT) * "\n" * 
		L"$\delta=$" * @sprintf("%.2f kJ/mol", δ), ha="left", va="center")
	decorate_fig(L"$\frac{1}{K}\frac{d K}{d (k_BT)}|_{T_0}-\frac{1}{K_□}\frac{d K_□}{d (k_BT)}|_{T_0}$ [1/(kJ/mol)]")
	tight_layout()
	savefig("dK_dT.pdf", format="pdf")
	gcf()
end

# ╔═╡ 8febf99a-fc93-11ea-2437-434f7dd1b923
begin
	figure()
	pcolor(ϵ_range, ϵ_range, ∂E_∂N_diffs, cmap=plt.cm.bwr_r, 
		vmax=maximum(abs.(∂E_∂N_diffs)), vmin=-maximum(abs.(∂E_∂N_diffs)))
	plot(ϵ_range, ϵ_range .+ 2 * δ, linestyle="--", color="0.5")
	plot(ϵ_range, ϵ_range, linestyle="--", color="0.5")
	text(-5, -5, L"$\epsilon_\bigtriangleup=\epsilon_\heartsuit$", ha="center", va="center", rotation=45, bbox=props)
	text(-8, -8+2*δ, L"$\epsilon_\bigtriangleup=\epsilon_\heartsuit+2\delta$", ha="center", va="center", rotation=45, bbox=props)
	text(-5, -9, L"$k_BT_0=$" * @sprintf("%.2f kJ/mol", kT) * "\n" * L"$\delta=$" * @sprintf("%.2f kJ/mol", δ), ha="left", va="center")
	decorate_fig(L"$\frac{\partial \langle E\rangle }{\partial \langle n \rangle}-\frac{\partial \langle E\rangle_L }{\partial \langle n \rangle_L}$ [kJ/mol]")
	tight_layout()
	savefig("dE_dNs.pdf", format="pdf")
	gcf()
end

# ╔═╡ Cell order:
# ╠═86b87b78-fc7f-11ea-2e71-13a4099f7ca2
# ╠═d30e022c-fc7f-11ea-2b57-775cc0121f5e
# ╟─d28bb768-fc7f-11ea-29c2-235136a39b07
# ╠═66d4b94a-fc85-11ea-0654-f72ec2b2e49f
# ╠═ed4675d6-fc7f-11ea-0afa-3de8c3953620
# ╟─596e3acc-fc89-11ea-1fe3-cb9bf6d05345
# ╠═73c59f08-fc89-11ea-153a-7987634c790f
# ╟─4b8617e4-fc8a-11ea-0cd1-2ffba5c42303
# ╠═517c173e-fc8a-11ea-28e6-ab0c43713125
# ╟─c5d406d2-fc8a-11ea-378d-7388e9ad53d4
# ╠═dfdda3b0-fc8a-11ea-3777-c7fc194bd218
# ╠═fa8835e0-fc8a-11ea-00cf-25aa4e5c3ffb
# ╟─409daf48-fc8c-11ea-0811-d3bbaa86a15b
# ╠═5b5c269a-fc8c-11ea-3eee-c3767d47d678
# ╟─8f1e2f58-fc8c-11ea-1fde-9d1d72b5d631
# ╠═9a054b92-fc8c-11ea-0aaf-f90ed61e9586
# ╠═323ced64-fc8d-11ea-0ffe-7d5a015476ed
# ╠═2a47bede-fc8d-11ea-2234-517b57146be1
# ╠═dc88b1ec-fc8d-11ea-0de4-f51f94533fea
# ╠═45f2e688-fc8d-11ea-169a-fd11486da805
# ╟─6b68afd0-fc8d-11ea-30f2-cd181b6d5fcf
# ╠═5719a1d6-fc8d-11ea-000b-ddcb985f72f0
# ╟─551fba04-fc8e-11ea-2058-2bda590970f8
# ╠═4658aca4-fc8e-11ea-04b5-b78a558415a3
# ╠═20c6b734-fc8f-11ea-3244-55e8ff942305
# ╠═7f6dbbdc-fc91-11ea-1212-63b9919b72f5
# ╠═392ad090-fc93-11ea-109b-5772e9d1ead3
# ╠═8febf99a-fc93-11ea-2437-434f7dd1b923
