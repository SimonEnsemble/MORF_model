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
	δ = 1.0
	ϵ♡s = range(-5.0, 0.0, length=10)
	ϵΔs = range(-5.0, 0.0, length=10)
	
	fig, ax = subplots(figsize=(4, 4))
	plot(ϵ♡s, ϵΔs, color=colorz["peace"], lw=3, label="peace")
	plot(ϵ♡s, ϵ♡s .+ δ, color="k", linestyle="--")
	
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
	
	plot([-δ, -δ], [-0.15, 0.15], color="k", clip_on=false)
	text(-δ, 0.4, L"$-\delta$", ha="center", va="center")
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

# ╔═╡ Cell order:
# ╠═86b87b78-fc7f-11ea-2e71-13a4099f7ca2
# ╠═d30e022c-fc7f-11ea-2b57-775cc0121f5e
# ╟─d28bb768-fc7f-11ea-29c2-235136a39b07
# ╠═66d4b94a-fc85-11ea-0654-f72ec2b2e49f
# ╠═ed4675d6-fc7f-11ea-0afa-3de8c3953620
# ╟─596e3acc-fc89-11ea-1fe3-cb9bf6d05345
# ╠═73c59f08-fc89-11ea-153a-7987634c790f
