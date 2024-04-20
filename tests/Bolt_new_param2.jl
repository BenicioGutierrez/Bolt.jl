### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 12da0ac6-8a54-40f6-92a5-bde72b5e4370
import Pkg

# ╔═╡ 1b4aefd2-77e9-4cee-999e-35b0513ded02
begin
	Pkg.add(url="https://github.com/BenicioGutierrez/Bolt.jl")
	Pkg.add("ForwardDiff")
	Pkg.add("Plots")
	Pkg.add("DataInterpolations")
end

# ╔═╡ fd436b90-809f-4b6d-9b88-f0201bde111d
using Bolt, Plots, ForwardDiff, LaTeXStrings, ThreadTools, Base.Threads, DataInterpolations, Interpolations

# ╔═╡ a48c8d20-febe-11ee-011b-8b62bce5a872
md"""
# Adding a new parameter (Part 2)
"""

# ╔═╡ 13fe63fb-2335-47b7-b0c1-b8cafcb8f96c
md"""
Here, we run some tests to observe how the new energy density affects results
"""

# ╔═╡ 8f150cf5-0cda-4cd3-b77b-f89cd892ead3
md"""
Starting with the base case, where Ω_new = 0 (as before) :
"""

# ╔═╡ 7c11853a-5e2b-432a-a296-6c3bba21101c
# Assign cosmological parameters

𝕡 = CosmoParams(Ω_c = 0.26)

# ╔═╡ 2ec9bcd4-b0c1-434c-a95d-145ae73c5c66
function FRW_setup(𝕡)
    # Compute expansion history quantities
    bg = Background(𝕡)
    # Compute ionization history (via RECFAST)
#    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
	𝕣 = Bolt.RECFAST(bg; Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b,OmegaG=𝕡.Ω_r)  # Check RECFAST
    ih = IonizationHistory(𝕣, 𝕡, bg)
    return bg, ih
end

# ╔═╡ 74a2f5b5-46bb-455b-84fb-ad949964d704
begin
	a_min = 10.0^(-11)
	a_max = 1.00
	n_a = 100
	a = log_a(a_min, a_max, n_a) # a grid
end

# ╔═╡ 4086fc84-e831-4ef0-8a0e-6bfe3eee266c
bg, ih = FRW_setup(𝕡);

# ╔═╡ ce21d45a-8804-4b2e-92a1-de184f9279ae
# Matter power spectrum

begin
	kmin1 = 1bg.H₀
	kmax1 = 5000bg.H₀
	nk1 = 32
	ks1 = log10_k(kmin1, kmax1, nk1) # k grid
end

# ╔═╡ e32b2579-a2f5-443f-8aa2-0c1a6f215678
pL = tmap(k -> plin(k, 𝕡, bg, ih), ks1)

# ╔═╡ b0928a06-5c68-4c58-91dd-d3bd1468d060
# CMB Cᵀᵀ(ℓ)
begin 
	ℓmin, ℓmax, nℓ = 2, 20, 1200
	ℓs = ℓmin:ℓmax:nℓ
	kmin, kmax, nk = 0.1bg.H₀, 1000bg.H₀, 100
	ks = quadratic_k(kmin, kmax, nk)
	sf = source_grid(𝕡, bg, ih, ks, BasicNewtonian()) # set up LOS source function interpolator
end

# ╔═╡ 93a7acd7-37c9-4879-ae3c-7d449e793fe1
Cᵀᵀ = cltt(ℓs, 𝕡, bg, ih, sf)

# ╔═╡ 530cf8b5-e7e5-4b05-b69b-c52e144ff92e
md"""
## Test 4: Ω_new as a spline
"""

# ╔═╡ 4c2bd4a5-ba03-4cf4-aebe-d10c53c2a383
md"""
Ω_new is modeled as a gaussian:
Peak: 5
x_peak: a = 0.000833

The standard deviation will be varied
"""

# ╔═╡ 44b6ddcb-53a5-42ca-a589-254e788507c9
begin
	gauss(x::Float64, peak, x_peak, std) = peak * exp(-(x - x_peak)^2 / (2 * std^2))

	x_peak = 0.000833
	peak = 5.00
end

# ╔═╡ bf328343-ae35-42be-b241-02cf144162eb
begin
	
a_grid = 1e-11:1e-6:1.0
#a_grid = [10^x for x in a_log]

#a_start = log10(1e-15)
#a_end = log10(1)
#num_points = 500

#a_log = range(a_start, a_end, length=num_points)

#a_grid = 10 .^ a_log


end

# ╔═╡ 4596a9c2-1046-45bb-8b28-394fe3e5ddc5
begin

	Ω_new_1_nodes = [gauss(a, peak, x_peak, 1e-4) for a in a_grid]
	Ω_new_1 = cubic_spline_interpolation(a_grid, Ω_new_1_nodes, extrapolation_bc=Flat())

	Ω_new_2 = cubic_spline_interpolation(a_grid, 
		[gauss(x, peak, x_peak, 1e-3) for x in a_grid], extrapolation_bc=Flat())

	Ω_new_3 = cubic_spline_interpolation(a_grid, 
		[gauss(x, peak, x_peak, 1e-2) for x in a_grid], extrapolation_bc=Flat())

	Ω_new_4 = cubic_spline_interpolation(a_grid, 
		[gauss(x, peak, x_peak, 1e-1) for x in a_grid], extrapolation_bc=Flat())
end 

# ╔═╡ 76986aa4-275c-4fd8-bbec-face050e774f
begin
	plot(a_grid, Ω_new_1(a_grid), xscale=:log10, label = "Std = 1e-4", 
	title = "Gaussian Energy Density", legend=:topleft, xlabel="a", ylabel="Ω_new")

	plot!(a_grid, Ω_new_2(a_grid), xscale=:log10, label = "Std = 1e-3")
	plot!(a_grid, Ω_new_3(a_grid), xscale=:log10, label = "Std = 1e-2")
	plot!(a_grid, Ω_new_4(a_grid), xscale=:log10, label = "Std = 0.1")

	xlims!(1e-10, 1)

end

# ╔═╡ 1706e1fc-7e3f-4a0a-9aef-9c828f18d061
begin
	𝕡4 = CosmoParams(Ω_c = 0.26, Ω_new = a -> Ω_new_1(a))

	𝕡5 = CosmoParams(Ω_c = 0.26, Ω_new = a -> Ω_new_2(a))

	𝕡6 = CosmoParams(Ω_c = 0.26, Ω_new = a -> Ω_new_3(a))

	𝕡7 = CosmoParams(Ω_c = 0.26, Ω_new = a -> Ω_new_4(a))
end

# ╔═╡ dded3ee1-7f99-45b5-8014-21eeff73883c
begin
	bg4, ih4 = FRW_setup(𝕡4);

	bg5, ih5 = FRW_setup(𝕡5);

	bg6, ih6 = FRW_setup(𝕡6);

	bg7, ih7 = FRW_setup(𝕡7);
end

# ╔═╡ ac5bc3ff-acda-4ac4-86f0-0b8a61bef38a
begin
	pL4 = tmap(k -> plin(k, 𝕡4, bg4, ih4), ks1)

	pL5 = tmap(k -> plin(k, 𝕡5, bg5, ih5), ks1)

	pL6 = tmap(k -> plin(k, 𝕡6, bg6, ih6), ks1)

	pL7 = tmap(k -> plin(k, 𝕡7, bg7, ih7), ks1)
end

# ╔═╡ 96e1307c-d048-45ee-bc6f-0360a39368a2
begin
	plot(
    ks1, vcat(pL...),   
    xscale=:log10, yscale=:log10, label="Ω_new = 0",
    xlabel=L"k \ [h/\mathrm{Mpc}]", ylabel=L"P_L(k) \ [\mathrm{Mpc}/h]^3",
	legend=:bottomleft
	)

	plot!(
	    title = "Matter Power Spectrum", ks1, vcat(pL4...),  
	    xscale=:log10, yscale=:log10, label = "std = 1e-4",
		legend=:bottomleft,
	    xlabel=L"k \ [h/\mathrm{Mpc}]", ylabel=L"P_L(k) \ [\mathrm{Mpc}/h]^3"
	)

	plot!(ks1, vcat(pL5...), label = "std = 1e-3", linestyle=:dash, linewidth=2)
	plot!(ks1, vcat(pL6...), label = "std = 1e-2")
	plot!(ks1, vcat(pL7...), label = "std = 0.1", linestyle=:dash)
end

# ╔═╡ eba68a57-4a9f-48ab-b81d-88eaab7ac2cc
begin
	sf4 = source_grid(𝕡4, bg4, ih4, ks, BasicNewtonian())
	Cᵀᵀ4 = cltt(ℓs, 𝕡4, bg4, ih4, sf4)

	sf5 = source_grid(𝕡5, bg5, ih5, ks, BasicNewtonian())
	Cᵀᵀ5 = cltt(ℓs, 𝕡5, bg5, ih5, sf5)

	sf6 = source_grid(𝕡6, bg6, ih6, ks, BasicNewtonian())
	Cᵀᵀ6 = cltt(ℓs, 𝕡6, bg6, ih6, sf6)

	sf7 = source_grid(𝕡7, bg7, ih7, ks, BasicNewtonian())
	Cᵀᵀ7 = cltt(ℓs, 𝕡7, bg7, ih7, sf7)
	
end

# ╔═╡ 29f30d4d-1e7e-4ff4-aa82-0ba18a650b9b
begin
	plot(ℓs, (@. ℓs^2 * Cᵀᵀ), label="Ω_new = 0", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)", legend=:topright)
	
	p_4 = plot!(ℓs, (@. ℓs^2 * Cᵀᵀ4), label="std = 1e-4", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)")

	# Overlapping
	p_5 = plot!(ℓs, (@. ℓs^2 * Cᵀᵀ5), label="std = 1e-3", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)", linestyle=:dash, linewidth=2)

	p_6 = plot!(ℓs, (@. ℓs^2 * Cᵀᵀ6), label="std = 1e-2", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)")

	p_7 = plot!(ℓs, (@. ℓs^2 * Cᵀᵀ7), label="std = 0.1", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)", linestyle=:dash)
end

# ╔═╡ 6ada68b7-da71-4936-b329-9a5826a00bc6
begin
# Hubble parameter
	
plot(bg.x_grid, bg.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="Ω_new = 0", yscale=:log10, legend=:topright)
	
plot!(bg4.x_grid, bg4.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="std = 1e-4", yscale=:log10, legend=:topright)

plot!(bg5.x_grid, bg5.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="std = 1e-3", yscale=:log10)

plot!(bg6.x_grid, bg6.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="std = 1e-2", yscale=:log10)

plot!(bg7.x_grid, bg7.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="std = 0.1", yscale=:log10, linestyle=:dash)		

xlims!(-15, 0)

ylims!(1e-4, 1e2)
end

# ╔═╡ Cell order:
# ╟─a48c8d20-febe-11ee-011b-8b62bce5a872
# ╟─13fe63fb-2335-47b7-b0c1-b8cafcb8f96c
# ╟─8f150cf5-0cda-4cd3-b77b-f89cd892ead3
# ╠═12da0ac6-8a54-40f6-92a5-bde72b5e4370
# ╠═1b4aefd2-77e9-4cee-999e-35b0513ded02
# ╠═fd436b90-809f-4b6d-9b88-f0201bde111d
# ╠═7c11853a-5e2b-432a-a296-6c3bba21101c
# ╠═2ec9bcd4-b0c1-434c-a95d-145ae73c5c66
# ╠═74a2f5b5-46bb-455b-84fb-ad949964d704
# ╠═4086fc84-e831-4ef0-8a0e-6bfe3eee266c
# ╠═ce21d45a-8804-4b2e-92a1-de184f9279ae
# ╠═e32b2579-a2f5-443f-8aa2-0c1a6f215678
# ╠═93a7acd7-37c9-4879-ae3c-7d449e793fe1
# ╠═b0928a06-5c68-4c58-91dd-d3bd1468d060
# ╟─530cf8b5-e7e5-4b05-b69b-c52e144ff92e
# ╠═4c2bd4a5-ba03-4cf4-aebe-d10c53c2a383
# ╠═44b6ddcb-53a5-42ca-a589-254e788507c9
# ╠═bf328343-ae35-42be-b241-02cf144162eb
# ╠═4596a9c2-1046-45bb-8b28-394fe3e5ddc5
# ╠═76986aa4-275c-4fd8-bbec-face050e774f
# ╠═1706e1fc-7e3f-4a0a-9aef-9c828f18d061
# ╠═dded3ee1-7f99-45b5-8014-21eeff73883c
# ╠═ac5bc3ff-acda-4ac4-86f0-0b8a61bef38a
# ╠═96e1307c-d048-45ee-bc6f-0360a39368a2
# ╠═eba68a57-4a9f-48ab-b81d-88eaab7ac2cc
# ╠═29f30d4d-1e7e-4ff4-aa82-0ba18a650b9b
# ╠═6ada68b7-da71-4936-b329-9a5826a00bc6
