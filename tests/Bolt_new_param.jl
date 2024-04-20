### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 28f3f7b0-eb37-11ee-3618-791d0ed69bcb
import Pkg

# ╔═╡ a4bf84ce-85ad-4d87-88bc-23c4b001f364
begin
	Pkg.add(url="https://github.com/BenicioGutierrez/Bolt.jl")
	Pkg.add("ForwardDiff")
	Pkg.add("Plots")
	Pkg.add("DataInterpolations")
end

# ╔═╡ 6df97ffd-fff0-4732-adaa-fb301e711383
using Bolt, Plots, ForwardDiff, LaTeXStrings, ThreadTools, Base.Threads, DataInterpolations, Interpolations

# ╔═╡ e04e40bb-cb03-473c-9725-0e602fe7e5fd
md"""
# Adding a new parameter
"""

# ╔═╡ 41c2c4de-b759-4d37-8a90-3774f3ca5ec5
md"""
The goal is to add a new energy density that is time dependent. This approach inserts the energy density into CosmoParams
"""

# ╔═╡ d678de54-ca81-4251-ae2d-bf0ef3e6588f
md"""
# Test 1: Zero
"""

# ╔═╡ 2bfe72f2-d2df-4dd2-be65-3bddcabfedf4
md"""
Original case, where there is no contribution from the new energy density
"""

# ╔═╡ 37498b7d-1aed-49fa-8181-a0d2e581b978
# Assign cosmological parameters
begin
	𝕡 = CosmoParams(Ω_c = 0.26) # set kwargs like so to change the default values
end

# ╔═╡ 3786ea80-640b-4c66-9bbd-1f10896dcfb0
function FRW_setup(𝕡)
    # Compute expansion history quantities
    bg = Background(𝕡)
    # Compute ionization history (via RECFAST)
#    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
	𝕣 = Bolt.RECFAST(bg; Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b,OmegaG=𝕡.Ω_r)  # Check RECFAST
    ih = IonizationHistory(𝕣, 𝕡, bg)
    return bg, ih
end

# ╔═╡ a183ddc1-0ba9-433f-bf75-79122104f0d6
begin
	a_min = 10.0^(-11)
	a_max = 1.00
	n_a = 100
	a = log_a(a_min, a_max, n_a) # a grid
end

# ╔═╡ 7bb554d4-5a01-40f1-b78b-c8521483444d
bg, ih = FRW_setup(𝕡);

# ╔═╡ 2f52fca6-f847-4367-8ca0-49f0eb292d25
# Matter power spectrum

begin
	kmin1 = 1bg.H₀
	kmax1 = 5000bg.H₀
	nk1 = 32
	ks1 = log10_k(kmin1, kmax1, nk1) # k grid
end

# ╔═╡ 493235a6-d653-4871-a213-391a45493351
md"""
## Power Matter Spectrum
"""

# ╔═╡ 5d4b54cd-8f65-402c-b99b-a45611379615
pL = tmap(k -> plin(k, 𝕡, bg, ih), ks1)

# ╔═╡ 905e5daf-878b-427a-9b5f-91c0b5912675
md"""
### Gradient in Relation to Ωc
"""

# ╔═╡ 7c5926d0-c989-42a4-9965-3089254a1254
# Gradient wrt Ω_c
# Define a function that changes 𝕡 - need to recompute background components, as they depend on Ω_c
function pL_Ω_c(Ω_c::T) where T
    𝕡 = CosmoParams{T}(Ω_c=Ω_c)
    bg, ih = FRW_setup(𝕡)
    return tmap(k -> plin(k, 𝕡, bg, ih)[1], ks1)
end

# ╔═╡ 3e11edaf-f2bc-48b2-b1bb-52519291cea9
∂pL_∂Ω_c = ForwardDiff.derivative(pL_Ω_c, 0.26);

# ╔═╡ ba91b8aa-2e60-4be3-a1a6-ede4e3372893
plot(
    ks1, abs.(∂pL_∂Ω_c), 
    xscale=:log10, yscale=:log10, label=false,
    xlabel=L"k \ [h/\mathrm{Mpc}]", ylabel=L"\vert \partial_{\Omega_c} P_L(k) \vert"
)

# ╔═╡ f232fb88-53d6-4a96-b7de-ffa69f8f6bf6
md"""
##  Angular power spectrum of the temperature fluctuations
"""

# ╔═╡ 6e4ba72c-2afb-4db6-be85-61b9c0445132
# CMB Cᵀᵀ(ℓ)
begin 
	ℓmin, ℓmax, nℓ = 2, 20, 1200
	ℓs = ℓmin:ℓmax:nℓ
	kmin, kmax, nk = 0.1bg.H₀, 1000bg.H₀, 100
	ks = quadratic_k(kmin, kmax, nk)
	sf = source_grid(𝕡, bg, ih, ks, BasicNewtonian()) # set up LOS source function interpolator
end

# ╔═╡ 1435e76f-40f5-49d0-9259-5b0300d61f6d
Cᵀᵀ = cltt(ℓs, 𝕡, bg, ih, sf)

# ╔═╡ 1160e265-0f75-4dcc-8958-7438c2a53131
md"""
### Gradient with relation to Ωb
"""

# ╔═╡ 70254a89-e340-41fd-8c0b-afd99ffd9e9d
# gradient wrt Ω_b
function Cᵀᵀ_Ω_b(Ω_b::T) where T # type-stable wrapper
    𝕡 = CosmoParams{T}(Ω_b=Ω_b)
    bg, ih = FRW_setup(𝕡)
    sf = source_grid(𝕡, bg, ih, ks, BasicNewtonian())
    return cltt(ℓs, 𝕡, bg, ih, sf)
end

# ╔═╡ e1688dc9-ed22-47ae-b9c1-310444e5281a
∂Cᵀᵀ_∂Ω_b = ForwardDiff.derivative(Cᵀᵀ_Ω_b,0.045);

# ╔═╡ 1c53378b-ad8d-4a9a-bd2a-90af603c233e
plot(
    ℓs, (@. ℓs^2 * ∂Cᵀᵀ_∂Ω_b), 
    label=false, xlabel=L"\ell", ylabel=L"\ell^2 \partial_{\Omega_b} C^{TT}(\ell)"
)

# ╔═╡ 7a82f4f2-5623-445b-97d8-378b6e8ba7df
md"""
### Conformal Time
"""

# ╔═╡ 99899ddc-9184-46aa-a5b6-22f74143ea4f
# Conformal time
plot(bg.x_grid, bg.η, xlabel=L"\log(a)", ylabel=L"\eta", label=false, yscale=:log10)

# ╔═╡ d8999b26-c710-4e37-9dd2-9837fe697f0c
md"""
### Hubble Parameter
"""

# ╔═╡ d0387f8e-6b91-429b-abc5-e361f676da44
md"""
### Free electron fraction
"""

# ╔═╡ 77e40eb9-1830-428b-b0b5-22f1f55b75b6
# Free electron fraction
plot(bg.x_grid, ih.Xₑ, xlabel=L"\log(a)", ylabel=L"X_{e}", label=false)

# ╔═╡ 62903351-98dd-4354-8379-51ccb8542470
md"""
### Visibiliy Function
"""

# ╔═╡ 507b0cd2-5068-46b1-b232-1fbbb1e206fc
begin
	# Visibility function
	plot(bg.x_grid, ih.g̃, xlabel=L"\log(a)", ylabel=L"g", label=false)
	xlims!(-10, 0) 
end

# ╔═╡ 2931a76b-79b5-4c43-a5b4-bddf4f2d9499
md"""
# Test 2: Step Function
"""

# ╔═╡ 4c5338f5-f15d-4f3e-9433-abb0e098a1e2
md"""
Ω_new is zero with redshift z < 1200, 0.3 if greater
"""

# ╔═╡ 1e253c9e-2ccc-4f1d-afd0-4e06eca67542
# Assign new cosmological parameters
begin
	𝕡1 = CosmoParams(Ω_c = 0.26, Ω_new = a -> (a > 0.000833) ? 0.3 : 0.0)  
	# set kwargs like so to change the default values

	# Heaviside step function: Ω_new = a -> (a > 0.000833) ? 0.3 : 0.0 
	# Step function: Ω_new = a -> (a > 0.000833 && a < 0.0099) ? 0.3 : 0.0
end

# ╔═╡ c3f8368b-013a-42ab-801f-8a4aca035962
begin 
	# Generate values of a
	a_values = exp10.(range(log10(1e-5), log10(1), length=100))

	# Calculate corresponding values of Ω_new
	Ω_values = 𝕡1.Ω_new.(a_values)

	# Compute ln(a)
	ln_a_values = log.(a_values)

	# Plot Ω_new vs ln(a)
	plot(ln_a_values, 𝕡1.Ω_new.(a_values), xlabel="ln(a)", ylabel="Ω_new", label="Ω_new vs ln(a)", legend=:topleft)
end

# ╔═╡ 12ab3afd-df46-48ce-b922-774f6e5defe8
bg1, ih1 = FRW_setup(𝕡1);

# ╔═╡ ea7a4b40-4a55-4864-bf6a-92f0aea6f008
md"""
## Power Matter Spectrum
"""

# ╔═╡ 8f84cc40-649a-461c-9d56-581738f0a216
pL1 = tmap(k -> plin(k, 𝕡1, bg1, ih1), ks1)

# ╔═╡ e107e2be-f5c2-43ac-9b7c-ace9e29abb36
# CMB Cᵀᵀ(ℓ)
begin 
	kmin3, kmax3, nk3 = 0.1bg1.H₀, 1000bg1.H₀, 100
	ks3 = quadratic_k(kmin3, kmax3, nk3)
	sf1 = source_grid(𝕡1, bg1, ih1, ks3, BasicNewtonian()) # set up LOS source function interpolator
end

# ╔═╡ e8af4c77-36f9-4226-98f4-bae2271227bf
md"""
##  Angular power spectrum of the temperature fluctuations
"""

# ╔═╡ 4f577a38-4f2a-4b8c-939d-fcc8753f0c03
Cᵀᵀ1 = cltt(ℓs, 𝕡1, bg1, ih1, sf1)

# ╔═╡ cf1e0ed5-3de9-43b6-bc11-b4f1c862cfc4
md"""
### Conformal Time
"""

# ╔═╡ c99886fb-61de-440e-b87b-f6ee5f5e5442
# Conformal time
plot(bg1.x_grid, bg1.η, xlabel=L"\log(a)", ylabel=L"\eta", label=false, yscale=:log10)

# ╔═╡ 22705027-626f-405f-990c-16b1eac09372
md"""
### Hubble Parameter
"""

# ╔═╡ 8b6d081b-0d96-49d2-9145-29cdda98f522
md"""
### Free electron fraction
"""

# ╔═╡ 99ebc0b7-c4e7-426c-9048-17b39d8abb17
# Free electron fraction
plot(bg1.x_grid, ih1.Xₑ, xlabel=L"\log(a)", ylabel=L"X_{e}", label=false)

# ╔═╡ c7bbd62b-87f0-4fd4-8441-05dfba63efa2
md"""
### Visibiliy Function
"""

# ╔═╡ bbd4bfbf-9d33-42ff-9e11-1f8330f256f6
begin
	# Visibility function
	plot(bg1.x_grid, ih1.g̃, xlabel=L"\log(a)", ylabel=L"g", label=false)
	xlims!(-10, 0) 
end

# ╔═╡ 99e3ba4e-0d70-492a-b46f-383cec9ec578
md"""
# Test 3: Huge Step Function
"""

# ╔═╡ 7f7af48c-0b06-492f-82e0-cf706c768428
# Assign new cosmological parameters
begin
	𝕡2 = CosmoParams(Ω_c = 0.26 , Ω_new = a -> (a >= 0.000833) ? 1 : 0.0)  
	# set kwargs like so to change the default values

	# Heaviside step function: Ω_new = a -> (a > 0.000833) ? 1 : 0.0 
	# Step function: Ω_new = a -> (a > 0.000833 && a < 0.0099) ? 1 : 0.0
end

# ╔═╡ 028c677a-dbe5-4422-9507-1c62c8a15c2b
bg2, ih2 = FRW_setup(𝕡2);

# ╔═╡ 21df2a22-c182-463f-9b97-37b19f8c3814
pL2 = tmap(k -> plin(k, 𝕡2, bg2, ih2), ks1)

# ╔═╡ 733cb93e-6956-4506-97e5-1d77f9b1dce4
# CMB Cᵀᵀ(ℓ)
begin 
	sf2 = source_grid(𝕡2, bg2, ih2, ks, BasicNewtonian()) # set up LOS source function interpolator
end

# ╔═╡ af570fa3-7160-4643-9a5f-f7a78ebc6f97
Cᵀᵀ2 = cltt(ℓs, 𝕡2, bg2, ih2, sf2)

# ╔═╡ fbd7d59d-b8df-4351-bdb2-0fda7b378096
md"""
# Plots - Step Function
"""

# ╔═╡ 83f04d4a-5a40-42e0-a3d8-c13722c092e6
md"""
## Power Matter Spectrum
"""

# ╔═╡ ea085472-205f-4abd-91c6-dac749b1f976
begin

	
# Plotting the first data
plot(
    ks1, vcat(pL...),   
    xscale=:log10, yscale=:log10, label="Ω_new = 0",
    xlabel=L"k \ [h/\mathrm{Mpc}]", ylabel=L"P_L(k) \ [\mathrm{Mpc}/h]^3",
	legend=:bottomleft
)

# Overlapping the second data
plot!(
    ks1, vcat(pL1...),  
    xscale=:log10, yscale=:log10, label="Ω_new = 0.3, ln(a) = -7.09"
)

plot!(
    ks1, vcat(pL2...),  
    xscale=:log10, yscale=:log10, label="Ω_new = 1, ln(a) = -7.09"
)
		
end

# ╔═╡ 407808e8-18ea-4916-b9cd-17a183ed0f45
md"""
## Angular power spectrum of the temperature fluctuations
"""

# ╔═╡ c6e4fcf0-c9fa-4a28-ab19-8616cfdf23bc
begin

	# Plotting the first data
	p2 = plot(ℓs, (@. ℓs^2 * Cᵀᵀ), label="Ω_new = 0", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)", legend=:topright)

	# Overlapping the second data
	p4 = plot!(ℓs, (@. ℓs^2 * Cᵀᵀ1), label="Ω_new = 0.3, ln(a) = -7.09", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)")

	p5 = plot!(ℓs, (@. ℓs^2 * Cᵀᵀ2), label="Ω_new = 1, ln(a) = -7.09", xlabel=L"\ell", ylabel=L"\ell^2 C^{TT}(\ell)")
		
end

# ╔═╡ 9fa20b01-7469-46eb-a755-72f59be4db8d
md"""
## Hubble Parameter
"""

# ╔═╡ 20a648ba-d366-48a2-96b8-e6d12e2f2bbc
begin

	
# Plotting the first data
# Hubble parameter
plot(bg.x_grid, bg.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="Ω_new = 0", yscale=:log10, legend=:topright)

# Overlapping the second data
# Hubble parameter
plot!(bg1.x_grid, bg1.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="Ω_new = 0.3, ln(a) = -7.09", yscale=:log10)

plot!(bg2.x_grid, bg2.ℋ, xlabel=L"\log(a)", ylabel=L"\mathcal{H}", label="Ω_new = 1, ln(a) = -7.09", yscale=:log10)
		
end

# ╔═╡ Cell order:
# ╟─e04e40bb-cb03-473c-9725-0e602fe7e5fd
# ╟─41c2c4de-b759-4d37-8a90-3774f3ca5ec5
# ╠═28f3f7b0-eb37-11ee-3618-791d0ed69bcb
# ╠═a4bf84ce-85ad-4d87-88bc-23c4b001f364
# ╠═6df97ffd-fff0-4732-adaa-fb301e711383
# ╟─d678de54-ca81-4251-ae2d-bf0ef3e6588f
# ╟─2bfe72f2-d2df-4dd2-be65-3bddcabfedf4
# ╠═37498b7d-1aed-49fa-8181-a0d2e581b978
# ╠═3786ea80-640b-4c66-9bbd-1f10896dcfb0
# ╠═a183ddc1-0ba9-433f-bf75-79122104f0d6
# ╠═7bb554d4-5a01-40f1-b78b-c8521483444d
# ╠═2f52fca6-f847-4367-8ca0-49f0eb292d25
# ╟─493235a6-d653-4871-a213-391a45493351
# ╠═5d4b54cd-8f65-402c-b99b-a45611379615
# ╟─905e5daf-878b-427a-9b5f-91c0b5912675
# ╠═7c5926d0-c989-42a4-9965-3089254a1254
# ╠═3e11edaf-f2bc-48b2-b1bb-52519291cea9
# ╠═ba91b8aa-2e60-4be3-a1a6-ede4e3372893
# ╟─f232fb88-53d6-4a96-b7de-ffa69f8f6bf6
# ╠═6e4ba72c-2afb-4db6-be85-61b9c0445132
# ╠═1435e76f-40f5-49d0-9259-5b0300d61f6d
# ╟─1160e265-0f75-4dcc-8958-7438c2a53131
# ╠═70254a89-e340-41fd-8c0b-afd99ffd9e9d
# ╠═e1688dc9-ed22-47ae-b9c1-310444e5281a
# ╠═1c53378b-ad8d-4a9a-bd2a-90af603c233e
# ╟─7a82f4f2-5623-445b-97d8-378b6e8ba7df
# ╠═99899ddc-9184-46aa-a5b6-22f74143ea4f
# ╟─d8999b26-c710-4e37-9dd2-9837fe697f0c
# ╟─d0387f8e-6b91-429b-abc5-e361f676da44
# ╠═77e40eb9-1830-428b-b0b5-22f1f55b75b6
# ╟─62903351-98dd-4354-8379-51ccb8542470
# ╠═507b0cd2-5068-46b1-b232-1fbbb1e206fc
# ╟─2931a76b-79b5-4c43-a5b4-bddf4f2d9499
# ╟─4c5338f5-f15d-4f3e-9433-abb0e098a1e2
# ╠═1e253c9e-2ccc-4f1d-afd0-4e06eca67542
# ╠═c3f8368b-013a-42ab-801f-8a4aca035962
# ╠═12ab3afd-df46-48ce-b922-774f6e5defe8
# ╟─ea7a4b40-4a55-4864-bf6a-92f0aea6f008
# ╠═8f84cc40-649a-461c-9d56-581738f0a216
# ╠═e107e2be-f5c2-43ac-9b7c-ace9e29abb36
# ╟─e8af4c77-36f9-4226-98f4-bae2271227bf
# ╠═4f577a38-4f2a-4b8c-939d-fcc8753f0c03
# ╟─cf1e0ed5-3de9-43b6-bc11-b4f1c862cfc4
# ╠═c99886fb-61de-440e-b87b-f6ee5f5e5442
# ╟─22705027-626f-405f-990c-16b1eac09372
# ╟─8b6d081b-0d96-49d2-9145-29cdda98f522
# ╠═99ebc0b7-c4e7-426c-9048-17b39d8abb17
# ╟─c7bbd62b-87f0-4fd4-8441-05dfba63efa2
# ╠═bbd4bfbf-9d33-42ff-9e11-1f8330f256f6
# ╟─99e3ba4e-0d70-492a-b46f-383cec9ec578
# ╠═7f7af48c-0b06-492f-82e0-cf706c768428
# ╠═028c677a-dbe5-4422-9507-1c62c8a15c2b
# ╠═21df2a22-c182-463f-9b97-37b19f8c3814
# ╠═733cb93e-6956-4506-97e5-1d77f9b1dce4
# ╠═af570fa3-7160-4643-9a5f-f7a78ebc6f97
# ╟─fbd7d59d-b8df-4351-bdb2-0fda7b378096
# ╟─83f04d4a-5a40-42e0-a3d8-c13722c092e6
# ╠═ea085472-205f-4abd-91c6-dac749b1f976
# ╟─407808e8-18ea-4916-b9cd-17a183ed0f45
# ╠═c6e4fcf0-c9fa-4a28-ab19-8616cfdf23bc
# ╟─9fa20b01-7469-46eb-a755-72f59be4db8d
# ╠═20a648ba-d366-48a2-96b8-e6d12e2f2bbc
