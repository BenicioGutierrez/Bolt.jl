using DelimitedFiles
using PyPlot
cd("test")

recfastdata = readdlm("data/test_recfast_1.dat", ',', Float64, '\n', header=true)[1]
z⃗, Xe = recfastdata[:,1], recfastdata[:,2]

# clf()
# plot(z, Xe, "-", label=raw"$X_e$")
# legend()
# gcf()

##
using Bolt

𝕡 = CosmoParams(Σm_ν = 0.0, N_ν = 3.0)
bg = Background(𝕡)
𝕣 = Bolt.RECFAST(bg=bg)

##
Nz = 1000
xe_bespoke, Tmat = Bolt.recfast_xe(𝕣; Nz=Nz, zinitial=10000., zfinal=0.);

# z⃗ = 10000.0-10.0:-10.0:0.0
dz = (0. - 10000.)/float(Nz)
z⃗ = (10000. + dz):(dz):0.0

clf()
plot(z⃗, Tmat, "-", label=raw"$T_{\mathrm{mat}}$")
plot(z⃗, 𝕣.Tnow .* (1 .+ z⃗), "--", label=raw"$T_{\mathrm{rad}}$")

yscale("log")
xscale("log")
legend()
ylim(1, 2e4)
xlim(10, 10000)
xlabel("redshift")
ylabel("temperature [K]")
gcf()

##
clf()
plot(z⃗, Tmat ./ (𝕣.Tnow .* (1 .+ z⃗)), "-")
xscale("log")
legend()
xlim(10, 10000)
xlabel("redshift")
ylabel(raw"$T_{\mathrm{mat}} \, / \, T_{\mathrm{rad}}$")
gcf()

##

clf()
plot(z⃗, Xe ./ xe_bespoke , "-", label=raw"RECFAST / recfast.jl")
ylim(1 - 0.01, 1 + 0.01)

# plot(z, Xe , "-")
# plot(z, xe_bespoke, "--")
xlabel(raw"redshift")
legend()
gcf()

##
using UnitfulAstro, NaturallyUnitful
x₁ = let z = 100.0
    Hz = 𝕣.HO * sqrt((1+z)^4/(1+𝕣.z_eq)*𝕣.OmegaT + 𝕣.OmegaT*(1+z)^3 + 𝕣.OmegaK*(1+z)^2 + 𝕣.OmegaL)
    (𝕣.HO^2 /2/Hz)*(4*(1+z)^3/(1+𝕣.z_eq)*𝕣.OmegaT + 3*𝕣.OmegaT*(1+z)^2 + 2*𝕣.OmegaK*(1+z))
end

##
x₂ = let z = 100.0
    H0_natural_unit_conversion = ustrip(u"s", unnatural(u"s", 1u"eV^-1"))

    a = 1 / (1+z)  # scale factor
    x_a = a2x(a)
	Hz = 𝕣.bg.ℋ(x_a) / a / 𝕣.H0_natural_unit_conversion
	dHdz = (-𝕣.bg.ℋ′(x_a) + 𝕣.bg.ℋ(x_a)) / 𝕣.H0_natural_unit_conversion
end

##
let z = 0.0
    a = 1.0
    x_a = a2x(a)
    H0_natural_unit_conversion = ustrip(u"s", unnatural(u"s", 1u"eV^-1"))
	Hz = 𝕣.bg.H₀ / 𝕣.H0_natural_unit_conversion

end
##
𝕣.HO
