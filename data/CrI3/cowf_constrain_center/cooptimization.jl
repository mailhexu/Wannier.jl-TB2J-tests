using WannierIO
using Wannier
using Wannier.Datasets
#using WannierPlots

model_up = read_w90(dataset"CrI3/CrI3_up")
model_dn = read_w90(dataset"CrI3/CrI3_dn")

omega(model_up)
omega(model_dn)

Mud = read_amn(dataset"CrI3/CrI3_updn.mud");

model = Wannier.MagModel(model_up, model_dn, Mud)

r₀ = [zeros(3) for _ in 1:n_wannier(model.up)]
# the first 6 WFs are the `4s,3d` orbitals of the 1st `Cr` atom
r₀[1:6] .= Ref(model_up.atom_positions[1])
#=
!!! note
    Here I am using julia's `Ref` to treat the `atom_positions[1]` as a scaler
    when broadcasting to the 1th to 6th elements of vector `r₀`. The result is
    the same as writing a `for` loop over the first six elements.
=#
# the next 6 WFs are the `4s,3d` orbitals of the 2nd `Cr` atom
r₀[7:12] .= Ref(model_up.atom_positions[2])
# the next 4 WFs are the `5s,5p` orbitals of `I` atom
r₀[13:16] .= Ref(model_up.atom_positions[3])
# and similarly for the remaining 5 `I` atoms
r₀[17:20] .= Ref(model_up.atom_positions[4])
r₀[21:24] .= Ref(model_up.atom_positions[5])
r₀[25:28] .= Ref(model_up.atom_positions[6])
r₀[29:32] .= Ref(model_up.atom_positions[7])
r₀[33:36] .= Ref(model_up.atom_positions[8])

# convert to Cartesian coordinates
r₀ = map(r₀) do v
    model.up.lattice * v
end

λc = 10.0
λs = 10.0
U_up, U_dn = Wannier.disentangle_center(model, r₀, λc, λs);

# the final centers and spreads are
omega(model, U_up, U_dn, r₀, λc, λs)

# as a comparison, the spreads of independent Wannierizations are
omega(model, U_up_mlwf, U_dn_mlwf, r₀, λc, λs)


