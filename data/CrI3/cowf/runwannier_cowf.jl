# ## Preparation
# Load the packages
using WannierIO
using Wannier
#using WannierPlots
#

win = read_win("up/cri3_up.win")
model_up = read_w90("up/cri3_up")
model_dn = read_w90("dn/cri3_dn")

#The spin-up and down overlap matrices is written in the same format as `amn`
#
Mud = read_amn("cri3_updn.mud");

# then assemble into a [`MagModel`](@ref)
model = Wannier.MagModel(model_up, model_dn, Mud)

#=
Now let's disentangle with spin overlap constraint.
Here `λs` is the Lagrange multiplier for the constraint.
=#
λs = 100.0
U_up, U_dn = disentangle(model, λs;max_iter=4000);
#=
The resulting spin-up and spin-down WFs have very similar centers and spreads,
however, their centers drift from the original positions which were centered
on atoms.
=#
omega(model, U_up, U_dn, λs)

# position operator
pos_up = Wannier.TBPosition(model_up, gauges=U_up)
pos_dn = Wannier.TBPosition(model_dn, gauges=U_dn)

# Hamiltonian
ham_up = Wannier.TBHamiltonian(model_up, U_up)
ham_dn = Wannier.TBHamiltonian(model_dn, U_dn)

# write to tb files
Wannier.write_w90_tb("up/cri3_up", ham_up, pos_up)
Wannier.write_w90_tb("dn/cri3_dn", ham_dn, pos_dn)

