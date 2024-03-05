
# ## Preparation
# Load the packages
using WannierIO
using Wannier
#using WannierPlots

win = read_win("up/cri3_up.win")
model_up = read_w90("up/cri3_up")
model_dn = read_w90("dn/cri3_dn")

omega(model_up)
omega(model_dn)

# position operator
pos_up = Wannier.TBPosition(model_up)
pos_dn = Wannier.TBPosition(model_dn)

# Hamiltonian
#ham_up = Wannier.TBHamiltonian(model_up, gauges=model_up.gauges)
ham_up = Wannier.TBHamiltonian(model_up)
ham_dn = Wannier.TBHamiltonian(model_dn)

# write to tb files
Wannier.write_w90_tb("up/cri3_up", ham_up, pos_up)
Wannier.write_w90_tb("dn/cri3_dn", ham_dn, pos_dn)

