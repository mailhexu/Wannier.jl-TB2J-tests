# # 10. Co-optimizing spin-up and spin-down WFs

#=
```@meta
CurrentModule = Wannier
```
=#

#=
Usually, for spin-polarized systems, we run two independent Wannierizations
for the spin-up and spin-down channels. However, in some cases, e.g., computing
magnetic exchange constants within an Heisenberg model, we need to make sure
the spin-up and spin-down WFs have as similar as possible real-space shapes,
and are centered on each atoms.

In this tutorial, we will show how to co-optimize the spin-up and spin-down
WFs, with two constraints:
1. WF center constraint to construct atom-centered WFs
2. spin-up and down overlap constraint to make sure each pair of
    spin-up and spin-down WFs are as similar as possible

We will Wannierize a 2D ``CrI_3`` system, using pseudo-atomic-orbital projections
as the starting guess (computed by QE).

## Outline

1. plot QE band structure as a reference
2. run two independent Wannierizations of spin-up and spin-down channels
3. construct a [`MagModel`](@ref) that merges the two spin channels
4. disentangle, with overlap constraint
5. disentangle, with both WF center and overlap constraints

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`tutorial.ipynb`](./tutorial.ipynb)
    - Julia script: [`tutorial.jl`](./tutorial.jl)
=#

# ## Preparation
# Load the packages
using WannierIO
using Wannier
#using WannierPlots

#=
## Plot QE band structure

It's always good to check band-structure-interpolation qualities
of WFs against DFT bands.
To do this, let's first load QE band structure
=#
qe = WannierIO.read_qe_xml("qe_bands.xml")

#=
Note that I used Wannier90 generated kpath for QE bands calculation, to be
comparable with Wannier-interpolated bands.
To generate a kpath equivalent to Wannier90, here we use the [`get_kpath`](@ref)
function, with `unit_cell` and `kpoint_path` parsed from `win` file
by [`read_win`](@ref) function.
=#
win = read_win("up/cri3_up.win")
#kpath = Wannier.get_kpath(win.unit_cell, win.kpoint_path)

#=
then we construct an [`KPathInterpolant`](@ref) object which stores the exact
kpoint coordinates to be interpolated, using 100 points in the 1st kpath segment
(equivalent to Wannier90 input `bands_num_points`).
=#
#kpi = Wannier.generate_w90_kpoint_path(kpath, 100)

# Now we can plot the QE bands for two spin channels
#P = plot_band_diff(kpi, qe.E_up, qe.E_dn; fermi_energy=qe.fermi_energy)
#Main.HTMLPlot(P, 500)  # hide

#=
## Model construction

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct two [`Model`](@ref)s
for spin-up and spin-down channels.
Note the frozen windows for the two channels are set independently
according to the `dis_froz_max` inside the two `win` files, respectively;
in our case, they are both `dis_froz_max = -2 eV` since there are gaps
in both spin channels.
=#
model_up = read_w90("up/cri3_up")
model_dn = read_w90("dn/cri3_dn")

#=
## Projection-only WFs

The projection-only WFs, i.e., without any disentanglement or maximal
localization, are usually centered on each atom. This can be checked
by computing WF centers and spreads
=#
omega(model_up)

# and for spin-down channel
omega(model_dn)

pos_up = Wannier.TBPosition(model_up)
pos_dn = Wannier.TBPosition(model_dn)


#ham_up = Wannier.TBHamiltonian(model_up, gauges=model_up.gauges)
ham_up = Wannier.TBHamiltonian(model_up)
ham_dn = Wannier.TBHamiltonian(model_dn)


Wannier.write_w90_tb("up/cri3_up", ham_up, pos_up)
Wannier.write_w90_tb("dn/cri3_dn", ham_dn, pos_dn)


#=
However, their band interpolations are often not good, so they shouldn't
be used for magnetic exchange constants calculations.

In our case, the projection-only WFs are `Cr:4s,3d` and `I:5s,5p` orbitals
centered on atoms. We will use Wannier interpolation for band structures,
and compare them with QE bands.

We first construct two [`InterpModel`](@ref)s for spin-up and spin-down channels,
using the `kpoint_path` from `win` file (otherwise, by default the `InterpModel`
will use `Brillouin.jl` to auto generate a kpath, which might be different
from user's input)
=#
#interpModel_up = Wannier.InterpModel(model_up; kpath=kpath)
#interpModel_dn = Wannier.InterpModel(model_dn; kpath=kpath)

# then interpolate eigenvalues
#E_up_projonly = Wannier.interpolate(interpModel_up, kpi)
#E_dn_projonly = Wannier.interpolate(interpModel_dn, kpi)

# and plot the spin-up bands compared with QE
#P = plot_band_diff(kpi, qe.E_up, E_up_projonly; fermi_energy=qe.fermi_energy)
#Main.HTMLPlot(P, 500)  # hide

# and the spin-down bands
#P = plot_band_diff(kpi, qe.E_dn, E_dn_projonly; fermi_energy=qe.fermi_energy)
#Main.HTMLPlot(P, 500)  # hide

#=
As can be seen from the above two figures, the projection-only WFs do not
reproduce DFT bands, i.e., they do not correctly describe the electronic
structure of the system, thus should not be used for physical property
calculations.

As a side node, we can also plot the Wannier-interpolated spin-up and down bands
in one figure,
=#
#P = plot_band_diff(kpi, E_up_projonly, E_dn_projonly; fermi_energy=qe.fermi_energy)
#Main.HTMLPlot(P, 500)  # hide

#=
Save the MLWF gauge into `chk` files, to allow other codes, e.g., Wannier90,
to restart from the generated gauge.
Note we need to explicitly provide the `exclude_bands`, so as to be consistent
with `win` file.

!!! note

    The written `chk` are Fortran formatted files, which you can convert to
    binary file using Wannier90 `chk2chk.x` executable.
    There is an additional keyword argument `binary` for [`write_chk`](@ref),
    which is `false` by default, but can be set to `true` to write binary `chk`;
    however, Fortran binary format is compiler-dependent, so it's not gaurantted
    to be able to restart from the written binary `chk` file for other codes.
=#
exclude_bands = collect(1:8)
Wannier.write_chk("up/wjl_up_projonly.chk", model_up; exclude_bands)
Wannier.write_chk("dn/wjl_dn_projonly.chk", model_dn; exclude_bands)


