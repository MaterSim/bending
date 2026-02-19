units          real
dimension      3
atom_style     molecular
boundary       f p f

variable nlayers equal 10
read_data      ${nlayers}layer.data

variable e equal 0.5     # interlayer - moderate adhesion for slip with cohesion
variable s equal 4.8     # 4.5Å spacing, atoms near attractive minimum (not repulsive)
variable rcut   equal 20.0

variable k1 equal 50.0    # inplane strech (Morse)
variable k2 equal 0.5     # outplane strech (Morse)
variable k3 equal 3.0     # inplane bending (harmonic)
variable vdown equal 5e-7 # slower for better quasi-static behavior
variable dump_file string "bend-L${nlayers}-s${k1}-b${k3}-l${k2}.dump"

# ----------------------------
# Bonded interactions
# ----------------------------
bond_style     morse
angle_style    harmonic
bond_coeff     1 ${k1} 1.0 3.0
bond_coeff     2 ${k2} 2.0 5.0
angle_coeff    1 ${k3} 180.0
# ----------------------------
# Neighbor / timestep
# ----------------------------
neighbor       0.3 bin
neigh_modify   delay 0 every 1 check yes
timestep       0.5 #fs
variable       dt equal dt

pair_style     lj/cut ${rcut}
pair_modify    shift yes
#special_bonds  lj 0.0 0.0 0.0 # provide short range repulsion for all atoms
pair_coeff     * * ${e} ${s}
fix            ylock all setforce NULL 0.0 NULL

# ----------------------------
# minimization
# ----------------------------
thermo 10
thermo_style custom step pe etotal fmax
dump mindump all custom 1 minimize.dump id type x y z
min_modify dmax 1.0
minimize 1e-6 1e-8 2000 20000
undump mindump

# ----------------------------
# Thermostat / integrator
# ----------------------------
velocity       all create 100 12345 mom yes rot yes dist gaussian
velocity all set NULL 0.0 NULL
compute        t_xz all temp/partial 1 0 1
fix            lange all langevin 100 100 0.5 904297
fix_modify     lange temp t_xz
fix            int   all nve
fix            mom   all momentum 100 linear 1 0 1

compute        zmin all reduce min z
compute        zmax all reduce max z
thermo_style   custom step temp etotal pe ke c_zmin c_zmax 
thermo_modify  temp t_xz
thermo         2000

# ============================================================
# STAGE A: RELAXATION 
# ============================================================

dump eqsnap all custom 300 equi1.dump id type x y z
run 10000
undump eqsnap

# ----------------------------
# THREE-INDENTER BENDING
# ----------------------------
variable gap   equal 0.01

# supports (roller-like): fixed
variable Rind1 equal 8.0
variable Kind1 equal 0.0005

# top indenter: moving - larger radius, softer spring for distributed loading
variable Rind2 equal 50
variable Kind2 equal 0.2

# place based on box size (robust)
variable ys    equal 0.5*(ylo+yhi)
variable xtop  equal 0.5*(xlo+xhi)
variable xs1   equal xlo+0.2*(xhi-xlo)
variable xs2   equal xlo+0.8*(xhi-xlo)

# IMPORTANT: no initial overlap
variable zsup  equal c_zmin-${Rind1}-${gap}
variable ztop0 equal c_zmax+${Rind2}+${gap}

# move down slowly enough to be stable, fast enough to see bending
variable ztop  equal v_ztop0-v_vdown*step

fix sup1 all indent ${Kind1} cylinder y ${xs1} ${zsup} ${Rind1} units box
fix sup2 all indent ${Kind1} cylinder y ${xs2} ${zsup} ${Rind1} units box
dump eqin all custom 300 equi2.dump id type x y z
run 5000
undump eqin

# Freeze current COM drift before loading
velocity all zero linear
velocity all set NULL 0.0 NULL
reset_timestep 0
fix top  all indent ${Kind2} cylinder y ${xtop} v_ztop  ${Rind2} units box
thermo         10000
thermo_style   custom step temp etotal pe ke c_zmin c_zmax v_ztop
dump           d1 all custom 500000 ${dump_file} id type x y z
dump_modify    d1 sort id
run            500000000
