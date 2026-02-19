#!/usr/bin/env python3

def make_square_layers_molecular(
    nx=120,
    ny=8,
    nlayers=6,
    a=2.0,            # in-plane spacing (square lattice)
    dz=5.0,           # inter-layer spacing
    origin=(0.0, 0.0, 0.0),
    out_data="6layer.data",
    pad_x=10.0,       # vacuum padding in x
    pad_z=10.0,       # vacuum padding in z
    mass=40.0,         # mass for all atom types
):
    """
    Square lattice layers (4-coordinated in-plane), atom_style molecular.
    Periodic in y; vacuum in x and z.

    Atoms section:
      atom-ID  molecule-ID  atom-type   x   y   z

    Bonds section (within each layer only, distance < bond_cut):
      bond-ID  bond-type  atom1 atom2

    Angles section (within each layer only, 90-degree at center atom):
      angle-ID  angle-type  atom1 atom2 atom3

    Atom type = layer index (1..nlayers) so you can set:
      pair_coeff i i  ...  (strong in-plane)
      pair_coeff i j  ...  (weak interlayer)  for i!=j
    """

    ox, oy, oz = origin

    atoms = []  # (id, mol, type, x, y, z)
    idx = {}    # (k, i, j) -> atom_id

    atom_id = 0
    for k in range(nlayers):
        z = oz + k * dz
        atype = k + 1
        mol = k + 1

        for j in range(ny):
            y = oy + j * a + 1.0
            for i in range(nx):
                x = ox + i * a
                atom_id += 1
                atoms.append((atom_id, mol, atype, x, y, z))
                idx[(k, i, j)] = atom_id

    # Bonds: square net within each layer
    bonds = []
    bond_id = 0

    def add_bond(a1, a2, btype=1):
        nonlocal bond_id
        bond_id += 1
        bonds.append((bond_id, btype, a1, a2))

    # Add top/bot layer, 1-3 connection to prevent folding
    add_bond(idx[(0, 0, 0)], idx[(0, 2, 0)], btype=2)
    add_bond(idx[(0, nx-3, 0)], idx[(0, nx-1, 0)], btype=2)
    add_bond(idx[(nlayers-1, 0, 0)], idx[(nlayers-1, 2, 0)], btype=2)
    add_bond(idx[(nlayers-1, nx-3, 0)], idx[(nlayers-1, nx-1, 0)], btype=2)
    for k in range(nlayers):
        for j in range(ny):
            for i in range(nx):
                a0 = idx[(k, i, j)]

                # Bond to x-neighbor only if i < nx-1 (not periodic in x)
                if i + 1 < nx:
                    add_bond(a0, idx[(k, i + 1, j)], btype=1)
                if k < nlayers - 1:
                    print(f"Adding vertical bond: {a0} to {idx[(k+1, i, j)]}")
                    add_bond(a0, idx[(k+1, i, j)], btype=2)
                    if i > 1: add_bond(a0, idx[(k+1, i-1, j)], btype=2)
                    if i > 2: add_bond(a0, idx[(k+1, i-2, j)], btype=2)
                    if i > 3: add_bond(a0, idx[(k+1, i-3, j)], btype=2)
                    if i > 4: add_bond(a0, idx[(k+1, i-4, j)], btype=2)
                    if i + 1 < nx: add_bond(a0, idx[(k+1, i+1, j)], btype=2)
                    if i + 2 < nx: add_bond(a0, idx[(k+1, i+2, j)], btype=2)
                    if i + 3 < nx: add_bond(a0, idx[(k+1, i+3, j)], btype=2)
                    if i + 4 < nx: add_bond(a0, idx[(k+1, i+4, j)], btype=2)
                # Bond to y-neighbor only if j < ny-1 or j == ny-1 and ny > 1
                # To avoid double-counting with periodic y, only bond j -> j+1 (wrapping)
                if ny == 1:
                    pass  # Single row, no y-bonds
                else:
                    if j < ny - 1:
                        add_bond(a0, idx[(k, i, j + 1)])
                    elif j == ny - 1:
                        # Periodic: only bond if we haven't already bonded j=0 to j=ny-1
                        # Bond j=ny-1 to j=0 only from j=ny-1
                        add_bond(a0, idx[(k, i, 0)])

    # Angles: linear angles in 1D chain
    angles = []
    angle_id = 0
    atype = 1

    def add_angle(a1, a2, a3, atype):
        nonlocal angle_id
        angle_id += 1
        angles.append((angle_id, atype, a1, a2, a3))

    for k in range(nlayers):
        for j in range(ny):
            # For a 1D chain: create angle i-1, i, i+1 at each interior atom
            for i in range(nx):  # atoms with both neighbors
                center = idx[(k, i, j)]
                if 1 < i < nx -1: #i == 0 or i == nx - 1:
                #    if 0 < k < nlayers - 1:
                #        add_angle(idx[(k-1, i, j)], center, idx[(k+1, i, j)], atype=2)
                #        print(f"Added vertical angle: {idx[(k-1, i, j)]} - {center} - {idx[(k+1, i, j)]}")
                #else:
                    left = idx[(k, i - 1, j)]
                    right = idx[(k, i + 1, j)]
                    add_angle(left, center, right, atype=1)
    print(f"Generated {len(atoms)} atoms, {len(bonds)} bonds, {len(angles)} angles")

    n_atoms = len(atoms)
    n_bonds = len(bonds)
    n_angles = len(angles)
    #n_dihedrals = len(dihedrals)

    # Simulation box bounds
    xlo = ox - pad_x
    xhi = ox + (nx - 1) * a + pad_x
    ylo = oy
    yhi = oy + ny * a + 20 # periodic along y
    zlo = oz - pad_z
    zhi = oz + (nlayers - 1) * dz + 30 #pad_z

    with open(out_data, "w") as f:
        f.write("LAMMPS data file: square layers (atom_style molecular)\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write(f"{n_bonds} bonds\n")
        f.write(f"{n_angles} angles\n")
        f.write(f"{nlayers} atom types\n")
        f.write("2 bond types\n")
        f.write("1 angle types\n\n")
        f.write("1 dihedral types\n\n")
        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n\n")

        # Masses
        f.write("Masses\n\n")
        for t in range(1, nlayers + 1):
            f.write(f"{t} {mass:.2f}\n")
        f.write("\n")

        # Atoms
        f.write("Atoms # molecular\n\n")
        for (aid, mol, atype, x, y, z) in atoms:
            f.write(f"{aid} {mol} {atype} {x:.6f} {y:.6f} {z:.6f}\n")
        f.write("\n")

        # Bonds
        f.write("Bonds\n\n")
        for (bid, bt, a1, a2) in bonds:
            f.write(f"{bid} {bt} {a1} {a2}\n")
        f.write("\n")

        # Angles
        f.write("Angles\n\n")
        for (aid, at, a1, a2, a3) in angles:
            f.write(f"{aid} {at} {a1} {a2} {a3}\n")


    print(f"Wrote: {out_data}")
    print(
        "Square lattice layers with bonds (r0=2.0, k=10), angles "
        "(theta=90, k=0.5), and dihedrals (phi=180, k=0.5), "
        "periodic along y; vacuum in x and z."
    )


if __name__ == "__main__":
    nlayers = 20
    make_square_layers_molecular(
        nx=150,
        ny=1,
        nlayers=nlayers,
        a=3.0,
        dz=5.0,
        pad_x=30.0,
        pad_z=200.0,
        out_data=f"{nlayers}layer.data",
    )

