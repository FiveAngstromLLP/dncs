import numpy as np

def parse_pdb_data(data):
    """Parses the provided data and extracts atom coordinates."""
    atoms = []
    for line in data.strip().split("\n"):
        parts = line.split()
        if parts[0] == 'ATOM':
            atom_id = int(parts[1])
            atom_type = parts[2]
            res_name = parts[3]
            res_id = parts[4]
            x = float(parts[5])
            y = float(parts[6])
            z = float(parts[7])
            atoms.append((atom_id, atom_type, res_name, res_id, np.array([x, y, z])))
    return atoms


def rotation_matrix_from_euler_angles(phi, theta, psi):
    """Creates a rotation matrix from Euler angles (Z-X-Z convention)."""
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(theta), -np.sin(theta)],
                   [0, np.sin(theta), np.cos(theta)]])

    Rz1 = np.array([[np.cos(phi), -np.sin(phi), 0],
                    [np.sin(phi), np.cos(phi), 0],
                    [0, 0, 1]])

    Rz2 = np.array([[np.cos(psi), -np.sin(psi), 0],
                   [np.sin(psi), np.cos(psi), 0],
                   [0, 0, 1]])

    return np.dot(Rz2, np.dot(Rx, Rz1))


def apply_transform(v, R, t):
    """Applies rotation matrix and translation vector to the vector."""
    return np.dot(R, v) + t

def create_transformed_pdb(original_pdb, transformed_atoms):
    """Creates a string representation of the transformed PDB structure."""
    transformed_lines = []
    for i, (atom_id, atom_type, res_name, res_id, transformed_coord) in enumerate(transformed_atoms):
        transformed_lines.append(f"ATOM  {atom_id:5d} {atom_type:>2}  {res_name:>3} {res_id:<4}   {transformed_coord[0]:8.3f}{transformed_coord[1]:8.3f}{transformed_coord[2]:8.3f}  1.00  0.00           {atom_type:>2}    ")

    return "\n".join(transformed_lines)


# --- Parameters (These should be the ones you obtained from previous analysis) ---
phi_degrees = 118.13975977238037
theta_degrees = -0.0
psi_degrees = 180.0
translation_vector = np.array([ 7.00000000e+00, -5.42427897e-03, 1.25435651e-02])

# Convert the Euler angles to radians
phi_radians = np.deg2rad(phi_degrees)
theta_radians = np.deg2rad(theta_degrees)
psi_radians = np.deg2rad(psi_degrees)

# --- Original PDB Data (This should be your original PDB data) ---
original_pdb = """
ATOM      1  N   ALA A   2       0.000   0.000   0.000  1.00  0.00           N   14.007  -0.4630 N
ATOM      6  H   ALA A   2       0.547  -0.931  -0.011  1.00  0.00           H    1.008   0.1260 H
ATOM      2  CA  ALA A   2       0.677   1.294   0.013  1.00  0.00           C   12.011   0.0350 CT
ATOM      7  HA  ALA A   2       0.408   1.843  -0.909  1.00  0.00           H    1.008   0.0480 HC
ATOM      5  CB  ALA A   2       0.169   2.072   1.240  1.00  0.00           C   12.011  -0.0980 CT
ATOM      8  HB3 ALA A   2       0.625   3.077   1.307  1.00  0.00           H    1.008   0.0380 HC
ATOM      9  HB1 ALA A   2      -0.927   2.224   1.206  1.00  0.00           H    1.008   0.0380 HC
ATOM     10  HB2 ALA A   2       0.394   1.548   2.189  1.00  0.00           H    1.008   0.0380 HC
ATOM      3  C   ALA A   2       2.177   1.120   0.010  1.00  0.00           C   12.011   0.6160 C
ATOM      4  O   ALA A   2       2.707   0.003  -0.014  1.00  0.00           O   15.999  -0.5040 O
DUMM     11  N   ALA A   3       2.974   2.366   0.023  1.00  0.00           N   14.007  -0.4630 N
"""

# --- Transformation Calculation ---

# Calculate the rotation matrix from the euler angles.
rotation_matrix = rotation_matrix_from_euler_angles(phi_radians, theta_radians, psi_radians)

#parse the original coordinates
atoms = parse_pdb_data(original_pdb)

# Transform each atom
transformed_atoms = []
for atom_id, atom_type, res_name, res_id, original_coord in atoms:
    transformed_coord = apply_transform(original_coord, rotation_matrix, translation_vector)
    transformed_atoms.append((atom_id, atom_type, res_name, res_id, transformed_coord))

# Create the transformed PDB string
transformed_pdb = create_transformed_pdb(original_pdb, transformed_atoms)


# Print the transformed PDB
print(transformed_pdb)
