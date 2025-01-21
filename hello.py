import numpy as np

def parse_pdb_data(data):
    """Parses the provided data and extracts atom coordinates."""
    atoms1 = {}
    atoms2 = {}
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
            if int(parts[1]) < 12:
                atoms1[atom_type] = (atom_type, res_name, res_id, np.array([x, y, z]))
            else:
                atoms2[atom_type] = (atom_type, res_name, res_id, np.array([x, y, z]))
    return atoms1, atoms2


def best_fit_transform(A, B):
    """Calculates the rotation and translation that best fits the point sets A and B."""

    # convert to numpy arrays
    A = np.array([coord for _,_,_,coord in A.values()]).T
    B = np.array([coord for _,_,_,coord in B.values()]).T

    # find the centroids
    centroid_A = np.mean(A, axis=1, keepdims = True)
    centroid_B = np.mean(B, axis=1, keepdims = True)

    # subtract the centroids
    AA = A - centroid_A
    BB = B - centroid_B

    # calculate the rotation matrix
    H = np.dot(AA, BB.T)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # handle the special reflection case
    if np.linalg.det(R) < 0:
       Vt[2,:] *= -1
       R = np.dot(Vt.T, U.T)

    # calculate the translation vector
    t = centroid_B - np.dot(R, centroid_A)

    return R, t

def euler_angles_from_rotation_matrix(R):
    """Extracts Euler angles (ZXZ convention) from a rotation matrix."""
    if R[2, 0] < 1 and R[2, 0] > -1:
        theta = -np.arcsin(R[2, 0])
        phi = np.arctan2(R[2, 1] / np.cos(theta), R[2, 2] / np.cos(theta))
        psi = np.arctan2(R[1, 0] / np.cos(theta), R[0, 0] / np.cos(theta))
    elif R[2, 0] == 1:
        theta = np.pi / 2
        psi = 0
        phi = np.arctan2(R[0, 1], R[0, 2])
    elif R[2, 0] == -1:
        theta = -np.pi / 2
        psi = 0
        phi = -np.arctan2(R[0, 1], R[0, 2])
    return (phi, theta, psi)


def apply_transform(v, R, t):
    """Applies the rotation matrix and the translation vector to the vector."""
    return np.dot(R, v) + t


# Sample Data
data = """
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
ATOM     12  N   ALA A   5       0.000   0.000   0.000  1.00  0.00           N   14.007  -0.4630 N
ATOM     17  H   ALA A   5      -0.562   0.921   0.043  1.00  0.00           H    1.008   0.1260 H
ATOM     13  CA  ALA A   5       1.460  -0.014  -0.029  1.00  0.00           C   12.011   0.0350 CT
ATOM     18  HA  ALA A   5       1.826  -0.536   0.875  1.00  0.00           H    1.008   0.0480 HC
ATOM     16  CB  ALA A   5       1.893  -0.796  -1.282  1.00  0.00           C   12.011  -0.0980 CT
ATOM     19  HB3 ALA A   5       2.993  -0.866  -1.364  1.00  0.00           H    1.008   0.0380 HC
ATOM     20  HB1 ALA A   5       1.508  -1.833  -1.273  1.00  0.00           H    1.008   0.0380 HC
ATOM     21  HB2 ALA A   5       1.527  -0.324  -2.215  1.00  0.00           H    1.008   0.0380 HC
ATOM     14  C   ALA A   5       2.016   1.389   0.007  1.00  0.00           C   12.011   0.6160 C
ATOM     15  O   ALA A   5       1.283   2.383   0.066  1.00  0.00           O   15.999  -0.5040 O
DUMM     22  N   ALA A   6       3.490   1.503  -0.019  1.00  0.00           N   14.007  -0.4630 N
"""
# Parse the data
atoms1, atoms2 = parse_pdb_data(data)

#find the best fit rotation matrix and the translation vector
R, t = best_fit_transform(atoms1, atoms2)

#get the euler angles
phi, theta, psi = euler_angles_from_rotation_matrix(R)

print("Rotation Matrix:\n", R)
print("Translation Vector:\n", t)
print("\nEuler Angles (Z-X-Z convention):\n",
      "phi (z):", np.degrees(phi), "degrees\n",
      "theta (x):", np.degrees(theta), "degrees\n",
      "psi (z):", np.degrees(psi), "degrees")

#Verify results
for atom_name, (atom_type, res_name, res_id, coords1) in atoms1.items():
    if atom_name in atoms2:
        atom_type2, res_name2, res_id2, coords2 = atoms2[atom_name]
        predicted_coords = apply_transform(coords1, R, t)
        error = predicted_coords - coords2
        print(f"Atom {atom_name} ({atom_type}):")
        print("  Original Coordinates:", coords1)
        print("  Transformed Coordinates:", coords2)
        print("  Predicted Coordinates:", predicted_coords)
        print("  Error:", error)
        print("-" * 30)
    else:
      print(f"No corresponding atom found for {atom_name}")
