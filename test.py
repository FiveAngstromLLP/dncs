# %%

pdb_data = """ATOM      1 C1   NAG B 401      22.391  48.976  10.140                       C
ATOM      2 C2   NAG B 401      21.289  48.570  11.119                       C
ATOM      3 C3   NAG B 401      19.913  48.893  10.538                       C
ATOM      4 C4   NAG B 401      19.837  50.356  10.125                       C
ATOM      5 C5   NAG B 401      20.974  50.682   9.162                       C
ATOM      6 C6   NAG B 401      21.024  52.143   8.777                       C
ATOM      7 C7   NAG B 401      21.933  46.712  12.588                       C
ATOM      8 C8   NAG B 401      21.948  45.225  12.772                       C
ATOM      9 N2   NAG B 401      21.385  47.157  11.452                       N
ATOM     10 O3   NAG B 401      18.918  48.627  11.518                       O
ATOM     11 O4   NAG B 401      18.589  50.627   9.497                       O
ATOM     12 O5   NAG B 401      22.233  50.368   9.775                       O
ATOM     13 O6   NAG B 401      21.880  52.353   7.661                       O
ATOM     14 O7   NAG B 401      22.397  47.475  13.429                       O
"""

def subtract_coordinates(pdb_string):
    """Subtracts the coordinates of C1 from all other coordinates in a PDB string.

    Args:
      pdb_string: A string containing PDB data.

    Returns:
      A string containing the modified PDB data with subtracted coordinates.
    """

    lines = pdb_string.strip().split("\n")
    c1_coords = []
    modified_lines = []


    for line in lines:
        if line.startswith("ATOM") and line[12:17].strip() == "C1":
            c1_coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            modified_line = f"{line[0:30]}{0.000:8.3f}{0.000:8.3f}{0.000:8.3f}{line[54:]}"
            modified_lines.append(modified_line)
            continue


        if line.startswith("ATOM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            new_x = x - c1_coords[0]
            new_y = y - c1_coords[1]
            new_z = z - c1_coords[2]
            modified_line = f"{line[0:30]}{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}{line[54:]}"
            modified_lines.append(modified_line)
        else:
          modified_lines.append(line)

    return "\n".join(modified_lines)


modified_pdb = subtract_coordinates(pdb_data)
print(modified_pdb)
