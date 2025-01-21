#!/usr/bin/env python3
import math

def rotate_coordinates(x, y, z, x_deg, y_deg, z_deg):
    """
    Rotate a single point (x, y, z) about the x, y, and z axes by x_deg, y_deg, z_deg (in degrees).
    Returns the rotated coordinates (rx, ry, rz).
    """
    # Convert degrees to radians
    x_rad = x_deg * math.pi / 180.0
    y_rad = y_deg * math.pi / 180.0
    z_rad = z_deg * math.pi / 180.0

    # --- Rotate about X ---
    # x stays the same
    # y, z transform
    x1 = x
    y1 = y*math.cos(x_rad) - z*math.sin(x_rad)
    z1 = y*math.sin(x_rad) + z*math.cos(x_rad)

    # --- Rotate about Y ---
    # y stays the same
    # x1, z1 transform
    x2 = x1*math.cos(y_rad) + z1*math.sin(y_rad)
    y2 = y1
    z2 = -x1*math.sin(y_rad) + z1*math.cos(y_rad)

    # --- Rotate about Z ---
    # z stays the same
    # x2, y2 transform
    rx = x2*math.cos(z_rad) - y2*math.sin(z_rad)
    ry = x2*math.sin(z_rad) + y2*math.cos(z_rad)
    rz = z2

    return rx, ry, rz

def rotate_pdb(infile, outfile):
    """
    Reads a PDB file (infile), asks the user whether to rotate the coordinates or not,
    and writes a new PDB file (outfile) with either rotated or unchanged coordinates.
    """
    # Ask user if we should rotate or not
    inarg = int(input("Enter 1 to rotate coords, or 0 to leave them unchanged: "))

    if inarg == 1:
        x_rotang = float(input("Enter the angle to rotate about x (degrees): "))
        y_rotang = float(input("Enter the angle to rotate about y (degrees): "))
        z_rotang = float(input("Enter the angle to rotate about z (degrees): "))
    else:
        # No rotation
        x_rotang = 0.0
        y_rotang = 0.0
        z_rotang = 0.0

    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for line in fin:
            # Only rotate lines that start with ATOM or HETATM
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Parse x, y, z columns (PDB columns 31–38, 39–46, 47–54 in 1-based indexing)
                x_str = line[30:38]
                y_str = line[38:46]
                z_str = line[46:54]

                # Convert to float
                x_val = float(x_str)
                y_val = float(y_str)
                z_val = float(z_str)

                # Rotate coordinates
                rx, ry, rz = rotate_coordinates(x_val, y_val, z_val,
                                                x_rotang, y_rotang, z_rotang)

                # Format the new coordinates back into the line
                # We'll keep the rest of the line intact from column 0 to 30 and 54 onward,
                # and only overwrite columns 31-54 with new coordinates.
                new_line = (
                    line[:30] +
                    f"{rx:8.3f}" +
                    f"{ry:8.3f}" +
                    f"{rz:8.3f}" +
                    line[54:]
                )
                fout.write(new_line)
            else:
                # For other lines (TER, END, REMARK, etc.), just copy
                fout.write(line)

if __name__ == "__main__":
    # Example usage
    input_pdb = "test.pdb"
    output_pdb = "rotated_output.pdb"

    rotate_pdb(input_pdb, output_pdb)
    print(f"Finished. Rotated PDB written to {output_pdb}")
