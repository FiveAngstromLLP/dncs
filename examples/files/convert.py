#!/usr/bin/env python3
"""
rotate_lib.py  <infile> <outfile>

Reads a residue-library file, rotates every residue 180° around
the Cα–C bond (N-Cα-C torsion flipped), and writes the new file.
"""

import sys
import numpy as np

def read_records(fname):
    """Return a list of blocks, each block = list of lines."""
    with open(fname) as fh:
        raw = fh.read()

    lines = [ln.rstrip('\n') for ln in raw.splitlines()]

    blocks, cur = [], []
    for ln in lines:
        if ln.startswith(('ATOM', 'HETATM', 'DUMM')) or not ln.strip():
            cur.append(ln)
        else:  # header line like "ALA A 11"
            if cur:
                blocks.append(cur)
            cur = [ln]
    if cur:
        blocks.append(cur)
    return blocks


def parse_atom_line(ln):
    """Return (x,y,z) in Å as float numpy array."""
    # PDB-like fixed columns
    x = float(ln[30:38])
    y = float(ln[38:46])
    z = float(ln[46:54])
    return np.array([x, y, z])


def write_atom_line(ln, xyz):
    """Replace x,y,z and return new line string."""
    # keep formatting: 30-38, 38-46, 46-54
    new_x = f"{xyz[0]:8.3f}"
    new_y = f"{xyz[1]:8.3f}"
    new_z = f"{xyz[2]:8.3f}"
    new_ln = ln[:30] + new_x + new_y + new_z + ln[54:]
    return new_ln

def rotate_block(block):
    """
    Flip the residue 180° around the N→Cα bond while keeping
    the N atom exactly at (0,0,0).
    """
    atom_lines = [ln for ln in block if ln.startswith(('ATOM','HETATM','DUMM'))]
    if not atom_lines:
        return block

    # build a dict of atomic positions
    atoms = {}
    for ln in atom_lines:
        aname = ln[12:16].strip()
        atoms[aname] = parse_atom_line(ln)

    if 'N' not in atoms or 'CA' not in atoms:
        return block

    n_xyz  = atoms['N']          # should be (0,0,0) already
    ca_xyz = atoms['CA']

    # rotation axis N→CA
    axis = ca_xyz - n_xyz
    axis = axis / np.linalg.norm(axis)

    # 180° rotation matrix around this axis through N
    u = axis.reshape(3, 1)
    R = 2 * (u @ u.T) - np.eye(3)

    new_block = []
    for ln in block:
        if ln.startswith(('ATOM','HETATM','DUMM')):
            pos = parse_atom_line(ln)
            aname = ln[12:16].strip()
            if aname == 'N':
                # keep N fixed
                new_block.append(ln)
            else:
                # rotate around N
                new_pos = n_xyz + R @ (pos - n_xyz)
                new_block.append(write_atom_line(ln, new_pos))
        else:
            new_block.append(ln)
    return new_block

def main():
    if len(sys.argv) != 3:
        print("Usage: rotate_lib.py <infile> <outfile>")
        sys.exit(1)
    infile, outfile = sys.argv[1], sys.argv[2]

    blocks = read_records(infile)
    new_blocks = [rotate_block(b) for b in blocks]

    with open(outfile, 'w') as fh:
        for b in new_blocks:
            for ln in b:
                print(ln, file=fh)
            # print('', file=fh)   # blank line between blocks, keeps format

if __name__ == '__main__':
    main()