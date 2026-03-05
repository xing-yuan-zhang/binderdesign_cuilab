def parse_pdb_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append({
                    "line": line.rstrip("\n"),
                    "rec":  line[0:6],
                    "serial": int(line[6:11]),
                    "atom": line[12:16],
                    "resn": line[17:20],
                    "chain": line[21],
                    "resi": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                })
    return atoms

def write_pdb_lines(path, lines):
    import os
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for l in lines:
            f.write(l + "\n")
