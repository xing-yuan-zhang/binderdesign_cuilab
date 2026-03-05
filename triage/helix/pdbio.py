import numpy as np

def parse_pdb_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atoms.append({
                        "line": line.rstrip("\n"),
                        "atom": line[12:16],
                        "resn": line[17:20],
                        "chain": line[21],
                        "resi": int(line[22:26]),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                    })
                except ValueError:
                    continue
    return atoms

def atoms_to_array(atoms):
    if not atoms:
        return np.zeros((0, 3), dtype=float)
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)

def get_chain_atoms(atoms, chain_id):
    return [a for a in atoms if a["chain"] == chain_id]

def get_chain_CA_atoms(atoms, chain_id):
    return [a for a in atoms if a["chain"] == chain_id and a["atom"].strip() == "CA"]

def get_hotspot_atoms(atoms, hotspot_res):
    chain = hotspot_res[0]
    resi  = int(hotspot_res[1:])
    return [a for a in atoms if a["chain"] == chain and a["resi"] == resi]
