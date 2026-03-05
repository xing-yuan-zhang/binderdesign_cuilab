receptor_chain = "A"
binder_chain   = "B"

iface_cutoff        = 5.0
min_iface_contacts  = 40

hotspots = ["A436", "A176", "A177", "A178", "A179", "A253", "A254", "A255"]
hotspot_cutoff          = 8.0
min_hotspots_contacted  = 3

min_allowed_dist = 2.0

binder_len_min      = 70
binder_len_max      = 100
Rg_min              = 7.0
Rg_max              = 20.0
bb_contact_cutoff   = 8.0
min_internal_CA_contacts   = 80
min_neighbors_per_residue  = 2

use_helix_filter    = True
min_helices         = 3
min_helix_len       = 5

phi_helix_min, phi_helix_max = -120.0, -30.0
psi_helix_min, psi_helix_max = -80.0,  45.0
