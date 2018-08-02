load ref.pdb
load md_HydrationSites.pdb
save ref_md_HydrationSites.pdb, ref or md_HydrationSites

load ref_md_HydrationSites.pdb
load prot_only.pdb
align ref_md_HydrationSites, prot_only
save HydrationSites.pdb, ref_md_HydrationSites and resn WAT
