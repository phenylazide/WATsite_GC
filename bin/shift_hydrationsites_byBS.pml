load ref.pdb
load lig.pdb
load md_HydrationSites.pdb
select bs_atm, ref w. 5 of lig
select bs_res, byres bs_atm
save bs_res_ref_HydrationSites.pdb, bs_res or md_HydrationSites

load prot_only.pdb
load MyBindingSite.pdb
select bs_atm, prot_only w. 5 of MyBindingSite
select bs_res, byres bs_atm
remove bs_res and resn SOL
save bs_res.pdb, bs_res

load bs_res_ref_HydrationSites.pdb
load bs_res.pdb
align bs_res_ref_HydrationSites, bs_res
save HydrationSites_byBS.pdb, bs_res_ref_HydrationSites and resn WAT 

