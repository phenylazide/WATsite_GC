load prot_only.pdb
load MyBindingSite_formatted.pdb
save prot_only_MyBindingSite.pdb, prot_only or MyBindingSite_formatted

load ref.pdb
load prot_only_MyBindingSite.pdb
align prot_only_MyBindingSite, ref
save lig.pdb, prot_only_MyBindingSite and resn BDS
