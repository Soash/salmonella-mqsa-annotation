# Rename the objects to simpler names (recommended)
set_name rec.pdb, receptor
set_name lig.000.00.pdb, ligand

# Combine them into one complex
create complex, receptor or ligand

# Save the merged structure
save complex.pdb, complex

reinitialize

count_atoms all
count_atoms chain A
count_atoms elem C

remove solvent
remove hetatm

save chainA.pdb, chain A

remove chain C+D

iterate all, print(chain)
save new_model.pdb, all

align model_01, model_02

RMSD between selections
    align chain A and model_01, chain A and model_02




align chain A and complex, chain A and 3hi2

align chain A and complex, chain C and 3hi2

align chain B and complex, chain B and 3hi2
