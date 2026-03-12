import os
import glob
from pymol import cmd

def align_fragments_pipeline(folder_path, reference_code):

    ref_obj = f"{reference_code}_protein"
    if ref_obj not in cmd.get_names():
        print(f"Fetching reference: {reference_code} as {ref_obj}")
        cmd.fetch(reference_code, name=ref_obj)

    search_pattern = os.path.join(folder_path, "*.sdf")
    sdf_files = glob.glob(search_pattern)

    if not sdf_files:
        print(f"No SDF files found in {folder_path}")
        return

    for sdf_path in sdf_files:
        pdb_code = os.path.splitext(os.path.basename(sdf_path))[0]
        protein_obj_name = f"{pdb_code}_protein"
        fragment_obj_name = f"{pdb_code}_fragments"

        if protein_obj_name in cmd.get_names():
            print(f"Skipping {pdb_code}, {protein_obj_name} already exists.")
            continue

        print(f"--- Processing {pdb_code} ---")

        try:
         
            cmd.fetch(pdb_code, name=protein_obj_name)

  
            alignment = cmd.align(f"{protein_obj_name} and name CA", f"{ref_obj} and name CA")
            rmsd = alignment[0]

            if rmsd > 4.0:
                print(f"!!! RMSD too high ({rmsd:.2f}). Removing {protein_obj_name}.")
                cmd.delete(protein_obj_name)
                continue

 
            matrix = cmd.get_object_matrix(protein_obj_name)

     
            cmd.load(sdf_path, fragment_obj_name)


            num_states = cmd.count_states(fragment_obj_name)
            for state_idx in range(1, num_states + 1):
                cmd.transform_object(fragment_obj_name, matrix, state=state_idx)

            cmd.disable(protein_obj_name)
            cmd.remove(f"{protein_obj_name} and resn DMS")
            print(f"Successfully aligned {num_states} states for {pdb_code}")

        except Exception as e:
            print(f"Error processing {pdb_code}: {e}")

    cmd.zoom(ref_obj)
    print("\nPipeline Complete.")

my_folder = "/home/joan/posal/deliverables/SAFR_v1.0.1/P0DTD1"
ref = "5R8T"
align_fragments_pipeline(my_folder, ref)
