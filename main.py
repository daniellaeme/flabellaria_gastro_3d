'''
This master script coordinates Stage 1 to Stage 4.
It enforces directory structures, processes geometries, executes Vina,
and generates the visual PDB complex automatically.
'''


import os
import csv

import pandas as pd

from scripts.ligand_prep import prep_ligand
from scripts.receptor_prep import prep_receptor
from scripts.docking import prep_ligand_pdbqt, prep_receptor_pdbqt, run_vina_docking
from scripts.generate_complex import extract_best_pose_and_merge


if __name__ == '__main__':

    # 1
    # INITIALISE WORKSPACE STRUCTURED SUBDIRECTORIES
    dir_raw = os.path.join('.', 'data', 'raw')
    dir_ligand_prep = os.path.join('.', 'data', 'processed', 'ligand_prep')
    dir_receptor_prep = os.path.join('.', 'data', 'processed', 'receptor_prep')
    dir_docking_prep = os.path.join('.', 'data', 'processed', 'docking')
    dir_results = os.path.join('.', 'results')

    # Enforce directory creation and that it matches active workspace exactly
    for directory in [dir_ligand_prep, dir_receptor_prep, dir_docking_prep, dir_results]:
        os.makedirs(directory, exist_ok=True)

    # 2
    print('Loading databases...')
    targets_path = os.path.join(dir_raw, 'targets.csv')
    compounds_path = os.path.join(dir_raw, 'compounds.csv')

    # Read the CSVs into Pandas DataFrames
    df_targets = pd.read_csv(targets_path)
    df_compounds = pd.read_csv(compounds_path)

    # 3 LOAD MATRIX CHECKPOINT IF IT EXISTS
    csv_output_path = os.path.join(dir_results, 'polypharmacology_matrix.csv')

    if os.path.exists(csv_output_path):
        print('\nFound existing polypharmacology matrix. Loading checkpoints...')
        # Load the existing matrix, using the Compound names as the row index
        results_matrix = pd.read_csv(csv_output_path, index_col=0)

        # Ensure any newly added targets or compounds in your raw CSVs are added to the matrix
        for comp in df_compounds['Compound_Name']:
            if comp not in results_matrix.index:
                results_matrix.loc[comp] = pd.Series(dtype='float64')
        for targ in df_targets['Target_Name']:
            if targ not in results_matrix.columns:
                results_matrix[targ] = pd.NA
    else:
        # Create a blank table if it's the first run
        #                 PROTON_PUMP | UREASE | COX_1
        # Friedelin     |     NaN     |   NaN  |  NaN
        results_matrix = pd.DataFrame(index=df_compounds['Compound_Name'], columns=df_targets['Target_Name'])

        # 4 OUTER LOOP (For the Targets)
        # iterrows() goes through the targets.csv line by line
        for index, row in df_targets.iterrows():
            target_name = row['Target_Name']
            pdb_file = row['PDB_ID']
            native_ligand = row['Native_Ligand']

            # Handle cases where there is no cofactor (eg Proton pump)
            cofactor = str(row['Cofactor'])
            if cofactor.lower() in ['none', 'nan', '']:
                keep_cofactor = False
                cofactor = 'NONE'
            else:
                keep_cofactor = True

            print(f'\n{'=' * 50}')
            print(f'INITIATING PIPELINE FOR TARGET: {target_name}')
            print(f'\n{'=' * 50}')


            # 4-A
            # RECEPTOR PREPARATION
            print(f'/nPREPARING RECEPTOR: {target_name}')
            raw_pdb_path = os.path.join('.', 'data', 'raw', pdb_file)
            clean_receptor_dest = os.path.join(dir_receptor_prep, f'{target_name}_clean.pdb')
            rec_pdbqt_dest = os.path.join(dir_docking_prep, f'{target_name}_prepared.pdbqt')

            if not os.path.exists(raw_pdb_path):
                print(f'\nError: "{raw_pdb_path}" not found! Please download it from rcsb.org to data/raw folder')
                continue # Skip to the next target instead of crashing the whole script

            try:
                # Process and write clean receptor straight to receptor_prep subdirectory
                clean_path, (cx, cy, cz) = prep_receptor(
                    protein_name=target_name,
                    file_path=raw_pdb_path,
                    target_ligand_name=native_ligand,
                    cofactor_name=cofactor,
                    output_path=clean_receptor_dest
                )
                # Convert to PDBQT
                rec_pdbqt = prep_receptor_pdbqt(clean_path, rec_pdbqt_dest)
                print('\nSuccessfully prepared receptor!')
                print(f'Receptor File: {clean_path}')
                print(f'AutoDock Vina Grid Parameters:')
                print(f'  --center_x {cx:.3f}')
                print(f'  --center_y {cy:.3f}')
                print(f'  --center_z {cz:.3f}')
                print(f'  --size_x 22.0 --size_y 22.0 --size_z 22.0')
            except Exception as e:
                print(f'\nReceptor Prep Failed: {str(e)}')
                continue  # Skip to the next target instead of crashing the whole program


            # 4-B
            # THE INNER LOOP (Iterating through Compounds for this Receptor)
            for index, comp in df_compounds.iterrows():
                comp_name = comp['Compound_Name']
                smiles = comp['SMILES']
                safe_comp_name = comp_name.lower().replace(' ', '_')

                # 2D CHECKPOINT
                # Check the exact cell in the Pandas Matrix for this compound and target
                existing_score = results_matrix.at[comp_name, target_name]
                if pd.notna(existing_score) and existing_score != 'FAILED':
                    print(f'\n--- Skipping {comp_name.upper()} into {target_name.upper()} (Already docked: {existing_score} kcal/mol) ---')
                    continue

                print(f'\n--- Docking {comp_name.upper()} into {target_name.upper()} ---')

                # Define precise inputs and destinations within respective processed folders
                sdf_output = os.path.join(dir_ligand_prep, f'{safe_comp_name}_best.sdf')
                lig_pdbqt_dest = os.path.join(dir_docking_prep, f'{safe_comp_name}_prepared.pdbqt')
                docked_poses_dest = os.path.join(dir_docking_prep, f'{safe_comp_name}_{target_name}_poses.pdbqt')
                docking_log_dest = os.path.join(dir_results, f'{safe_comp_name}_{target_name}_vina.log')
                merged_complex_dest = os.path.join(dir_results, f'{safe_comp_name}_{target_name}_complex.pdb')

                try:
                    # A: Ligand 3D Prep
                    prep_ligand(smiles, output_path=sdf_output)
                    print('Ligand prepared...')

                    # B: Ligand PDBQT Conversion
                    lig_pdbqt = prep_ligand_pdbqt(sdf_output, lig_pdbqt_dest)
                    print('Ligand converted to PDBQT...')

                    # C: Docking
                    # Both outputs (docked poses and vina logs) are kept inside docking_prep
                    best_score = run_vina_docking(
                        receptor_pdbqt=rec_pdbqt_dest,
                        ligand_pdbqt=lig_pdbqt,
                        center_coords=[cx, cy, cz],
                        output_poses_path=docked_poses_dest,
                        log_path=docking_log_dest
                    )
                    # Proper pandas assignment to avoid chained indexing warnings
                    results_matrix.loc[comp_name, target_name] = best_score
                    print(f'\nSimulation successful. Top binding score: {best_score:.2f} kcal/mol')
                    print(f'Results stored under: {dir_docking_prep}')

                    # D: Complex Generation
                    extract_best_pose_and_merge(
                        receptor_path=clean_path,
                        docked_pdbqt_path=docked_poses_dest,
                        output_complex_pdb=merged_complex_dest
                    )
                    print(f'Operational Complex saved directly to: {merged_complex_dest}')
                    # print('Standalone PML script ready. Load via PyMOL command line: @style_complex.pml')

                except Exception as e:
                    print(f'\nPipeline failed for {comp_name} into {target_name}: {str(e)}')
                    results_matrix.loc[comp_name, target_name] = 'FAILED'


    # 5
    # FINAL REPORT

    print('\n=========================================')
    print('ALL SIMULATIONS COMPLETE. EXPORTING MATRIX.')
    print('===========================================')

    # Convert the hidden index into a visible column ie 'Compound'
    df_results = results_matrix.copy()
    df_results.index.name = 'Compound'
    df_results.reset_index(inplace=True)

    print('\nBATCH DOCKING RESULTS SUMMARY\n')
    print(df_results.to_string(index=False))

    df_results.to_csv(csv_output_path, index=False)
    print(f'Master Polypharmacology Matrix automatically updated at: {csv_output_path}')