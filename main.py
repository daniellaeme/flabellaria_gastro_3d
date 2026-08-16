'''
This master script coordinates Stage 1 to Stage 4.
It enforces directory structures, processes geometries, executes Vina,
and generates the visual PDB complex automatically.
'''


import os
import csv
from scripts.ligand_prep import prep_ligand
from scripts.receptor_prep import prep_receptor
from scripts.docking import prep_ligand_pdbqt, prep_receptor_pdbqt, run_vina_docking
from scripts.generate_complex import extract_best_pose_and_merge


if __name__ == '__main__':

    # 1
    # INITIALISE WORKSPACE STRUCTURED SUBDIRECTORIES
    dir_ligand_prep = os.path.join('.', 'data', 'processed', 'ligand_prep')
    dir_receptor_prep = os.path.join('.', 'data', 'processed', 'receptor_prep')
    dir_docking_prep = os.path.join('.', 'data', 'processed', 'docking')
    dir_results = os.path.join('.', 'results')

    # Enforce directory creation and that it matches active workspace exactly
    for directory in [dir_ligand_prep, dir_receptor_prep, dir_docking_prep, dir_results]:
        os.makedirs(directory, exist_ok=True)


    # 2
    # RECEPTOR PREPARATION
    print('/nPREPARING RECEPTOR (COX-2)')
    raw_pdb_path = os.path.join('.', 'data', 'raw', '6COX.pdb')
    protein_name = '6COX'
    native_ligand = 'S58'
    cofactor_name = 'HEM'
    clean_receptor_dest = os.path.join(dir_receptor_prep, 'receptor_clean.pdb')
    rec_pdbqt_dest = os.path.join(dir_docking_prep, 'receptor_prepared.pdbqt')

    if not os.path.exists(raw_pdb_path):
        print(f'\nError: "{raw_pdb_path}" not found! Please download it from rcsb.org')
        exit(1)

    try:
        # Process and write clean receptor straight to receptor_prep subdirectory
        clean_path, (cx, cy, cz) = prep_receptor(
            protein_name=protein_name,
            file_path=raw_pdb_path,
            target_ligand_name=native_ligand,
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
        exit(1)


    # 3
    # COMPOUND LIBRARY
    compounds = {
        # DANIELLI OLIVERI (Rolfe) Hutch. & Dalziel (Fabaceae)
        # Source Literature: Anti-inflammatory and antinociceptive activities of Daniellia oliveri (Fabaceace)
        # stem bark extract (Sofidiya et al., 2023, Journal of Ethnopharmacology)
        # Extraction/Identification Method: Methanolic stem bark extract, identified and quantified via HPLC-DAD-UV.

        # Phenolic Acids
        'Caffeic acid': 'C1=CC(=C(C=C1/C=C/C(=O)O)O)O',
        'p-Coumaric acid': 'C1=CC(=CC=C1/C=C/C(=O)O)O',
        'Ferulic acid': 'COC1=C(C=CC(=C1)/C=C/C(=O)O)O',

        # Flavonoids & Flavonoid Glycosides
        'Rutin': 'C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC(=C(C=C5)O)O)O)O)O)O)O)O',
        'Apigenin-7-glucoside': 'C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)O)O',
        'Quercetin': 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O',
        'Kaempferol': 'C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O',


        # FLABELLARIA PANICULATA Cav. (Malpighiaceae)
        # Source Literature: Gastroprotective Effect and Chemical Constituents of Flabellaria paniculata (Malpighiaceae)
        # (Sofidiya et al.)
        # Extraction/Isolation Method: Column chromatography of the Ethyl Acetate (EtOAc) fraction of methanolic leaf
        # extract, elucidated via NMR and Mass Spectrometry.

        # Triterpenoids:
        'Friedelin': 'C[C@H]1C(=O)CC[C@@H]2[C@@]1(CC[C@H]3[C@]2(CC[C@@]4([C@@]3(CC[C@@]5([C@H]4CC(CC5)(C)C)C)C)C)C)C',
        'Friedelinol': 'CC1C(CCC2C1(CCC3C2(CCC4(C3(CCC5(C4CC(CC5)(C)C)C)C)C)C)C)O',

        # Steroids
        'Beta-Sitosterol': 'CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C',
        'b-Sitosterol b-D-glucoside': 'CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O)C)C)C(C)C',
        # AKA Sitosterol-β-D-glucoside

        # Flavonoid Glycosides
        'Kaempferol-3-O-rutinoside': 'C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC=C(C=C5)O)O)O)O)O)O)O',

        # CONTROLS
        'SC-558': 'C1=CC(=CC=C1C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F)Br',  # Native 6COX Ligand
        'Celecoxib': 'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F',  # Selective COX-2
        'Indomethacin': 'CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=CC(=C2)OC)CC(=O)O',  # Non-selective NSAID
    }


    # 4
    # LOADS CSV CHECKPOINTS, IF THEY EXIST
    results_summary = {}
    csv_output_path = os.path.join(dir_results, 'master_docking_scores.csv')

    if os.path.exists(csv_output_path):
        print('Checkpoint file found. Loading prior docking results...')
        with open(csv_output_path, mode='r', encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader, None)      # Skip the header row
            for row in reader:
                if len(row) == 2:
                    name, score = row
                    try:
                        results_summary[name] = float(score)
                    except ValueError:
                        results_summary[name] = score


    # 5
    # HIGH-THROUGHPUT VIRTUAL SCREENING LOOP
    for name, smiles in compounds.items():
        if not smiles:
            continue

        # Checkpoint: Skip if the compound has already been docked
        if name in results_summary and results_summary[name] != 'FAILED':
            print(f'{name} already processed. Score: {results_summary[name]}\nSkipping...')
            continue

        print(f'\nProcessing Compound: {name.upper()}')

        safe_name = name.lower().replace(' ', '_')

        # Define precise inputs and destinations within respective processed folders
        sdf_output = os.path.join(dir_ligand_prep, f'{safe_name}_best.sdf')
        lig_pdbqt_dest = os.path.join(dir_docking_prep, f'{safe_name}_prepared.pdbqt')
        docked_poses_dest = os.path.join(dir_docking_prep, f'{safe_name}_docked_poses.pdbqt')
        docking_log_dest = os.path.join(dir_results, f'{safe_name}_vina.log')
        merged_complex_dest = os.path.join(dir_results, f'{safe_name}_{protein_name}_complex.pdb')

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
            results_summary[name] = best_score
            print(f'\nSimulation successful. Top binding score: {best_score:.2f} kcal/mol')
            print(f'Results stored under: {dir_docking_prep}')

            # D: Complex Generation
            extract_best_pose_and_merge(
                receptor_path=clean_path,
                docked_pdbqt_path=docked_poses_dest,
                output_complex_pdb=merged_complex_dest
            )
            print(f'Operational Complex saved directly to: {merged_complex_dest}')
            print('Standalone PML script ready. Load via PyMOL command line: @style_complex.pml')

        except Exception as e:
            print(f'\nPipeline failed for {name}: {str(e)}')
            results_summary[name] = 'FAILED'


    # 6
    # FINAL REPORT
    print('\nBATCH DOCKING RESULTS SUMMARY\n')
    sorted_results = sorted(
        [(k, v) for k, v in results_summary.items() if isinstance(v, float)],
        key=lambda item: item[1]
    )

    # Print to terminal
    for name, score in sorted_results:
        print(f'{name.ljust(20)} : {score:.2f} kcal/mol')

    for name, score in results_summary.items():
        if score == 'FAILED':
            print(f'{name.ljust(20)} : {score}')

    # Sort everything and rewrite to the CSV file
    with open(csv_output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Compound', 'Affinity_kcal_mol'])

        for name, score in sorted_results:
            writer.writerow([name, score])

        for name, score in results_summary.items():
            if score == 'FAILED':
                writer.writerow([name , score])

    print(f'Master results table automatically updated at: {csv_output_path}')