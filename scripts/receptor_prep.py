import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select


def prep_receptor(protein_name: str, file_path: str, target_ligand_name: str, cofactor_name: str,output_path=os.path.join('.', 'data', 'processed', 'receptor_prep', 'receptor_clean.pdb')):
    '''
    Parses, cleans, isolates, and prepares receptor structures.
    Preserves cofactors and outputs directly to receptor_prep subdirectory.
    :param pdb:
    :return:
    '''

    # 1
    # Parse the raw PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein_name', file_path)
    model_0 = structure[0]              # Access Model 0

    print(f"DEBUG: Number of chains found in model: {len(list(model_0))}")
    all_residues = [r.get_resname().strip() for c in model_0 for r in c]
    print(f"DEBUG: Total residues found: {len(all_residues)}")
    print(f"DEBUG: Unique residue names: {set(all_residues)}")

    print('\n[1/5]')
    print('PDB file successfully parsed')

    # 2
    # Dynamically locate the native ligand
    target_residue = None           # Before we strip the native ligand (LVG) away, we must find where it lives, inside the protein structure so we can calculate the center point of our docking box

    for chain in model_0:               # Search Structure -> Model -> Chain -> Residue to find 'LVG'
        for residue in chain:           # Search all chains
            if residue.get_resname().strip() == target_ligand_name:
                target_residue = residue
                host_chain_id = chain.id
                break
        if target_residue:              # Break outer loop once found
            break

    if target_residue is None:          # If not found, stop execution with error
        raise ValueError(f'Error: Target ligand {target_ligand_name} not found in any chain of {protein_name}')

    print('\n[2/5]')
    print(f'Located target ligand "{target_ligand_name}" inside Chain "{host_chain_id}"')

    # 3
    # Calculate the 3D Docking Centroid of the active site
    # Extract the [x, y, z] coordinates of every single atom belonging to 'LVG'
    all_atoms = [atom.get_coord() for atom in target_residue]       # Return 3D np array: [x, y, z]
    all_atoms_np = np.array(all_atoms)          # Converts to 2D np array

    # Compute mean along vertical columns for each coordinate
    center_x = float(all_atoms_np[:, 0].mean())
    center_y = float(all_atoms_np[:, 1].mean())
    center_z = float(all_atoms_np[:, 2].mean())

    print(f'\n[3/5]')
    print(f'Calculated active site Centroid based on native {target_ligand_name}:')
    print(center_x, center_y, center_z)

    # 4
    # Spatial Detection of Interfacial Partner Chain
    chains = list(model_0)
    num_chains = len(chains)                 # Find which other chain forms the binding cleft with the host chain
    partner_chain_id = None

    if num_chains > 1:
        print('\n[4/5]')
        print(f'Multi-chain structure detected ({num_chains} chains). Initiating spatial contact detection...')
        min_distance = float('inf')

        for other_chain in model_0:
            if other_chain.id == host_chain_id:
                continue

            # Collect coordinates of all structural protein atoms in this chain
            other_chain_coords = []
            for residue in other_chain:
                # Skip solvent, cofactors, or other ligands to avoid false proximity readings
                if residue.get_resname().strip() in ['HOH', target_ligand_name]:
                    continue
                for atom in residue:
                    other_chain_coords.append(atom.get_coord())

            if not other_chain_coords:
                continue

            other_chain_coords_np = np.array(other_chain_coords)

            # Calculate minimum distance between the ligand and the chain's atoms
            for lig_coord in all_atoms_np:
                diff = other_chain_coords_np - lig_coord
                dist_sq = np.sum(diff ** 2, axis=1)
                current_min = np.sqrt(np.min(dist_sq))
                if current_min < min_distance:
                    min_distance = current_min
                    partner_chain_id = other_chain.id

        if partner_chain_id is None:
            raise ValueError('Error: A physical partner chain could not be identified or the interface! Multiple chains are present in the structural model, '
                  'but no physical partner chain was close enough to the interfacial cleft.')
        else:
            print(f'Spatially identified interfacial partner chain: Chain "{partner_chain_id}".')
            print(f'Minimum inter-chain contact distance: {min_distance:.2f} Å.')


    else:
        print('\n[4/5]')
        print(f'Operating in Monomer Mode (only 1 chain found: Chain "{host_chain_id}"). '
              f'Skipping interfacial partner detection.')

    # 5
    # Custom Selector to inspect every chain and residue and determine acceptance criteria
    class CustomFilter(Select):
        """
        Custom structural filter that isolates only the functional chains,
        dehydrates the pocket, and purges non-functional crystallization ligands.
        """

        def __init__(self, host, partner, target_ligand, keep_cofactor=True, cofactor_name='None'):
            self.host = host
            self.partner = partner
            self.target_ligand = target_ligand
            self.keep_cofactor = keep_cofactor
            self.cofactor_name = cofactor_name

        def accept_chain(self, chain):                  # If monomeric, keep only the host chain. Otherwise, keep the host-partner dimer.
            return chain.id in [self.host, self.partner]        # Reject Chain C to keep only the functional interfacial dimer (A & B)

        # Blacklisting Method
        # def accept_residue(self, residue):
        #     res_name = residue.get_resname().strip()  # strip? # Read residue name
        #     if res_name in ['HOH', 'NAG', self.target_ligand]:
        #         return False  # Discard water molecules, crystallization lipids, and detergents
        #     if self.cofactor_name != 'None' and res_name == self.cofactor_name:
        #         return self.keep_cofactor  # Retain or remove cofactor only if the protein actually uses one
        #     return True  # Accept standard protein residues (amino acids)

        # Whitelisitng Method
        def accept_residue(self, residue):
            res_name = residue.get_resname().strip()

            # Always keep the assigned cofactor
            if self.cofactor_name != 'None' and res_name == self.cofactor_name:
                return self.keep_cofactor

            # Whitelist standard amino acids
            # We explicitly add 'KCX' because Urease requires this modified Lysine to function
            standard_amino_acids = [
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                'KCX'
            ]

            # Only accept atoms on the whitelist. Automatically delete everything else.
            if res_name in standard_amino_acids:
                return True

            return False


    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    # The 'keep_cofactor=True' default maintains the cofactor-bound state.
    # Set to 'False' if modeling an apo-state docking hypothesis.
    # Set to keep_cof if our CSV says 'None' for cofactor, tell the filter not to look for one
    keep_cof = False if cofactor_name == 'None' else True

    io.save(output_path,
            select=CustomFilter(host=host_chain_id, partner=partner_chain_id, target_ligand=target_ligand_name, keep_cofactor=keep_cof, cofactor_name=cofactor_name)
    )  # Pass custom filter to the select parameter
    print(f'\n[5/5]')
    if partner_chain_id:
        print(f'Interfacial dimer (Chains {host_chain_id} & {partner_chain_id}) successfully isolated and cleaned!')
    else:
        print(f'Monomer (Chain {host_chain_id}) successfully isolated and cleaned!')

    print(f'Saved cleaned structural coordinates to: "{output_path}"')

    return output_path, (center_x, center_y, center_z)



