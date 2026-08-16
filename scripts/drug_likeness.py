import os
import urllib.request
import sys
from rdkit import Chem
import pandas as pd

urls = {'utils.py': 'https://raw.githubusercontent.com/daniellaeme/My-Chem-Bioinformatics-Toolkit/refs/heads/main/%5B2%5D_drug-likeness-predictor/utils.py'}
for filename, url in urls.items():
    if os.path.exists(filename):
        os.remove(filename)
    urllib.request.urlretrieve(url, filename)
    print(f'Downloaded {filename}')

sys.path.append('.')

import utils

def smiles_to_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return utils.get_descriptors(mol)
    return {"MW":None, "LogP":None, "HBD":None, "HBA":None, "TPSA":None, "RotBonds":None, "QED":None}


def count_lipinski_violations(row):
    """Counts how many of Lipinski's Rule of 5 are broken"""
    violations = 0
    if pd.notna(row['MW']) and row['MW'] > 500: violations += 1
    if pd.notna(row['LogP']) and row['LogP'] > 5: violations += 1
    if pd.notna(row['HBD']) and row['HBD'] > 5: violations += 1
    if pd.notna(row['HBA']) and row['HBA'] > 10: violations += 1
    return violations


xlsx_input = '../data/raw/compound_library.xlsx'
if not os.path.exists(xlsx_input):
    print(f'Error: file {xlsx_input} does not exist')
    sys.exit()

df = pd.read_excel(xlsx_input)

print("Calculating descriptors...")
descriptors_df = df['canonical_smiles'].apply(smiles_to_descriptors).apply(pd.Series)
dl_df = pd.concat([df, descriptors_df], axis=1)
dl_df['Lipinski_Violations'] = dl_df.apply(count_lipinski_violations, axis=1)

print("\nResults:")
print(dl_df[['canonical_smiles', 'MW', 'LogP', 'HBD', 'HBA', 'Lipinski_Violations']].head())

os.makedirs('../adme_t_drug_likeness', exist_ok=True)
output_path = '../adme_t_drug_likeness/compound_drug_likeness.csv'
dl_df.to_csv(output_path, index=False)
print(f'\nDrug-Likeness profiles saved to: {output_path}')
