
<img width="1827" height="921" alt="flabellaria_network csv" src="https://github.com/user-attachments/assets/db97355d-e604-4e55-a36b-a8857e1c92d4" />

# The Polypharmacology of Flabellaria paniculata (Cav.)

Computational Validation of Traditional Gastroprotective & Anti-Inflammatory Synergy

Authors: Daniella E. Ene-Obong & Margaret O. Sofidiya

Institution: Department of Pharmacognosy, Faculty of Pharmacy, University of Lagos, Akoka, Nigeria.

To be presented at: 19th UNILAG Annual Research Conference

# Project Overview

African traditional medicine frequently relies on complex, multi-component aqueous decoctions to manage chronic conditions like gastric ulcers. The mechanistic validation of these polyherbal remedies remains a significant challenge. Flabellaria paniculata is extensively utilized in Nigerian traditional medicine to treat inflammatory and gastrointestinal disorders.

The Problem: Traditional synthetic NSAIDs efficiently inhibit COX-2 to treat inflammation, but indiscriminately block COX-1, destroying mucosal defenses and causing severe gastric ulcers. Furthermore, these ulcers are exacerbated by overactive Gastric Proton Pumps and the mucosal-damaging ammonia cloud produced by H. pylori Urease.

**Objective**: This project bridges indigenous knowledge with rational drug design by deploying an automated in silico pipeline (Network Pharmacology, Molecular Docking, and ADMET profiling) to elucidate the multi-targeted synergy of this indigenous remedy.

# Computational Pipeline & Methodology

1. Phytochemical Curation: Phytochemicals identified via existing literature and spectral data were curated, and 3D structures were geometrically optimized (PubChem/PDB).

2. Network Pharmacology: Putative targets were predicted (SwissTargetPrediction) and mapped topologically (Cytoscape v3.9) to uncover the shared genetic axis between systemic inflammation and gastric ulcers.

3. High-Throughput Molecular Docking: 52 rigorous thermodynamic binding simulations were executed via AutoDock Vina against critical targets (COX-1, COX-2, Urease, and Proton Pump).

4. Structural Validation & ADMET: Atom-to-atom noncovalent interactions were mathematically extracted using PLIP and visualized via Mol*/PyMOL. In vivo safety and absorption parameters were confirmed via predictive toxicology (ProTox-II/SwissADME).

# Results

| Compound                        | PROTON_PUMP | UREASE | COX_1  | COX_2  |
|---------------------------------|-------------|--------|--------|--------|
| Friedelin                       | -7.515      | -7.099 | 5.98   | -8.363 |
| Friedelinol                     | -4.384      | -7.125 | 8.696  | -5.241 |
| Beta-Sitosterol                 | -9.618      | -6.571 | -6.079 | -6.907 |
| b-Sitosterol b-D-glucoside      | -8.799      | -6.399 | -6.153 | -7.003 |
| Kaempferol-3-O-rutinoside       | -9.101      | -7.907 | -4.788 | -4.193 |
| Vonoprazan (HKT) - NL           | -9.165      | -6.537 | -6.488 | -7.697 |
| Omeprazole - SD                 | -8.29       | -6.852 | -7.45  | -6.965 |
| Acetohydroxamic acid (HAE) - NL | -4.143      | -4.273 | -3.76  | -3.729 |
| Thiourea                        | -3.244      | -3.478 | -3.157 | -3.151 |
| Ibuprofen (IBF) - NL            | -7.35       | -6.346 | -7.435 | -7.223 |
| Indomethacin - SD               | -9.076      | -7.085 | -7.172 | -8.898 |
| SC-558 (S58) - NL               | -10.09      | -7.627 | -6.211 | -11.34 |
| Celecoxib - SD                  | -10.67      | -7.71  | -6.429 | -11.36 |

<img width="1480" height="550" alt="b_sitosterol_b_d_glucoside_proton_pump_complex_labelled_pymol_image" src="https://github.com/user-attachments/assets/4eeaadb9-6459-4b5a-a37f-7d23567bf741" />

<img width="1480" height="550" alt="beta_sitosterol_proton_pump_complex_labelled_pymol_image" src="https://github.com/user-attachments/assets/30cdf01f-9b7a-497f-98c0-f1877aa05944" />

<img width="1480" height="550" alt="friedelin_COX_2_complex_labelled_pymol_image" src="https://github.com/user-attachments/assets/810dec12-6485-48ba-96ec-90ed5e309c72" />

<img width="1480" height="550" alt="kaempferol_3_o_rutinoside_urease_complex_labelled_pymol_image" src="https://github.com/user-attachments/assets/90df8aa0-0669-4748-a7d0-74f70fba93dc" />

# Key Findings

1. Network Pharmacology (The Shared Genetic Axis)

Topological mapping confirmed the plant's traditional polypharmacology is rooted in a shared genetic axis. Targets including PTGS2, TNF, and IL1B bridge both systemic inflammation and gastric ulcers.

Friedelin anchors the anti-inflammatory arm by inhibiting PTGS2.

Beta-sitosterol and Kaempferol-3-O-rutinoside directly modulate the ulcerogenic microenvironment by targeting the gastric proton pump (ATP4A) and H. pylori survival (Urease).

2. Resolving the NSAID Paradox (Molecular Docking)

The 3D binding thermodynamics computationally resolve the classic NSAID paradox:

Friedelin (COX-2 Selective): The massive triterpenoid scaffold exhibits deep pocket insertion into COX-2 (-8.36 kcal/mol) via hydrophobic enclosure (Val349, Leu352). Structural analysis predicts severe steric clashes with the narrower COX-1 pocket, explaining the decoction's ability to halt inflammation without triggering mucosal toxicity.

Beta-sitosterol (Proton Pump): Exhibits superior mechanical inhibition (-9.62 kcal/mol) of the pump mechanism via robust hydrophobic interactions (Cys813, Tyr799), outperforming Omeprazole's initial non-covalent affinity phase.

Kaempferol-3-O-rutinoside (Urease): Neutralizes the bacterial active site (-7.91 kcal/mol) via strong H-bonds and a critical salt bridge (His322), preventing mucosal-damaging ammonia cloud generation.

3. Formulation Dynamics (The Nanoparticle Hypothesis)

ADMET profiling reveals that the highly lipophilic aglycones (like Friedelin) are likely solubilized in the traditional aqueous decoction via amphiphilic glycosides acting as natural biosurfactants. Exhibiting poor GI absorption (Tox Class IV), these bulky glycosides function computationally as microbiome-cleavable prodrugs. They survive upper GI transit to deliver active, binding-competent aglycones directly to inflamed lower intestinal lesions.

# Repository Structure

/data: Contains raw interaction CSVs, target lists, and SMILES strings.

/scripts: Python scripts utilized for ligand/receptor preparation and Vina automation.

/results: Contains the final 52 docking logs, PLIP interaction tables, and SwissTargetPrediction outputs.

/poster: High-resolution digital copies of the A0 Conference Poster (PDF/PNG).

# Contact & Future Work

Future Directions: To transition these static docking poses into dynamic physiological models, 100ns Molecular Dynamics (MDS) simulations (GROMACS) coupled with MM-PBSA free energy calculations are planned to validate complex stability in solvated environments.

Contact the Authors:

Daniella E. Ene-Obong: https://www.linkedin.com/in/daniella-ene-obong/

Margaret O. Sofidiya

Acknowledgments: The authors acknowledge the University of Lagos for providing the academic and computational platform for this research.
