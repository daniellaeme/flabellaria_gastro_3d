# from rdkit import RDConfig
# print(RDConfig.RDDataDir)
# from rdkit.Chem import Descriptors
# print(len(Descriptors._descList))

from rdkit.Chem import Descriptors, Lipinski, QED
import plotly.graph_objects as go

def get_descriptors(mol):
    """
    Calculates molecular descriptors.
    """
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP" : Descriptors.MolLogP(mol),
        "HBD" : Lipinski.NumHDonors(mol),
        "HBA" : Lipinski.NumHAcceptors(mol),
        "TPSA" : Descriptors.TPSA(mol),
        "RotBonds" : Descriptors.NumRotatableBonds(mol),
        "QED" : QED.qed(mol)
    }

# Check VIolations
def check_violations(metrics):
    """
    Checks the individual descriptors against Lipinski snd Veber rules.
    """
    violations = []

    # Lipinski Rules
    if metrics["MW"] > 500:
        violations.append("MW > 500")
    if metrics["LogP"] > 5:
        violations.append("LogP > 5")
    if metrics["HBD"] > 5:
        violations.append("HBD > 5")
    if metrics["HBA"] > 10:
        violations.append("HBA > 10")

    # Veber Rules
    if metrics["TPSA"] > 140:
        violations.append("TPSA > 140")
    if metrics["RotBonds"] > 10:
        violations.append("RotBonds > 10")

    return violations

def plot_radar_chart(metrics):
    """
    Creates a Swiss-ADME-style radar chart using Plotly.
    Values ae normalised so the 'safe zone' is always uniform.
    """
    # Define upper limits for normalization
    limits = {
        "MW":600, "LogP":6, "HBD":6, "HBA":12, "TPSA":160, "RotBonds":12
    }

    ideal_max = {
        "MW": 500, "LogP": 5, "HBD": 5, "HBA": 10, "TPSA": 140, "RotBonds": 10
    }

    categories = list(limits.keys())

    # Normalise actual values to a 0-1 scale relative to the limits
    # eg if MW is 300 and limit is 600, value is 0.5
    values = [min(metrics[cat]/limits[cat], 1.0) for cat in categories]

    # Normalise the Safe Zone threshold
    safe_zone = [ideal_max[cat]/limits[cat] for cat in categories]

    # Close the loop for the plot
    categories += [categories[0]]
    values += [values[0]]
    safe_zone += [safe_zone[0]]

    fig = go.Figure()

    # Draw the Safe Zone (pink area)
    fig.add_trace(go.Scatterpolar(
        r=safe_zone,
        theta=categories,
        fill="toself",
        name="Safe Limit",
        line=dict(color="rgba(255, 99, 71, 0.5)"),
        fillcolor="rgba(255, 99, 71, 0.2)"
    ))

    # Draw the Molecule  (blue line)
    fig.add_trace(go.Scatterpolar(
        r=values,
        theta=categories,
        fill="toself",
        name="Molecule",
        line=dict(color="blue", width=2)
    ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(visible=True, range=[0,1], showticklabels=False)
        ),
        showlegend=True,
        margin=dict(l=40, r=40, t=40, b=40)
    )

    return fig
