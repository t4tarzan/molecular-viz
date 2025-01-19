import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def create_methane():
    # Create methane with explicit hydrogens
    mol = Chem.MolFromSmiles("[CH4]")
    mol = Chem.AddHs(mol)  # Add explicit hydrogens
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, randomSeed=42)  # Use seed for reproducibility
    AllChem.MMFFOptimizeMolecule(mol)  # Optimize the geometry
    
    return mol

def get_viewer_html(mol, size=(400, 400)):
    """Generate HTML for the 3D viewer with debug info."""
    if mol is None:
        return None
        
    # Convert molecule to PDB format
    pdb = Chem.MolToPDBBlock(mol)
    
    # Print debug info
    st.text("Number of atoms: " + str(mol.GetNumAtoms()))
    st.text("Atom details:")
    for atom in mol.GetAtoms():
        st.text(f"Atom {atom.GetIdx()}: {atom.GetSymbol()} with {atom.GetNumExplicitHs()} explicit H")
    
    html = f"""
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <div id="viewport" style="width: {size[0]}px; height: {size[1]}px; position: relative;"></div>
    <script>
        $(document).ready(function() {{
            let viewer = $3Dmol.createViewer($("#viewport"), {{
                backgroundColor: "black"
            }});
            
            let pdb = `{pdb}`;
            viewer.addModel(pdb, "pdb", {{keepH: true}});  // Explicitly keep hydrogens
            
            // Clear any existing styles
            viewer.setStyle({{}}, {{}});
            
            // Style carbon atoms
            viewer.setStyle({{"elem": "C"}}, {{
                "sphere": {{
                    "radius": 0.8,
                    "color": "gray"
                }},
                "stick": {{
                    "radius": 0.2,
                    "color": "gray"
                }}
            }});
            
            // Style hydrogen atoms - make them very visible
            viewer.setStyle({{"elem": "H"}}, {{
                "sphere": {{
                    "radius": 0.4,
                    "color": "white"
                }},
                "stick": {{
                    "radius": 0.2,
                    "color": "white"
                }}
            }});
            
            // Add labels for debugging
            viewer.addLabel("C", {{"position": {{"x": -0.000, "y": -0.000, "z": -0.000}}, "backgroundColor": "black", "fontColor": "white"}});
            viewer.addLabel("H", {{"position": {{"x": -0.675, "y": 0.854, "z": -0.085}}, "backgroundColor": "black", "fontColor": "white"}});
            viewer.addLabel("H", {{"position": {{"x": -0.395, "y": -0.836, "z": -0.581}}, "backgroundColor": "black", "fontColor": "white"}});
            viewer.addLabel("H", {{"position": {{"x": 0.085, "y": -0.294, "z": 1.049}}, "backgroundColor": "black", "fontColor": "white"}});
            viewer.addLabel("H", {{"position": {{"x": 0.986, "y": 0.276, "z": -0.382}}, "backgroundColor": "black", "fontColor": "white"}});
            
            viewer.zoomTo();
            viewer.render();
        }});
    </script>
    """
    return html

def main():
    st.title("Methane Visualization Test")
    
    # Create methane molecule
    mol = create_methane()
    
    # Show SMILES and molecular formula
    st.text("SMILES: [CH4]")
    st.text("Molecular Formula: CH4")
    
    # Show 3D visualization
    viewer_html = get_viewer_html(mol)
    st.components.v1.html(viewer_html, height=450)
    
    # Show PDB format for debugging
    st.text("PDB Format:")
    st.code(Chem.MolToPDBBlock(mol))

if __name__ == "__main__":
    main()
