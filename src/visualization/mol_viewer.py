from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import json

class MoleculeViewer:
    @staticmethod
    def generate_3d_coords(smiles):
        """Generate 3D coordinates for a molecule from SMILES."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Generate 3D coordinates
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            return mol
        except:
            return None
    
    @staticmethod
    def get_viewer_html(mol, size=(400, 400), style="stick"):
        """Generate HTML for an interactive 3D viewer."""
        if mol is None:
            return None
            
        # Convert molecule to PDB format
        pdb = Chem.MolToPDBBlock(mol)
        
        # Create a custom HTML with embedded JavaScript
        html = f"""
        <div id="viewport" style="width: {size[0]}px; height: {size[1]}px; position: relative;"></div>
        <script>
            var viewer = $3Dmol.createViewer("viewport", {{
                backgroundColor: "black"
            }});
            
            var pdb = `{pdb}`;
            viewer.addModel(pdb, "pdb");
            
            viewer.setStyle({{}}, {{"stick": {{"radius": 0.2}}}});
            viewer.addStyle({{}}, {{"sphere": {{"radius": 0.5}}}});
            
            viewer.zoomTo();
            viewer.render();
        </script>
        """
        
        # Add py3Dmol dependencies
        full_html = f"""
        <script src="https://3dmol.org/build/3Dmol-min.js"></script>
        <script src="https://3dmol.org/build/3Dmol.ui-min.js"></script>
        {html}
        """
        
        return full_html
    
    @staticmethod
    def get_viewer_component(mol, size=(400, 400), style="stick"):
        """Generate a py3Dmol viewer component."""
        if mol is None:
            return None
        
        viewer = py3Dmol.view(width=size[0], height=size[1])
        pdb = Chem.MolToPDBBlock(mol)
        viewer.addModel(pdb, "pdb")
        viewer.setStyle({}, {"stick": {"radius": 0.2}})
        viewer.addStyle({}, {"sphere": {"radius": 0.5}})
        viewer.setBackgroundColor("black")
        viewer.zoomTo()
        return viewer
