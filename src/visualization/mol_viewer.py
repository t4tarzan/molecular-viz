from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import json

class MoleculeViewer:
    @staticmethod
    def generate_3d_coords(smiles):
        """Generate 3D coordinates for a molecule from SMILES."""
        try:
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Failed to create molecule from SMILES: {smiles}")
                return None
                
            # Add explicit hydrogens
            mol = Chem.AddHs(mol)
            
            # Use ETKDG for better 3D conformer generation
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.enforceChirality = True
            params.useExpTorsionAnglePrefs = True
            params.useBasicKnowledge = True
            params.ETversion = 2
            
            # Generate conformer
            result = AllChem.EmbedMolecule(mol, params)
            if result == -1:
                # Try again with simpler parameters
                params = AllChem.ETKDG()
                result = AllChem.EmbedMolecule(mol, params)
                if result == -1:
                    print("Failed to generate conformer")
                    return None
            
            # Optimize the structure
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except:
                # If MMFF fails, try UFF
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            
            return mol
        except Exception as e:
            print(f"Error generating 3D coordinates: {e}")
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
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://3dmol.org/build/3Dmol-min.js"></script>
        <script src="https://3dmol.org/build/3Dmol.ui-min.js"></script>
        <div id="viewport" style="width: {size[0]}px; height: {size[1]}px; position: relative; background-color: black;"></div>
        <script>
            $(document).ready(function() {{
                let viewer = $3Dmol.createViewer($("#viewport"), {{
                    backgroundColor: "black"
                }});
                
                let pdb = `{pdb}`;
                viewer.addModel(pdb, "pdb", {{keepH: true}});  // Explicitly keep hydrogens
                
                // Clear any existing styles
                viewer.setStyle({{}}, {{}});
                
                if ("{style}" === "stick") {{
                    // Style carbon atoms
                    viewer.setStyle({{"elem": "C"}}, {{
                        "stick": {{
                            "radius": 0.2,
                            "color": "gray"
                        }},
                        "sphere": {{
                            "radius": 0.5,
                            "color": "gray"
                        }}
                    }});
                    
                    // Style hydrogen atoms
                    viewer.setStyle({{"elem": "H"}}, {{
                        "stick": {{
                            "radius": 0.1,
                            "color": "white"
                        }},
                        "sphere": {{
                            "radius": 0.3,
                            "color": "white"
                        }}
                    }});
                    
                    // Style oxygen atoms
                    viewer.setStyle({{"elem": "O"}}, {{
                        "stick": {{
                            "radius": 0.2,
                            "color": "red"
                        }},
                        "sphere": {{
                            "radius": 0.5,
                            "color": "red"
                        }}
                    }});
                    
                    // Style nitrogen atoms
                    viewer.setStyle({{"elem": "N"}}, {{
                        "stick": {{
                            "radius": 0.2,
                            "color": "blue"
                        }},
                        "sphere": {{
                            "radius": 0.5,
                            "color": "blue"
                        }}
                    }});
                    
                }} else if ("{style}" === "sphere") {{
                    // Style carbon atoms
                    viewer.setStyle({{"elem": "C"}}, {{
                        "sphere": {{
                            "radius": 0.8,
                            "color": "gray"
                        }}
                    }});
                    
                    // Style hydrogen atoms
                    viewer.setStyle({{"elem": "H"}}, {{
                        "sphere": {{
                            "radius": 0.4,
                            "color": "white"
                        }}
                    }});
                    
                    // Style oxygen atoms
                    viewer.setStyle({{"elem": "O"}}, {{
                        "sphere": {{
                            "radius": 0.8,
                            "color": "red"
                        }}
                    }});
                    
                    // Style nitrogen atoms
                    viewer.setStyle({{"elem": "N"}}, {{
                        "sphere": {{
                            "radius": 0.8,
                            "color": "blue"
                        }}
                    }});
                    
                }} else if ("{style}" === "structure") {{
                    // Style all atoms with standard colors
                    viewer.setStyle({{}}, {{
                        "stick": {{
                            "radius": 0.2,
                            "colorscheme": "Jmol"
                        }},
                        "sphere": {{
                            "radius": 0.5,
                            "colorscheme": "Jmol"
                        }}
                    }});
                }}
                
                viewer.zoomTo();
                viewer.render();
            }});
        </script>
        """
        return html
    
    @staticmethod
    def get_viewer_component(mol, size=(400, 400), style="stick"):
        """Generate a py3Dmol viewer component."""
        if mol is None:
            return None
        
        viewer = py3Dmol.view(width=size[0], height=size[1])
        pdb = Chem.MolToPDBBlock(mol)
        viewer.addModel(pdb, "pdb")
        
        if style == "stick":
            viewer.setStyle({}, {
                "stick": {"radius": 0.15, "colorscheme": "Jmol"},
                "sphere": {"radius": 0.4, "colorscheme": "Jmol"}
            })
        elif style == "sphere":
            viewer.setStyle({}, {
                "sphere": {"radius": 0.8, "colorscheme": "Jmol"}
            })
        elif style == "structure":
            viewer.setStyle({}, {
                "cartoon": {"color": "spectrum"},
                "stick": {"radius": 0.15, "colorscheme": "Jmol"},
                "sphere": {"radius": 0.4, "colorscheme": "Jmol"}
            })
        
        viewer.setBackgroundColor("black")
        viewer.zoomTo()
        return viewer
