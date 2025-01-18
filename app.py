import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
from pathlib import Path
import sys

# Add the src directory to Python path
sys.path.append(str(Path(__file__).parent))

from src.model.predictor import MolecularPropertyPredictor
from src.visualization.mol_viewer import MoleculeViewer
import py3Dmol

# Page config
st.set_page_config(
    page_title="AI Molecular Property Visualizer",
    page_icon="🧬",
    layout="wide"
)

# Initialize session state
if 'predictor' not in st.session_state:
    st.session_state.predictor = MolecularPropertyPredictor()
    # Load and train the model
    data = pd.read_csv(Path(__file__).parent.parent / 'data/molecules.csv')
    st.session_state.predictor.fit(
        data['SMILES'].values,
        data['Polarity'].values,
        data['Solubility'].values
    )

# Title and description
st.title("🧬 AI Molecular Property Visualizer")
st.markdown("""
This application uses AI to predict and visualize molecular properties. Upload or select a molecule to see its 3D structure and predicted properties.
""")

# Sidebar
with st.sidebar:
    st.header("Molecule Input")
    input_type = st.radio(
        "Choose input method:",
        ["Select from database", "Enter SMILES"]
    )
    
    if input_type == "Select from database":
        data = pd.read_csv(Path(__file__).parent.parent / 'data/molecules.csv')
        molecule_name = st.selectbox(
            "Select a molecule:",
            data['Molecule'].tolist()
        )
        smiles = data[data['Molecule'] == molecule_name]['SMILES'].iloc[0]
        st.write(f"SMILES: {smiles}")
    else:
        smiles = st.text_input(
            "Enter SMILES string:",
            value="CCO"  # Ethanol as default
        )

# Main content
col1, col2 = st.columns([2, 1])

with col1:
    st.header("3D Visualization")
    
    try:
        # Generate 3D structure
        mol_viewer = MoleculeViewer()
        mol_3d = mol_viewer.generate_3d_coords(smiles)
        
        if mol_3d:
            # Create viewer HTML
            viewer_html = mol_viewer.get_viewer_html(
                mol_3d,
                size=(800, 500)
            )
            
            if viewer_html:
                st.components.v1.html(
                    viewer_html,
                    height=550,
                    width=800
                )
            else:
                st.error("Failed to generate 3D viewer.")
        else:
            st.error("Could not generate 3D structure for the given SMILES string.")
    except Exception as e:
        st.error(f"Error in 3D visualization: {str(e)}")

with col2:
    st.header("Property Predictions")
    
    if smiles:
        # Get predictions
        polarity, solubility = st.session_state.predictor.predict(smiles)
        confidence_pol, confidence_sol = st.session_state.predictor.predict_proba(smiles)
        
        # Display predictions with improved styling
        st.markdown("### Predicted Properties")
        
        # Polarity prediction
        st.markdown(
            f"""
            <div style='padding: 1.5em; border-radius: 8px; background-color: #1E1E1E; margin-bottom: 1em;'>
                <h4 style='color: #FFFFFF; margin: 0;'>Polarity</h4>
                <p style='font-size: 1.4em; color: #FFFFFF; margin: 0.5em 0;'>{polarity}</p>
                <p style='color: #B0B0B0;'>Confidence: {confidence_pol:.1%}</p>
            </div>
            """,
            unsafe_allow_html=True
        )
        
        # Solubility prediction
        st.markdown(
            f"""
            <div style='padding: 1.5em; border-radius: 8px; background-color: #1E1E1E; margin-bottom: 1em;'>
                <h4 style='color: #FFFFFF; margin: 0;'>Solubility</h4>
                <p style='font-size: 1.4em; color: #FFFFFF; margin: 0.5em 0;'>{solubility}</p>
                <p style='color: #B0B0B0;'>Confidence: {confidence_sol:.1%}</p>
            </div>
            """,
            unsafe_allow_html=True
        )
    
    # Visualization Options
    st.markdown("### Visualization Options")
    view_style = st.radio(
        "Molecule style:",
        ["stick", "sphere", "cartoon"],
        key="viz_style"
    )
    
    if mol_3d and view_style != "stick":
        with col1:
            viewer_html = mol_viewer.get_viewer_html(
                mol_3d,
                size=(800, 500),
                style=view_style
            )
            st.components.v1.html(viewer_html, height=550, scrolling=False)
