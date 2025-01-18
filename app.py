import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
from pathlib import Path
import sys

# Page config must be the first Streamlit command
st.set_page_config(
    page_title="AI Molecular Property Visualizer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS to adjust column widths and property cards
st.markdown("""
<style>
    [data-testid="column"] {
        width: calc(20% - 1rem) !important;
        flex: 1 1 calc(20% - 1rem) !important;
        min-width: 200px !important;
    }
    [data-testid="column"]:first-child {
        width: 80% !important;
        flex: 1 1 80% !important;
    }
    .property-card {
        background-color: #1E1E1E;
        padding: 15px;
        border-radius: 8px;
        margin: 8px 0;
    }
    .property-title {
        color: #B0B0B0;
        font-size: 0.9em;
        margin: 0;
        font-weight: 600;
    }
    .property-value {
        color: #FFFFFF;
        font-size: 1.2em;
        margin: 4px 0;
        font-weight: 500;
    }
    .confidence-bar {
        width: 100%;
        height: 4px;
        background-color: #2D2D2D;
        border-radius: 2px;
        margin: 8px 0;
    }
    .confidence-fill {
        height: 100%;
        background-color: #4CAF50;
        border-radius: 2px;
        transition: width 0.3s ease;
    }
    .confidence-text {
        color: #B0B0B0;
        font-size: 0.8em;
        margin: 2px 0;
    }
</style>
""", unsafe_allow_html=True)

# Add src directory to path
root_dir = Path(__file__).parent
sys.path.append(str(root_dir))

from src.model.predictor import MolecularPropertyPredictor
from src.visualization.mol_viewer import MoleculeViewer
import py3Dmol

# Load data
data_path = root_dir / 'data' / 'molecules.csv'
try:
    data = pd.read_csv(data_path)
except FileNotFoundError:
    # Fallback data for both local and deployment
    print(f"Could not find {data_path}, using fallback data")
    data = pd.DataFrame({
        'Molecule': ['Methane', 'Ethanol', 'Benzene', 'Acetone'],
        'SMILES': ['C', 'CCO', 'C1=CC=CC=C1', 'CC(=O)C'],
        'Polarity': ['Non-polar', 'Polar', 'Non-polar', 'Polar'],
        'Solubility': ['Insoluble', 'Soluble', 'Insoluble', 'Soluble']
    })

# Initialize predictor with fallback data
if 'predictor' not in st.session_state:
    st.session_state.predictor = MolecularPropertyPredictor()
    # Always fit with the current data, whether from file or fallback
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

# Sidebar for molecule input
with st.sidebar:
    st.title("Molecule Input")
    
    # Input method selection
    input_type = st.radio(
        "Choose input method:",
        ["Select from database", "Enter SMILES"],
        key="input_type"
    )

    if input_type == "Select from database":
        molecule_name = st.selectbox(
            "Select a molecule:",
            data['Molecule'].tolist()
        )
        smiles = data[data['Molecule'] == molecule_name]['SMILES'].iloc[0]
    else:
        smiles = st.text_input(
            "Enter SMILES string:",
            value="C1=CC=CC=C1",  # Benzene as default
            help="Enter a valid SMILES string for your molecule"
        )

    st.markdown("---")
    st.markdown("**SMILES:**")
    st.code(smiles)

# Main content area with two columns
col1, col2 = st.columns([3, 1])

with col1:
    st.header("3D Visualization")
    
    # Visualization options
    viz_style = st.radio(
        "Molecule style:",
        ["stick", "sphere", "cartoon"],
        horizontal=True,
        key="viz_style"
    )
    
    try:
        # Generate 3D structure
        mol_viewer = MoleculeViewer()
        mol_3d = mol_viewer.generate_3d_coords(smiles)
        
        if mol_3d:
            # Create viewer HTML
            viewer_html = mol_viewer.get_viewer_html(
                mol_3d,
                size=(800, 500),
                style=viz_style
            )
            
            if viewer_html:
                components.html(
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
    st.header("Properties")
    
    # Make predictions
    polarity, solubility = st.session_state.predictor.predict(smiles)
    confidence_pol, confidence_sol = st.session_state.predictor.predict_proba(smiles)
    
    # Property cards with styling
    st.markdown(
        f"""
        <div class="property-card">
            <p class="property-title">POLARITY</p>
            <p class="property-value">{polarity}</p>
            <div class="confidence-bar">
                <div class="confidence-fill" style="width: {confidence_pol}%;"></div>
            </div>
            <p class="confidence-text">Confidence: {confidence_pol:.1f}%</p>
        </div>
        
        <div class="property-card">
            <p class="property-title">SOLUBILITY</p>
            <p class="property-value">{solubility}</p>
            <div class="confidence-bar">
                <div class="confidence-fill" style="width: {confidence_sol}%;"></div>
            </div>
            <p class="confidence-text">Confidence: {confidence_sol:.1f}%</p>
        </div>
        """,
        unsafe_allow_html=True
    )
