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

# Add custom CSS to adjust layout
st.markdown("""
<style>
    .main > div {
        padding-top: 1rem;
    }
    .stApp > header {
        height: 0px;
    }
    div[data-testid="stSidebarNav"] {
        display: none;
    }
    .css-1544g2n.e1fqkh3o4 {
        padding-top: 0rem;
    }
    .block-container {
        padding-top: 1rem;
        max-width: 100%;
        padding-right: 0;
        padding-left: 0;
    }
    .css-1y4p8pa {
        padding-top: 0rem;
    }
    .css-1r6slb0.e1tzin5v2 {
        gap: 0rem;
    }
    /* Remove padding from columns */
    div.css-12w0qpk.e1tzin5v2 {
        padding: 0 !important;
        gap: 0 !important;
    }
    div.css-1r6slb0.e1tzin5v2 {
        padding: 0 !important;
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
        'Solubility': ['Insoluble', 'Soluble', 'Insoluble', 'Soluble'],
        'Reactivity': ['Low', 'Medium', 'High', 'Low'],
        'Stability': ['High', 'Medium', 'Low', 'High']
    })

# Initialize predictor with fallback data
if 'predictor' not in st.session_state:
    st.session_state.predictor = MolecularPropertyPredictor()
    # Always fit with the current data, whether from file or fallback
    st.session_state.predictor.fit(
        data['SMILES'].values,
        {
            'polarity': data['Polarity'].values,
            'solubility': data['Solubility'].values,
            'reactivity': data['Reactivity'].values,
            'stability': data['Stability'].values
        }
    )

# Title and description
st.markdown("# 🧬 AI Molecular Property Visualizer")
st.markdown("""
This application uses AI to predict and visualize molecular properties. Upload or select a molecule to see its 3D structure and predicted properties.
""")

# Function to convert SMILES to molecular formula
def get_molecular_formula(smiles):
    """Convert SMILES to molecular formula for display."""
    formulas = {
        'C': 'CH4',
        'CC': 'CH3CH3',
        'C=C': 'CH2=CH2',
        'CCO': 'CH3CH2OH',
        'c1ccccc1': 'C6H6',
        'CC(=O)O': 'CH3COOH',
        'C#C': 'HC≡CH',
        'O': 'H2O',
        'N': 'NH3',
        'CCC': 'CH3CH2CH3',
        'ClC(Cl)Cl': 'CHCl3',
        'c1ccc(cc1)O': 'C6H5OH',
        'CC=O': 'CH3CHO'
    }
    return formulas.get(smiles, smiles)

# Create three columns with exact proportions
left_col, middle_col, right_col = st.columns([20, 60, 20])

with left_col:
    # Molecule Input Section
    st.subheader("Molecule Input")
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
    st.subheader("Molecular Formula:")
    st.code(get_molecular_formula(smiles), language=None)
    st.markdown("---")
    st.subheader("SMILES:")
    st.code(smiles)

with middle_col:
    # 3D Visualization Section
    st.header("3D Visualization")
    
    # Visualization options
    viz_style = st.radio(
        "Molecule style:",
        ["stick", "sphere", "structure"],
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

with right_col:
    # Properties Section
    st.header("Properties")
    
    # Make predictions
    predictions = st.session_state.predictor.predict(smiles)
    confidence_scores = st.session_state.predictor.predict_proba(smiles)
    
    # Property cards with styling
    properties_order = ['polarity', 'solubility', 'reactivity', 'stability']
    property_names = {
        'polarity': 'POLARITY',
        'solubility': 'SOLUBILITY',
        'reactivity': 'REACTIVITY',
        'stability': 'STABILITY'
    }
    
    for prop in properties_order:
        st.markdown(
            f"""
            <div class="property-card">
                <p class="property-title">{property_names[prop]}</p>
                <p class="property-value">{predictions[prop]}</p>
                <div class="confidence-bar">
                    <div class="confidence-fill" style="width: {confidence_scores[prop]}%;"></div>
                </div>
                <p class="confidence-text">Confidence: {confidence_scores[prop]:.1f}%</p>
            </div>
            """,
            unsafe_allow_html=True
        )
