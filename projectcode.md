# Building a Molecular Visualization and Property Prediction App

This document provides a detailed explanation of how the molecular visualization and property prediction app was built, combining physical models, visual charts, and AI/ML capabilities.

## Project Overview

This project combines three key components:
1. Physical clay models of molecular structures
2. Visual reference chart of molecular representations
3. Interactive AI-powered web application

## Technical Components

### 1. Data Management (`molecules.csv`)
```python
# Sample data structure
data = {
    'Molecule': ['Methane', 'Ethanol', 'Benzene', ...],
    'SMILES': ['C', 'CCO', 'C1=CC=CC=C1', ...],
    'Polarity': ['Non-Polar', 'Polar', 'Non-Polar', ...],
    'Solubility': ['Insoluble', 'Soluble', 'Insoluble', ...]
}
```

The data file contains:
- Molecule names
- SMILES notation (chemical structure representation)
- Polarity classification
- Solubility classification

### 2. Machine Learning Model (`predictor.py`)

```python
class MolecularPropertyPredictor:
    def __init__(self):
        # Initialize RandomForest classifiers
        self.polarity_model = RandomForestClassifier(n_estimators=100)
        self.solubility_model = RandomForestClassifier(n_estimators=100)
```

Key features:
- Uses RandomForest for classification
- Generates Morgan fingerprints for molecular features
- Predicts both polarity and solubility
- Provides confidence scores

Example usage:
```python
predictor = MolecularPropertyPredictor()
predictor.fit(smiles_list, polarity_labels, solubility_labels)
polarity, solubility = predictor.predict("CCO")  # Ethanol
confidence_pol, confidence_sol = predictor.predict_proba("CCO")
```

### 3. 3D Visualization (`mol_viewer.py`)

```python
class MoleculeViewer:
    def generate_3d_coords(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        return mol
```

Visualization styles:
1. **Stick Mode**:
```javascript
viewer.setStyle({}, {
    "stick": {
        "radius": 0.2,
        "color": "lightgray"
    },
    "sphere": {
        "radius": 0.3
    }
});
```

2. **Sphere Mode**:
```javascript
viewer.setStyle({}, {
    "sphere": {
        "radius": 1.0,
        "colorscheme": "Jmol"
    }
});
```

3. **Cartoon Mode**:
```javascript
viewer.setStyle({}, {
    "cartoon": {
        "color": "spectrum",
        "thickness": 0.8,
        "opacity": 0.8
    }
});
```

### 4. Web Interface (`app.py`)

The Streamlit interface is organized into three sections:

1. **Input Section (Left Sidebar)**:
```python
with st.sidebar:
    input_type = st.radio(
        "Choose input method:",
        ["Select from database", "Enter SMILES"]
    )
```

2. **Visualization Section (Main Area)**:
```python
viewer_html = mol_viewer.get_viewer_html(
    mol_3d,
    size=(800, 500),
    style=viz_style
)
```

3. **Properties Section (Right Column)**:
```python
st.markdown(
    f"""
    <div class="property-card">
        <p class="property-title">POLARITY</p>
        <p class="property-value">{polarity}</p>
        <div class="confidence-bar">
            <div class="confidence-fill" style="width: {confidence_pol}%;"></div>
        </div>
    </div>
    """
)
```

## CSS Styling

The app uses custom CSS for styling:
```css
.property-card {
    background-color: #1E1E1E;
    padding: 15px;
    border-radius: 8px;
    margin: 8px 0;
}

.confidence-bar {
    width: 100%;
    height: 4px;
    background-color: #2D2D2D;
    border-radius: 2px;
}
```

## Integration with Physical Models

The app complements physical clay models by:
1. Providing accurate 3D representations
2. Allowing rotation and manipulation
3. Showing different visualization styles
4. Adding property predictions

## Dependencies

```
streamlit>=1.28.0
rdkit>=2023.3.1
scikit-learn>=1.0.2
pandas>=1.3.0
py3Dmol>=2.0.0
matplotlib>=3.4.3
```

## Running the Application

```bash
streamlit run app.py
```

## Code Organization

```
molecular-viz/
├── app.py                 # Main application
├── data/
│   └── molecules.csv      # Molecule database
├── src/
│   ├── model/
│   │   └── predictor.py   # ML model
│   └── visualization/
│       └── mol_viewer.py  # 3D visualization
└── requirements.txt       # Dependencies
```

## Future Enhancements

Potential improvements:
1. Add more molecules to the database
2. Implement more property predictions
3. Add export functionality for 3D models
4. Include animation of molecular interactions
5. Add educational tooltips and explanations
