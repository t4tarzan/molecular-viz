# AI-Powered Molecular Property Visualization

A Streamlit-based web application that uses AI to predict and visualize molecular properties.

## Features
- 3D molecular visualization using RDKit and 3Dmol.js
- AI-powered prediction of molecular properties (polarity, solubility, reactivity, and stability)
- Interactive molecule selection and SMILES input
- Real-time property predictions with confidence scores
- Multiple visualization styles: stick, sphere, and structure views

## Installation
```bash
git clone https://github.com/t4tarzan/molecular-viz.git
cd molecular-viz
pip install -r requirements.txt
```

## Usage
```bash
streamlit run app.py
```

## Dependencies
- streamlit>=1.28.0
- rdkit>=2023.3.1
- scikit-learn>=1.0.2
- pandas>=1.3.0
- py3Dmol>=2.0.0
- matplotlib>=3.4.3

## Project Structure
```
molecular-viz/
├── app.py                 # Main Streamlit application
├── data/
│   └── molecules.csv      # Dataset of molecules and properties
├── src/
│   ├── model/
│   │   └── predictor.py   # ML model for property prediction
│   └── visualization/
│       └── mol_viewer.py  # 3D molecule visualization
├── projectcode.md         # Detailed code documentation
├── FAQ.md                 # Frequently Asked Questions
└── requirements.txt       # Project dependencies
```

## Recent Updates

### Molecular Visualization Improvements (January 19, 2025)

#### Key Changes
1. **Enhanced Hydrogen Visualization**
   - Added `keepH: true` option when loading molecules in 3Dmol.js
   - Implemented separate styling for hydrogen atoms with improved visibility
   - Fixed issues with hydrogen display in methane and other simple molecules

2. **Atom-Specific Styling**
   - Carbon atoms: Gray (sphere radius: 0.5-0.8, stick radius: 0.2)
   - Hydrogen atoms: White (sphere radius: 0.3-0.4, stick radius: 0.1)
   - Oxygen atoms: Red (sphere radius: 0.5-0.8, stick radius: 0.2)
   - Nitrogen atoms: Blue (sphere radius: 0.5-0.8, stick radius: 0.2)

3. **Visualization Modes**
   - Stick mode: Shows all atoms with bonds
   - Sphere mode: Emphasizes atomic positions
   - Structure mode: Uses Jmol color scheme for standard representation

4. **Testing and Validation**
   - Created `test_methane.py` for isolated testing of molecular visualization
   - Verified correct display of simple molecules (methane, water)
   - Confirmed proper visualization of complex molecules (glucose)

#### Technical Details
- Using 3Dmol.js for molecular visualization
- RDKit for molecular processing and 3D coordinate generation
- SMILES notation for molecular representation
- PDB format for 3D structure visualization

#### Known Issues Resolved
- Fixed missing hydrogen atoms in visualization
- Improved bond visibility between atoms
- Enhanced contrast with black background
- Corrected atomic radii for better proportions

## Rollback Versions
Latest stable version: commit `fcaa34c` (2025-01-19)
Features:
- Enhanced visualization styles:
  - Stick view (bonds with small atoms)
  - Sphere view (space-filling)
  - Structure view (combined stick and sphere)
- Fixed 3D coordinate generation
- Improved property predictions:
  - Polarity
  - Solubility
  - Reactivity
  - Stability
- Proper sidebar width and layout
- Working molecule visualization
- Property cards with confidence bars

Previous stable version: commit `363911f` (2025-01-19)
Features:
- Basic visualization modes
- Initial property predictions
- Sidebar layout

## Contributing
Feel free to open issues or submit pull requests. Please ensure that any new features maintain compatibility with the existing visualization and prediction functionality.

## License
MIT License - see LICENSE file for details
