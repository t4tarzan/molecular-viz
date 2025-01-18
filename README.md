# AI-Powered Molecular Property Visualization

A Streamlit-based web application that uses AI to predict and visualize molecular properties.

## Features
- 3D molecular visualization using RDKit and 3Dmol.js
- AI-powered prediction of molecular properties (polarity and solubility)
- Interactive molecule selection and SMILES input
- Real-time property predictions with confidence scores

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
└── requirements.txt       # Project dependencies
```

## Rollback Version
Latest working version: commit `363911f` (2025-01-19)
Features:
- Proper sidebar width and layout
- Working 3D visualization for stick and sphere modes
- Improved property predictions UI:
  - Compact cards with progress bars
  - Visual confidence indicators
  - Optimized spacing and typography
- Clean code organization

To rollback to this version:
```bash
git reset --hard 363911f
git clean -fd
```

## Development
The application is built using:
- Streamlit for the web interface
- RDKit for molecular processing
- scikit-learn for ML predictions
- 3Dmol.js for 3D visualization

## Deployment
The app is deployed on Streamlit Cloud and automatically updates when changes are pushed to the main branch.

## Contributing
1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request
