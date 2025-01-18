# Molecular Property Visualizer

An AI-powered molecular visualization tool that predicts and visualizes molecular properties, specifically polarity and solubility. Built with Streamlit, RDKit, and py3Dmol.

## Features

- **3D Molecular Visualization**: Interactive 3D visualization of molecular structures
- **Property Prediction**: AI-powered prediction of molecular properties
  - Polarity (Polar/Non-polar)
  - Solubility (Soluble/Insoluble)
- **Multiple Input Methods**:
  - Select from pre-defined molecules
  - Input custom SMILES strings
- **Visualization Options**:
  - Stick model
  - Space-filling (sphere) model
  - Cartoon representation

## Live Demo

Visit the live application at: [Your Streamlit App URL]

## Local Development

1. Clone the repository:
   ```bash
   git clone https://github.com/t4tarzan/molecular-viz.git
   cd molecular-viz
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Run the application:
   ```bash
   streamlit run app.py
   ```

## Tech Stack

- **Frontend**: Streamlit
- **Molecular Processing**: RDKit
- **3D Visualization**: py3Dmol
- **Machine Learning**: scikit-learn
- **Data Processing**: pandas

## Contributing

Feel free to open issues or submit pull requests for any improvements.

## License

MIT License
