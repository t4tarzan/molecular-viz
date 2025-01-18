from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
import numpy as np
import warnings

class MolecularPropertyPredictor:
    def __init__(self):
        # Initialize models for each property
        self.models = {
            'polarity': RandomForestClassifier(n_estimators=100, random_state=42),
            'solubility': RandomForestClassifier(n_estimators=100, random_state=42),
            'reactivity': RandomForestClassifier(n_estimators=100, random_state=42),
            'stability': RandomForestClassifier(n_estimators=100, random_state=42)
        }
        # Initialize label encoders
        self.encoders = {
            'polarity': LabelEncoder(),
            'solubility': LabelEncoder(),
            'reactivity': LabelEncoder(),
            'stability': LabelEncoder()
        }
    
    def _get_fingerprint(self, smiles):
        """Generate Morgan fingerprint for a molecule."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Suppress deprecation warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
    
    def fit(self, smiles_list, labels_dict):
        """Train the models on the provided data."""
        # Generate fingerprints
        fingerprints = []
        valid_indices = []
        
        for i, smiles in enumerate(smiles_list):
            fp = self._get_fingerprint(smiles)
            if fp is not None:
                fingerprints.append(fp)
                valid_indices.append(i)
        
        if not fingerprints:
            raise ValueError("No valid molecules found in the training data")
        
        # Train models for each property
        for prop in self.models.keys():
            if prop in labels_dict:
                # Convert labels
                labels = self.encoders[prop].fit_transform(
                    [labels_dict[prop][i] for i in valid_indices]
                )
                # Train model
                self.models[prop].fit(fingerprints, labels)
    
    def predict(self, smiles):
        """Predict all properties for a single SMILES string."""
        fingerprint = self._get_fingerprint(smiles)
        if fingerprint is None:
            return {prop: None for prop in self.models.keys()}
        
        predictions = {}
        for prop in self.models.keys():
            pred = self.encoders[prop].inverse_transform(
                self.models[prop].predict([fingerprint])
            )[0]
            predictions[prop] = pred
        
        return predictions
    
    def predict_proba(self, smiles):
        """Get prediction probabilities for a single SMILES string."""
        fingerprint = self._get_fingerprint(smiles)
        if fingerprint is None:
            return {prop: 0.0 for prop in self.models.keys()}
        
        probabilities = {}
        for prop in self.models.keys():
            proba = np.max(self.models[prop].predict_proba([fingerprint])) * 100
            probabilities[prop] = proba
        
        return probabilities
