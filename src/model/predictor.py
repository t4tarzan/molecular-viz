from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
import numpy as np

class MolecularPropertyPredictor:
    def __init__(self):
        self.polarity_model = RandomForestClassifier(n_estimators=100, random_state=42)
        self.solubility_model = RandomForestClassifier(n_estimators=100, random_state=42)
        self.polarity_encoder = LabelEncoder()
        self.solubility_encoder = LabelEncoder()
        # Initialize Morgan fingerprint generator
        self.morgan_gen = rdMolDescriptors.GetMorganGenerator(2, 1024)
    
    def _get_fingerprint(self, smiles):
        """Generate Morgan fingerprint for a molecule."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Use Morgan Generator instead of direct fingerprint calculation
        return list(self.morgan_gen.GetFingerprint(mol))
    
    def fit(self, smiles_list, polarity_labels, solubility_labels):
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
        
        # Convert labels
        polarity_labels = self.polarity_encoder.fit_transform(polarity_labels[valid_indices])
        solubility_labels = self.solubility_encoder.fit_transform(solubility_labels[valid_indices])
        
        # Train models
        self.polarity_model.fit(fingerprints, polarity_labels)
        self.solubility_model.fit(fingerprints, solubility_labels)
    
    def predict(self, smiles):
        """Predict properties for a single SMILES string."""
        fingerprint = self._get_fingerprint(smiles)
        if fingerprint is None:
            return None, None
        
        # Make predictions
        polarity_pred = self.polarity_encoder.inverse_transform(
            self.polarity_model.predict([fingerprint])
        )[0]
        solubility_pred = self.solubility_encoder.inverse_transform(
            self.solubility_model.predict([fingerprint])
        )[0]
        
        return polarity_pred, solubility_pred
    
    def predict_proba(self, smiles):
        """Get prediction probabilities for a single SMILES string."""
        fingerprint = self._get_fingerprint(smiles)
        if fingerprint is None:
            return 0.0, 0.0
        
        # Get probabilities
        polarity_proba = np.max(self.polarity_model.predict_proba([fingerprint])) * 100
        solubility_proba = np.max(self.solubility_model.predict_proba([fingerprint])) * 100
        
        return polarity_proba, solubility_proba
