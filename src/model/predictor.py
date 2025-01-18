from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
import numpy as np

class MolecularPropertyPredictor:
    def __init__(self):
        self.polarity_model = RandomForestClassifier(random_state=42)
        self.solubility_model = RandomForestClassifier(random_state=42)
        self.polarity_encoder = LabelEncoder()
        self.solubility_encoder = LabelEncoder()
        
    def _calculate_fingerprints(self, smiles):
        """Calculate Morgan fingerprints for a molecule."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048))
    
    def fit(self, X_smiles, y_polarity, y_solubility):
        """Train the models on the provided data."""
        # Prepare features
        X = [self._calculate_fingerprints(smiles) for smiles in X_smiles]
        X = [x for x in X if x is not None]
        
        # Encode labels
        y_polarity_encoded = self.polarity_encoder.fit_transform(y_polarity)
        y_solubility_encoded = self.solubility_encoder.fit_transform(y_solubility)
        
        # Train models
        self.polarity_model.fit(X, y_polarity_encoded)
        self.solubility_model.fit(X, y_solubility_encoded)
    
    def predict(self, smiles):
        """Predict properties for a given SMILES string."""
        fingerprint = self._calculate_fingerprints(smiles)
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
        """Predict probabilities for properties."""
        fingerprint = self._calculate_fingerprints(smiles)
        if fingerprint is None:
            return None, None
        
        # Get prediction probabilities
        polarity_proba = np.max(self.polarity_model.predict_proba([fingerprint]))
        solubility_proba = np.max(self.solubility_model.predict_proba([fingerprint]))
        
        return polarity_proba, solubility_proba
