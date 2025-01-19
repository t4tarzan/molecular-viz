from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import py3Dmol

# Test with methane
smiles = 'C'  # Methane
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
print(f"Number of atoms before 3D: {mol.GetNumAtoms()}")

# Generate 3D coordinates
result = AllChem.EmbedMolecule(mol, randomSeed=42)
print(f"Embedding result: {result}")
print(f"Number of conformers: {mol.GetNumConformers()}")

if mol.GetNumConformers() > 0:
    print("\nAtom coordinates:")
    conf = mol.GetConformer()
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        print(f"{atom.GetSymbol()}: ({pos.x:.3f}, {pos.y:.3f}, {pos.z:.3f})")

    # Save a 2D depiction
    img = Draw.MolToImage(mol)
    img.save("molecule_2d.png")
    print("\nSaved 2D structure as 'molecule_2d.png'")

    # Print PDB format
    print("\nPDB format:")
    print(Chem.MolToPDBBlock(mol))
