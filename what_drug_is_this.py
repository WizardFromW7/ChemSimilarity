"""
What Drug Is This Most Like? 🧪
--------------------------------
Give it a SMILES string and it tells you which approved drugs
your molecule most resembles, using Tanimoto similarity on
Morgan fingerprints.

Requirements:
    pip install rdkit

Usage:
    python what_drug_is_this.py "CC(=O)Oc1ccccc1C(=O)O"
    (that's aspirin, if you want a giveaway example)


What the Hell is Morgan Fingerprints and Tanimoto Similarity, Anyway?
--------------------------------
Morgan fingerprints are a way of representing molecules as a series of
bits (0s and 1s) that encode the presence or absence of certain
substructures. They are widely used in cheminformatics for tasks like
similarity searching and machine learning.

The Tanimoto similarity is a metric that compares two fingerprints and
gives a score between 0 and 1, where 1 means the molecules are identical
in terms of the features captured by the fingerprints.


How to Make Your Own SMILES?
---------------------------------
Quick rules to make one up:

    - Capital letters = atoms: C carbon, N nitrogen, O oxygen, F/Cl/Br halogens
    - Lowercase = aromatic ring atoms (c, n, o)
    - =  double bond, # triple bond
    - () = branch off the main chain
    - Numbers = ring closure (two atoms sharing the same number get
      bonded to close a ring)
    - No letter between atoms = single bond, assumed

Before we kick off, make sure you have the RDKit library installed.
You can do this via pip:
    pip install rdkit
"""

import sys
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors

# ---------------------------------------------------------
# A small, fun reference library of well-known drugs.
# Feel free to add your own favourites / weird ones.
# ---------------------------------------------------------
DRUG_LIBRARY = {
    "Paracetamol":  ("CC(=O)Nc1ccc(O)cc1", "The one everyone has in their cupboard."),
    "Aspirin":      ("CC(=O)Oc1ccccc1C(=O)O", "The OG painkiller, been around since 1897."),
    "Ibuprofen":    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "NSAID, your gym-injury best friend."),
    "Caffeine":     ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "The reason you're functional this morning."),
    "Nicotine":     ("CN1CCC[C@H]1c1cccnc1", "Highly addictive alkaloid, obviously."),
    "Metformin":    ("CN(C)C(=N)NC(=N)N", "First-line type 2 diabetes drug."),
    "Sildenafil":   ("CCCc1nn(C)c2c1nc(nc2=O)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1", "Yes, that one 👀"),
    "Diazepam":     ("CN1c2ccc(Cl)cc2C(=NCC1=O)c1ccccc1", "Valium — classic benzodiazepine."),
    "Morphine":     ("CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "Opioid, the reference point for all others."),
    "Amoxicillin":  ("CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O", "Common penicillin-type antibiotic."),
    "Cetirizine":   ("OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1", "Antihistamine, hayfever's nemesis."),
    "Omeprazole":   ("COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C", "PPI, stomach acid's off switch."),
    "Salbutamol":   ("CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", "Blue inhaler, asthma rescue drug."),
    "Warfarin":     ("CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "Blood thinner, also rat poison lol."),
    "Fluoxetine":   ("CNCCC(c1ccc(cc1)C(F)(F)F)Oc1ccccc1", "Prozac — SSRI antidepressant."),
    "MDMA":         ("CC(NC)Cc1ccc2c(c1)OCO2", "Not medically prescribed, but chemically interesting."),
    "THC":          ("CCCCCc1cc(O)c2c(c1)OC(C)(C)[C@@H]1CCC(C)=CC12", "The active compound in cannabis."),
    "Sertraline":   ("CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21", "Another SSRI, brand name Zoloft."),
    "Levothyroxine":("Ic1cc(CC(N)C(=O)O)cc(I)c1Oc1cc(I)c(O)c(I)c1", "Thyroid hormone replacement."),
    "Codeine":      ("COc1ccc2c3c1O[C@H]1[C@@H](O)C=C[C@H](C3)[C@H](N(C)CC1)2", "Mild opioid, cough syrup staple."),
}


def fingerprint(smiles: str):
    """Turn a SMILES string into a Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return mol, fp


def find_closest_drugs(query_smiles: str, top_n: int = 5):
    query_mol, query_fp = fingerprint(query_smiles)
    if query_mol is None:
        print("❌ That doesn't look like a valid SMILES string. Try again!")
        return

    mw = Descriptors.MolWt(query_mol)
    logp = Descriptors.MolLogP(query_mol)

    results = []
    for name, (smi, fun_fact) in DRUG_LIBRARY.items():
        mol, fp = fingerprint(smi)
        if mol is None:
            continue
        similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
        results.append((name, similarity, fun_fact))

    results.sort(key=lambda x: x[1], reverse=True)

    print(f"\n🔬 Your molecule — MW: {mw:.1f}, LogP: {logp:.2f}\n")
    print(f"Top {top_n} closest known drugs:\n")
    for i, (name, sim, fact) in enumerate(results[:top_n], start=1):
        bar = "█" * int(sim * 20)
        print(f"{i}. {name:<15} {sim*100:5.1f}%  {bar}")
        print(f"   ↳ {fact}\n")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        smiles_input = sys.argv[1]
    else:
        smiles_input = input("Paste a SMILES string: ").strip()

    find_closest_drugs(smiles_input)
