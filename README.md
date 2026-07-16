# 🧪 What Drug Is This Most Like?

A small tool that takes a SMILES string (a text representation of a molecule)
and tells you which well-known drugs it structurally resembles most, using
Tanimoto similarity on Morgan fingerprints (RDKit).

Built as a fun summer side project to play around with cheminformatics before
starting an MSc — not a rigorous screening tool, just a way to explore how
structural similarity search works.

## Example

```
Paste a SMILES string: CC(=O)Nc1ccc(O)cc1

🔬 Your molecule — MW: 151.2, LogP: 1.35

Top 5 closest known drugs:

1. Paracetamol     100.0%  ████████████████████
   ↳ The one everyone has in their cupboard.
2. Aspirin          22.2%  ████
   ↳ The OG painkiller, been around since 1897.
3. Amoxicillin      21.8%  ████
   ↳ Common penicillin-type antibiotic.
4. Warfarin         18.8%  ███
   ↳ Blood thinner, also rat poison lol.
5. Ibuprofen        18.4%  ███
   ↳ NSAID, your gym-injury best friend.
```

## How it works

1. Converts each molecule's SMILES string into a **Morgan fingerprint**
   (a bit-vector representation of its structural features).
2. Compares the input molecule's fingerprint against a small library of
   ~20 well-known drugs using **Tanimoto similarity**.
3. Ranks and returns the closest matches.

## Setup

```bash
pip install rdkit
```

## Usage

**Command line:**
```bash
python what_drug_is_this.py "CC(=O)Nc1ccc(O)cc1"
```

Or run it without an argument and paste a SMILES string when prompted:
```bash
python what_drug_is_this.py
```

**Web app (Streamlit):**
```bash
pip install streamlit
streamlit run app.py
```

## Try it with these SMILES strings

| Drug | SMILES |
|---|---|
| Aspirin | `CC(=O)Oc1ccccc1C(=O)O` |
| Caffeine | `Cn1cnc2c1c(=O)n(C)c(=O)n2C` |
| Paracetamol | `CC(=O)Nc1ccc(O)cc1` |

Or find SMILES for any other drug on [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
— search a drug name and copy the "Canonical SMILES" field.

## Notes

- Similarity scores are based purely on 2D structural fingerprints — they
  don't account for 3D shape, binding affinity, or actual pharmacological
  effect. A high score means "structurally similar," not "does the same thing."
- The drug reference library is small and hand-picked for fun/recognisability,
  not comprehensive.

## Built with

- [RDKit](https://www.rdkit.org/) — cheminformatics toolkit
- [Streamlit](https://streamlit.io/) — web app framework (optional UI layer)
