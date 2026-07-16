"""
What Drug Is This Most Like? 🧪 — Web App Edition
----------------------------------------------------
Same logic as the terminal version, just wrapped in Streamlit
so it's an actual clickable web app instead of a script.

Requirements:
    pip install rdkit streamlit

Run it with:
    streamlit run app.py

(NOT "python app.py" — Streamlit apps are launched differently!)
"""

import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Draw

# ---------------------------------------------------------
# Same drug library as before
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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return mol, fp


def find_closest_drugs(query_smiles: str, top_n: int = 5):
    query_mol, query_fp = fingerprint(query_smiles)
    if query_mol is None:
        return None, None, None

    mw = Descriptors.MolWt(query_mol)
    logp = Descriptors.MolLogP(query_mol)

    results = []
    for name, (smi, fun_fact) in DRUG_LIBRARY.items():
        mol, fp = fingerprint(smi)
        if mol is None:
            continue
        similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
        results.append((name, similarity, fun_fact, smi))

    results.sort(key=lambda x: x[1], reverse=True)
    return query_mol, (mw, logp), results[:top_n]


# ===========================================================
# STREAMLIT UI STARTS HERE
# ===========================================================

st.set_page_config(page_title="What Drug Is This?", page_icon="🧪")

st.title("🧪 What Drug Is This Most Like?")
st.write(
    "Paste a SMILES string and find out which approved drugs your molecule "
    "most resembles, based on structural similarity (Tanimoto / Morgan fingerprints)."
)

# A few clickable examples so people don't need to know SMILES to try it
st.write("**Try an example:**")
example_cols = st.columns(4)
examples = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Random weird one": "CCNOBrF",
}

if "smiles_input" not in st.session_state:
    st.session_state.smiles_input = ""

for col, (label, smi) in zip(example_cols, examples.items()):
    if col.button(label):
        st.session_state.smiles_input = smi

smiles_input = st.text_input(
    "Or paste your own SMILES string:",
    value=st.session_state.smiles_input,
    placeholder="e.g. CC(=O)Nc1ccc(O)cc1",
)

top_n = st.slider("How many matches to show?", min_value=3, max_value=10, value=5)

if smiles_input:
    query_mol, props, results = find_closest_drugs(smiles_input, top_n=top_n)

    if query_mol is None:
        st.error("❌ That doesn't look like a valid SMILES string. Check for typos!")
    else:
        mw, logp = props

        col1, col2 = st.columns([1, 2])
        with col1:
            st.image(Draw.MolToImage(query_mol, size=(250, 250)))
        with col2:
            st.metric("Molecular Weight", f"{mw:.1f} g/mol")
            st.metric("LogP", f"{logp:.2f}")

        st.subheader(f"Top {top_n} closest known drugs")

        top_score = results[0][1]

        for i, (name, sim, fact, smi) in enumerate(results, start=1):
            st.write(f"**{i}. {name}** — {sim*100:.1f}%")
            st.progress(sim)
            st.caption(fact)

        st.divider()
        if top_score > 0.8:
            st.warning("🚨 Basically the same drug, legally distinct.")
        elif top_score > 0.5:
            st.info("👀 Suspiciously similar.")
        else:
            st.success("🧬 A genuinely novel little guy.")
else:
    st.caption("👆 Click an example above or paste your own SMILES to get started.")
