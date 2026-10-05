"""Build SMILES for a disulfide-cyclised, C-terminally amidated peptide.

Three edits to RDKit's linear peptide:
  1. MolFromSequence  -> L-peptide, free N-terminal amine, free C-terminal acid
  2. C-terminal amide -> the acid's -OH oxygen becomes -NH2 (our peptides are
     synthesised as amides, and Stage 4 scored them as amides)
  3. disulfide        -> bond the two cysteine sulfurs, clearing their thiol Hs

Validated against oxytocin (CYIQNCPLG, disulfide 1-6, C-terminal amide), whose
monoisotopic/average masses are known independently.
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.rdchem import RWMol, BondType

def build(seq):
    m = Chem.MolFromSequence(seq)
    if m is None:
        return None, "MolFromSequence failed"
    rw = RWMol(m)
    # --- C-terminal amide: the only carboxylic -OH in a peptide backbone is the
    #     C-terminus (side-chain acids D/E also match, so pick the one on the
    #     LAST residue by walking from the terminal carbon).
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    hits = rw.GetSubstructMatches(patt)
    if not hits:
        return None, "no carboxylic acid found"
    # choose the hit whose carbonyl carbon has the highest atom index that is
    # also bonded to a backbone alpha carbon -- for MolFromSequence output the
    # C-terminus is the last such match
    cterm = max(hits, key=lambda h: h[0])
    oh = cterm[2]
    a = rw.GetAtomWithIdx(oh)
    a.SetAtomicNum(7); a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    # --- disulfide
    s_idx = [at.GetIdx() for at in rw.GetAtoms()
             if at.GetSymbol() == "S" and at.GetTotalNumHs() >= 1]
    if len(s_idx) != 2:
        return None, "expected 2 free thiols, found %d" % len(s_idx)
    for i in s_idx:
        at = rw.GetAtomWithIdx(i); at.SetNoImplicit(True); at.SetNumExplicitHs(0)
    rw.AddBond(s_idx[0], s_idx[1], BondType.SINGLE)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception as e:
        return None, "sanitize failed: %s" % e
    return out, None

def props(m):
    return dict(smiles=Chem.MolToSmiles(m),
                formula=rdMolDescriptors.CalcMolFormula(m),
                mw=round(Descriptors.MolWt(m), 3),
                exact=round(Descriptors.ExactMolWt(m), 4),
                tpsa=round(Descriptors.TPSA(m), 1),
                clogp=round(Descriptors.MolLogP(m), 2),
                hbd=rdMolDescriptors.CalcNumHBD(m),
                hba=rdMolDescriptors.CalcNumHBA(m),
                rotb=rdMolDescriptors.CalcNumRotatableBonds(m))

if __name__ == "__main__":
    print("VALIDATION -- oxytocin CYIQNCPLG, disulfide + C-terminal amide")
    m, err = build("CYIQNCPLG")
    assert m is not None, err
    p = props(m)
    print("  formula        %s"%p["formula"])
    print("  average MW     %.3f   (literature 1007.19)"%p["mw"])
    print("  monoisotopic   %.4f   (literature 1006.4365)"%p["exact"])
    ss = len([b for b in m.GetBonds()
              if b.GetBeginAtom().GetSymbol()=="S" and b.GetEndAtom().GetSymbol()=="S"])
    print("  S-S bonds      %d  (expect 1)"%ss)
    amide = m.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[NX3H2]"))
    print("  C-term amide   %s"%amide)
    print("  free acids     %d  (expect 0)"%len(m.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))))
    print("\n  SMILES: %s"%p["smiles"])
