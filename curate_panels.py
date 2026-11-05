"""
Curated target panels: one reference target per human protein kinase gene.

The released panels count assay targets, not kinases. Davis profiles mutant,
phosphorylation-state and single-domain constructs of the same gene (12 EGFR
constructs, BRAF and BRAF(V600E)) as well as lipid kinases and three non-human
kinases; Klaeger's Kinobeads matrix includes cyclins, other complex subunits and
non-kinase proteins, and four columns are protein groups (e.g. PRKCE;PRKCH);
Anastassiadis profiles some kinases in more than one complex or isoform
(CDK1/cyclin A and CDK1/cyclin B). Selectivity across kinases needs one entry per
kinase, chosen without reference to any compound's affinities:

  Davis          human protein kinases only (lipid and non-human kinases removed);
                 per gene, the wild-type construct with no domain or state label,
                 otherwise the wild-type catalytic (JH1) or N-terminal kinase
                 domain; genes profiled only as mutants are removed.
  Klaeger        protein kinases (klaeger_protein_kinases.py), without the four
                 protein-group columns.
  Anastassiadis  per gene, the first-listed complex or isoform.
  Metz           already one human protein kinase per column; unchanged.

Reads the released files (*_released.csv) and writes the curated panels under
the names every analysis script reads (davis_affinity.csv, klaeger_matrix.csv,
anastassiadis_matrix.csv).
"""
import script_logging; script_logging.capture(__file__)
import re
import pandas as pd

# ---------------------------------------------------------------- Davis
prot = pd.read_csv("davis_proteins.csv").set_index("Protein_Index")["Gene_Name"]
gene = prot.str.replace(r"\(.*$", "", regex=True).str.replace(r"-(phosphorylated|nonphosphorylated)$", "", regex=True)
non_human = prot[prot.str.contains(r"P\.falciparum|M\.tuberculosis")].index
lipid = prot[gene.str.match(r"PIK3|PIK4|PIP4K|PIP5K")].index
label = prot.str.extract(r"\((.*)\)")[0].fillna("")
mutant = label.str.contains(r"[A-Z]\d+[A-Z]|del|ins") & ~label.str.contains(r"Kin\.Dom|JH\d")
mutant |= label.str.contains(r"S808G")            # GCN2(Kin.Dom.2,S808G)

keep = []
for g, idx in gene.groupby(gene).groups.items():
    idx = [i for i in idx if i not in non_human and i not in lipid and not mutant[i]]
    if not idx:
        continue
    bare = [i for i in idx if label[i] == ""]
    pref = bare or [i for i in idx if re.search(r"JH1|Kin\.Dom\.1", label[i])] or idx
    keep.append(sorted(pref)[0])
keep = sorted(keep)
aff = pd.read_csv("davis_affinity_released.csv")
aff[aff["Protein_Index"].isin(keep)].to_csv("davis_affinity.csv", index=False)
dropped = sorted(set(prot.index) - set(keep))
print(f"Davis: {len(prot)} constructs -> {len(keep)} reference targets "
      f"({len(non_human)} non-human, {len(lipid)} lipid-kinase constructs removed; "
      f"{len(dropped) - len(non_human) - len(lipid)} further mutant, state or domain constructs not used)")
print("  references for genes with several constructs: "
      + ", ".join(prot[i] for i in keep if (gene == gene[i]).sum() > 1))
print("  genes with no wild-type construct: "
      + ", ".join(sorted(set(gene[mutant]) - set(gene[keep]) - set(gene[lipid]))))

# ---------------------------------------------------------------- Klaeger
src = open("klaeger_protein_kinases.py").read()
non_pk = re.findall(r'"([A-Z0-9;]+)"', src[src.index("NON_PK = {"):src.index("}\n", src.index("NON_PK = {"))])
km = pd.read_csv("klaeger_matrix_released.csv", index_col=0)
pk = [c for c in km.columns if c not in non_pk]
groups = [c for c in pk if ";" in c]
cur = [c for c in pk if ";" not in c]
km[cur].to_csv("klaeger_matrix.csv")
print(f"\nKlaeger: {km.shape[1]} targets -> {len(pk)} protein kinases -> {len(cur)} after removing "
      f"{len(groups)} protein-group columns ({', '.join(groups)})")

# ---------------------------------------------------------------- Anastassiadis
an = pd.read_csv("anastassiadis_matrix_released.csv", index_col=0)
def an_gene(name):
    n = name.split("/")[0].split(" cyclin")[0].strip()
    return {"LYN B": "LYN", "PKCb1": "PKCb", "PKCb2": "PKCb", "PKG1a": "PKG1", "PKG1b": "PKG1"}.get(n, n)
first = {}
for c in an.columns:
    first.setdefault(an_gene(c), c)
an_keep = [c for c in an.columns if first[an_gene(c)] == c]
an[an_keep].to_csv("anastassiadis_matrix.csv")
print(f"\nAnastassiadis: {an.shape[1]} targets -> {len(an_keep)} "
      f"(duplicates removed: {', '.join(c for c in an.columns if c not in an_keep)})")

print("\nMetz: unchanged (172 human protein kinases, one per column)")
