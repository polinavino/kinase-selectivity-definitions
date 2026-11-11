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
  Klaeger        the targets Klaeger et al. annotate as protein/lipid kinases
                 ("Target Annotation" sheet of aan4368_Table_S2.xlsx), without lipid
                 kinases, protein-group columns, and Q6ZSR9, which UniProt describes
                 as an uncharacterized protein with no kinase domain, encoded in the
                 3'-UTR of AAK1. This removes BCR and TRRAP, which the source does not
                 annotate as kinases and flags as possibly co-enriched with ABL1 and
                 CDK8. Measured values below the floor (Kd app above 10 uM) are set to
                 the floor, pKd = 5.0, the value given to pairs with no detected binding.
  Anastassiadis  per gene, the first-listed complex or isoform.
  Metz           already one human protein kinase per column; unchanged.

Reads the released files (*_released.csv) and writes the curated panels under
the names every analysis script reads (davis_affinity.csv, klaeger_matrix.csv,
anastassiadis_matrix.csv), and one table per dataset listing every released
target, whether it is kept, and why (target_curation_<dataset>.csv).
"""
import script_logging; script_logging.capture(__file__)
import re
import pandas as pd
from panel_curation import anastassiadis_gene as an_gene, curate_anastassiadis

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
def davis_reason(i):
    if i in keep:
        return "kept", "reference construct for " + gene[i]
    if i in non_human:
        return "removed", "non-human kinase"
    if i in lipid:
        return "removed", "lipid kinase"
    if mutant[i]:
        return "removed", "mutant construct"
    ref = [k for k in keep if gene[k] == gene[i]]
    return "removed", "another construct of " + gene[i] + (f" is the reference ({prot[ref[0]]})" if ref else "")
pd.DataFrame([(i, prot[i], gene[i]) + davis_reason(i) for i in prot.index],
             columns=["Protein_Index", "construct", "gene", "status", "reason"]).to_csv("target_curation_davis.csv", index=False)
print(f"Davis: {len(prot)} constructs -> {len(keep)} reference targets "
      f"({len(non_human)} non-human, {len(lipid)} lipid-kinase constructs removed; "
      f"{len(dropped) - len(non_human) - len(lipid)} further mutant, state or domain constructs not used)")
print("  references for genes with several constructs: "
      + ", ".join(prot[i] for i in keep if (gene == gene[i]).sum() > 1))
print("  genes with no wild-type construct: "
      + ", ".join(sorted(set(gene[mutant]) - set(gene[keep]) - set(gene[lipid]))))

# ---------------------------------------------------------------- Klaeger
FLOOR_K = 5.0
ann = pd.read_excel("aan4368_Table_S2.xlsx", sheet_name="Target Annotation").set_index("Protein")
km = pd.read_csv("klaeger_matrix_released.csv", index_col=0)
LIPID = {"PIP4K2A", "PIP4K2C", "PIP4K2A;PIP4K2B"}
NO_KINASE_DOMAIN = {"Q6ZSR9": "uncharacterized protein with no kinase domain, encoded in the 3'-UTR of AAK1 (UniProt)"}
def klaeger_reason(c):
    a = ann.loc[c]
    co = pd.notna(a["might be coenriched"])
    partner = f" (possibly co-enriched with {a['interacting kinases (BioGrid)']})" if co else ""
    if pd.isna(a["protein/lipid kinase"]):
        return "removed", "not annotated as a kinase" + partner
    if c in LIPID:
        return "removed", "lipid kinase"
    if ";" in c:
        return "removed", "protein group"
    if c in NO_KINASE_DOMAIN:
        return "removed", NO_KINASE_DOMAIN[c]
    return "kept", "annotated as a protein kinase" + partner
rows = [(c,) + klaeger_reason(c) for c in km.columns]
pd.DataFrame(rows, columns=["target", "status", "reason"]).to_csv("target_curation_klaeger.csv", index=False)
cur = [c for c, st, _ in rows if st == "kept"]
below = int((km[cur].values < FLOOR_K).sum())
km[cur].clip(lower=FLOOR_K).to_csv("klaeger_matrix.csv")
print(f"\nKlaeger: {km.shape[1]} targets -> {len(cur)} protein kinases annotated by the source")
n_ann = int(ann.loc[km.columns, "protein/lipid kinase"].notna().sum())
print(f"  annotated by the source as protein or lipid kinases: {n_ann}; not annotated as kinases: {km.shape[1] - n_ann}")
for reason in sorted({r for _, st, r in rows if st == "removed" and not r.startswith("not annotated as a kinase (")}):
    names = [c for c, st, r in rows if r == reason]
    print(f"  removed, {reason}: {len(names)}" + (f" ({', '.join(names)})" if len(names) <= 10 else ""))
co_removed = [(c, r) for c, st, r in rows if r.startswith("not annotated as a kinase (")]
print(f"  removed, not annotated as a kinase and flagged as possibly co-enriched: {len(co_removed)}; "
      f"of these, {', '.join(c for c, _ in co_removed if c in ('BCR', 'TRRAP'))} were in the previous curated panel")
print(f"  measured values below the floor set to {FLOOR_K}: {below} entries")

# ---------------------------------------------------------------- Anastassiadis
an = pd.read_csv("anastassiadis_matrix_released.csv", index_col=0)
curated_an = curate_anastassiadis(an)
an_keep = list(curated_an.columns)
first = {an_gene(c): c for c in an_keep}
curated_an.to_csv("anastassiadis_matrix.csv")
pd.DataFrame([(c, an_gene(c), "kept" if c in an_keep else "removed",
               "first-listed complex or isoform of " + an_gene(c) if c in an_keep
               else "another complex or isoform of " + an_gene(c) + f" is the reference ({first[an_gene(c)]})")
              for c in an.columns], columns=["target", "gene", "status", "reason"]).to_csv(
    "target_curation_anastassiadis.csv", index=False)
print(f"\nAnastassiadis: {an.shape[1]} targets -> {len(an_keep)} "
      f"(duplicates removed: {', '.join(c for c in an.columns if c not in an_keep)})")

print("\nMetz: unchanged (172 human protein kinases, one per column)")
