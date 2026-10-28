"""
Davis/Klaeger compound overlap, by PubChem CID -> synonym match against
Klaeger's drug-name index (revision for Molecular Informatics, Reviewer 1,
suggested platform-vs-definition benchmark).

Davis drugs are identified only by PubChem CID in this repository
(davis_drugs.csv); Klaeger drugs are identified by common/code name
(klaeger_matrix.csv index). This script resolves each Davis CID to its
PubChem synonym list via the PUG REST API and matches (case/punctuation
-insensitive) against the 222 Klaeger drug names, to find the compounds
profiled in both datasets.
"""
import script_logging; script_logging.capture(__file__)
import os
import re
import time
import json
import urllib.request
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def norm(s):
    return re.sub(r'[^A-Z0-9]', '', s.upper())


def pubchem_synonyms(cid, retries=3):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=15) as r:
                data = json.load(r)
            return data['InformationList']['Information'][0]['Synonym']
        except Exception as e:
            if attempt == retries - 1:
                print(f"  CID {cid}: FAILED ({e})")
                return []
            time.sleep(1.0)
    return []


davis = pd.read_csv(os.path.join(HERE, "davis_drugs.csv"))
klaeger_names = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0).index.tolist()
klaeger_norm = {norm(n): n for n in klaeger_names}

cache_path = os.path.join(HERE, "davis_pubchem_synonyms_cache.json")
cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}

matches = []
for _, row in davis.iterrows():
    cid = str(int(row['CID']))
    if cid in cache:
        syns = cache[cid]
    else:
        syns = pubchem_synonyms(row['CID'])
        cache[cid] = syns
        time.sleep(0.25)
    hit = None
    for syn in syns:
        n = norm(syn)
        if n in klaeger_norm:
            hit = klaeger_norm[n]
            break
    if hit:
        matches.append((row['Drug_Index'], row['CID'], hit))

json.dump(cache, open(cache_path, 'w'))

print(f"Resolved {len(davis)} Davis CIDs against {len(klaeger_names)} Klaeger drug names.")
print(f"Matches found: {len(matches)}")
for drug_idx, cid, klaeger_name in matches:
    print(f"  Davis Drug_Index={drug_idx} (CID {cid})  <->  Klaeger '{klaeger_name}'")

pd.DataFrame(matches, columns=['Davis_Drug_Index', 'Davis_CID', 'Klaeger_Name']).to_csv(
    os.path.join(HERE, "davis_klaeger_overlap.csv"), index=False)
print("\nSaved davis_klaeger_overlap.csv")
