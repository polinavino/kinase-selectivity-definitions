import requests
import pandas as pd
import time

fda_approved = {
    'Abemaciclib':'abemaciclib','Afatinib':'afatinib','Alectinib':'alectinib',
    'Axitinib':'axitinib','Baricitinib':'baricitinib','Binimetinib':'binimetinib',
    'Bosutinib':'bosutinib','Cabozantinib':'cabozantinib','Capmatinib':'capmatinib',
    'Ceritinib':'ceritinib','Cobimetinib':'cobimetinib','Crizotinib':'crizotinib',
    'Dabrafenib':'dabrafenib','Dacomitinib':'dacomitinib','Dasatinib':'dasatinib',
    'Encorafenib':'encorafenib','Entrectinib':'entrectinib','Erlotinib':'erlotinib',
    'Fedratinib':'fedratinib','Fostamatinib':'fostamatinib','Gefitinib':'gefitinib',
    'Gilteritinib':'gilteritinib','Ibrutinib':'ibrutinib','Imatinib':'imatinib',
    'Lapatinib':'lapatinib','Lenvatinib':'lenvatinib','Midostaurin':'midostaurin',
    'Momelotinib':'momelotinib','Neratinib':'neratinib','Nilotinib':'nilotinib',
    'Nintedanib':'nintedanib','Osimertinib':'osimertinib','Palbociclib':'palbociclib',
    'Pazopanib':'pazopanib','Pexidartinib':'pexidartinib','Ponatinib':'ponatinib',
    'Quizartinib':'quizartinib','Regorafenib':'regorafenib','Ribociclib':'ribociclib',
    'Ruxolitinib':'ruxolitinib','Selumetinib':'selumetinib','Sorafenib':'sorafenib',
    'Sunitinib':'sunitinib','Tofacitinib':'tofacitinib','Trametinib':'trametinib',
    'Vandetanib':'vandetanib','Vemurafenib':'vemurafenib',
}

base_url = "https://api.fda.gov/drug/event.json"
results = []

for klaeger_name, generic_name in fda_approved.items():
    row = {'drug': klaeger_name}
    base_search = f'patient.drug.openfda.generic_name:"{generic_name}"'
    queries = {
        'total_reports':           base_search,
        'serious_reports':         f'{base_search} AND serious:"1"',
        'death_reports':           f'{base_search} AND seriousnessdeath:"1"',
        'hospitalization_reports': f'{base_search} AND seriousnesshospitalization:"1"',
        'lifethreat_reports':      f'{base_search} AND seriousnesslifethreatening:"1"',
    }
    for label, search in queries.items():
        try:
            r = requests.get(f"{base_url}?search={search}&limit=1", timeout=15)
            row[label] = r.json().get('meta',{}).get('results',{}).get('total',0)
        except:
            row[label] = None
        time.sleep(0.25)
    results.append(row)
    print(f"{klaeger_name}: total={row['total_reports']}, serious={row['serious_reports']}")

df = pd.DataFrame(results)
df.to_csv("faers_counts.csv", index=False)
print("Saved faers_counts.csv")
