"""Check source/Word/figure/table consistency and inventory current artifacts."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

from PIL import Image
from docx_tools import NS, W, M, read_docx, validate_main

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tokens(element):
    return " ".join(n.text or "" for n in element.iter()
                    if n.tag in {f"{{{W}}}t", f"{{{M}}}t"})


def main():
    checks = json.loads((HERE / "submission_checks.json").read_text())
    for name, expected in checks["source_hashes"].items():
        assert digest(HERE / name) == expected, f"Build is stale: {name}"
    files = read_docx(HERE / "molecular_informatics_submission.docx")
    document = ET.fromstring(files["word/document.xml"])
    validate_main(HERE / "molecular_informatics_submission.docx")
    tables = document.findall('.//w:tbl', NS)
    for number, table in enumerate(tables, 1):
        standalone = ET.fromstring(read_docx(HERE / f"submission_tables/Table{number}.docx")["word/document.xml"])
        assert tokens(table) == tokens(standalone.find('.//w:tbl', NS)), f"Stale Table {number}"
    relationships = ET.fromstring(files['word/_rels/document.xml.rels'])
    rel = {r.get('Id'): r.get('Target') for r in relationships}
    drawing_ns = {**NS, 'pic': 'http://schemas.openxmlformats.org/drawingml/2006/picture',
                  'a': 'http://schemas.openxmlformats.org/drawingml/2006/main'}
    for drawing in document.findall('.//w:drawing', NS):
        picture = drawing.find('.//pic:cNvPr', drawing_ns)
        source = HERE / picture.get('descr')
        embedded_id = drawing.find('.//a:blip', drawing_ns).get(
            '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed')
        assert hashlib.sha256(files['word/' + rel[embedded_id]]).hexdigest() == digest(source), source
    mapping = {'metric_fidelity': 'Figure1', 'binding_profiles': 'Figure2', 'panel_recovery': 'Figure3'}
    figure_checks = {}
    for basename, number in mapping.items():
        for suffix in ('.png', '.pdf'):
            a, b, c = ROOT / (basename + suffix), HERE / 'figures' / (basename + suffix), HERE / 'submission_figures' / (number + suffix)
            assert digest(a) == digest(b) == digest(c), f"Stale {number}{suffix}"
        with Image.open(HERE / f'submission_figures/{number}.tif') as im:
            dpi = [float(v) for v in im.info['dpi']]
            assert min(dpi) >= 600
            with Image.open(HERE / f'figures/{basename}.png') as png:
                assert im.size == png.size and im.tobytes() == png.tobytes()
            figure_checks[number] = dict(pixels=list(im.size), dpi=dpi)
    for basename in ('cross_dataset_correlations', 'failure_modes', 'instability_by_family', 'panel_size_stability'):
        for suffix in ('.png', '.pdf'):
            assert digest(ROOT / (basename + suffix)) == digest(HERE / 'figures' / (basename + suffix))
    source = (HERE / 'main.tex').read_text()
    for name in re.findall(r'\\input\{([^}]+)\}', source):
        source += '\n' + (HERE / (name + '.tex')).read_text()
    labels = re.findall(r'\\label\{([^}]+)\}', source)
    assert len(labels) == len(set(labels)), 'Duplicate manuscript labels'
    assert set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', source)) <= set(labels)
    bibliography = (HERE / 'references.bib').read_text()
    keys = set(re.findall(r'@\w+\{([^,]+),', bibliography))
    for text in (source, (HERE / 'supporting_information.tex').read_text()):
        for group in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', text):
            assert set(group.split(',')) <= keys, group
    for name in ('main', 'supporting_information'):
        log = (HERE / (name + '.log')).read_text(errors='replace')
        assert not re.search(r'undefined|Overfull|^! ', log, flags=re.M), f'Check {name}.log'
    inventory = []
    for directory, dirs, names in os.walk(ROOT):
        dirs[:] = sorted(d for d in dirs if d not in {'.git', '__pycache__', '.vscode', '.claude', '.pytest_cache'})
        for name in sorted(names):
            path = Path(directory) / name
            relative = path.relative_to(ROOT)
            if name == '.DS_Store' or name == 'submission_manifest.json' or path.suffix in {'.aux','.log','.blg','.bbl','.fls','.fdb_latexmk','.out','.synctex.gz'}:
                continue
            kind = 'current'
            if 'archive' in relative.parts or name in {'AUDIT_FINDINGS.md','findings.txt','contribution_review.txt','verification.json'}:
                kind = 'historical_provenance'
            elif name in {'faers_counts.csv','selectivity_outcomes_merged.csv','faers_pull.txt','clinical_safety_data.csv','label_texts.json','extended_labels.json'}:
                kind = 'frozen_or_retired_data_provenance'
            elif name == 'reviewer1_round2_review.txt':
                kind = 'review_input'
            inventory.append(dict(path=str(relative), role=kind, bytes=path.stat().st_size, sha256=digest(path)))
    output = dict(document_checks=validate_main(HERE/'molecular_informatics_submission.docx'),
                  figure_checks=figure_checks, files=inventory)
    (HERE / 'submission_manifest.json').write_text(json.dumps(output, indent=2) + '\n')
    print(f"Consistency checks passed; inventoried {len(inventory)} files. Historical files are explicitly classified.")


if __name__ == '__main__':
    main()
