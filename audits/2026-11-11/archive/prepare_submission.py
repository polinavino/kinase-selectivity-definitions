"""Build reviewable submission assets without changing manuscript sources.

Output is a DRAFT package: author metadata, catalogue keywords, the revision
redline and scientific audit findings still require resolution.
"""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from zipfile import ZipFile, ZIP_DEFLATED
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "submission_prepared"


def finalize_document(path):
    """Keep editable content, mark correspondence and collect figure legends."""
    with ZipFile(path) as archive:
        files = {name: archive.read(name) for name in archive.namelist()}
    document = files["word/document.xml"].decode()
    document = document.replace('>Polina Vinogradova</w:t>', '>Polina Vinogradova*</w:t>', 1)
    for acronym, expansion in [
        ("FDA", "U.S. Food and Drug Administration (FDA)"),
        ("ATP", "adenosine triphosphate (ATP)"),
        ("CATDS", "concentration- and target-dependent selectivity (CATDS)"),
        ("CID", "compound identifier (CID)"),
        ("FAERS", "FDA Adverse Event Reporting System (FAERS)"),
    ]:
        if expansion not in document:
            document, count = re.subn(rf'\b{acronym}\b', expansion, document, count=1)
            assert count == 1, f"Could not locate first use of {acronym}"
    paragraphs = re.findall(r'<w:p(?:\s[^>]*)?>.*?</w:p>', document, flags=re.S)
    captions = [p for p in paragraphs if 'w:val="ImageCaption"' in p]
    if captions:
        for caption in captions:
            document = document.replace(caption, "", 1)
        heading = ('<w:p><w:pPr><w:pStyle w:val="Heading1"/></w:pPr>'
                   '<w:r><w:t>Figure legends</w:t></w:r></w:p>')
        document = document.replace(heading, "")
        position = document.rfind("<w:sectPr")
        if position < 0:
            position = document.rfind("</w:body>")
        document = document[:position] + heading + "".join(captions) + document[position:]
    files["word/document.xml"] = document.encode()
    with ZipFile(path, "w", ZIP_DEFLATED) as archive:
        for name, data in files.items():
            archive.writestr(name, data)


def number_table(path, number):
    with ZipFile(path) as archive:
        files = {name: archive.read(name) for name in archive.namelist()}
    document = files["word/document.xml"].decode()
    for paragraph in re.findall(r'<w:p(?:\s[^>]*)?>.*?</w:p>', document, flags=re.S):
        if 'w:val="TableCaption"' in paragraph and f"Table {number}." not in paragraph:
            numbered = re.sub(r'(<w:t\b[^>]*>)', rf'\g<1>Table {number}. ', paragraph, count=1)
            document = document.replace(paragraph, numbered, 1)
            break
    files["word/document.xml"] = document.encode()
    with ZipFile(path, "w", ZIP_DEFLATED) as archive:
        for name, data in files.items():
            archive.writestr(name, data)


def validate_assets():
    from PIL import Image
    figures = {}
    for number in ("Figure1", "Figure2", "Figure3"):
        with Image.open(OUT / (number + ".tif")) as im:
            dpi = [float(value) for value in im.info["dpi"]]
            figures[number] = dict(pixels=list(im.size), dpi=dpi,
                                   width_cm=im.width / dpi[0] * 2.54)
            assert dpi[0] >= 600 and im.width > 1800
    manuscript = OUT / "molecular_informatics_submission_CURRENT_DRAFT.docx"
    with ZipFile(manuscript) as archive:
        tree = ET.fromstring(archive.read("word/document.xml"))
        w = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
             "m": "http://schemas.openxmlformats.org/officeDocument/2006/math"}
        text = " ".join(node.text or "" for node in tree.findall(".//w:t", w))
        checks = dict(current_title="Analysis Choices and Ranking Stability" in text,
                      obsolete_D2_absent="bounded gap sensitivity" not in text,
                      corresponding_author_marked="Polina Vinogradova*" in text,
                      legends_collected_at_end=text.find("Figure legends") > text.find("References"),
                      tables=len(tree.findall(".//w:tbl", w)),
                      figures=len(tree.findall(".//w:drawing", w)),
                      editable_equations=len(tree.findall(".//m:oMath", w)))
    assert checks["current_title"] and checks["obsolete_D2_absent"]
    assert checks["corresponding_author_marked"] and checks["legends_collected_at_end"]
    assert checks["tables"] == 4 and checks["figures"] == 3
    for number in range(1, 5):
        with ZipFile(OUT / f"Table{number}.docx") as archive:
            tree = ET.fromstring(archive.read("word/document.xml"))
            caption = "".join(node.text or "" for node in tree.findall(".//w:t", w))
            assert f"Table {number}." in caption
    (OUT / "asset_checks.json").write_text(json.dumps(dict(document=checks, figures=figures), indent=2))


def main():
    OUT.mkdir(exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix="kinase-submission-assets-"))
    for source in ROOT.iterdir():
        if source.is_file() and source.suffix in {".py", ".csv", ".xls", ".xlsx", ".json"}:
            shutil.copy2(source, stage / source.name)
    shutil.copytree(ROOT / "paper", stage / "paper")
    env = {**os.environ, "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(stage / "mpl"),
           "KMP_DUPLICATE_LIB_OK": "TRUE"}
    jobs = [("metric_fidelity_robustness.py", "Figure1", "metric_fidelity"),
            ("binding_profiles.py", "Figure2", "binding_profiles"),
            ("panel_recovery.py", "Figure3", "panel_recovery")]
    # Re-render from data and plotting code, rather than upsampling old PNGs.
    wrapper = r'''
import runpy, sys
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image
script, outdir, number = sys.argv[1:]
original = plt.savefig
def savefig(*args, **kwargs):
    original(*args, **kwargs)
    fig = plt.gcf()
    out = Path(outdir)
    fig.savefig(out / (number + '.pdf'), bbox_inches='tight')
    fig.savefig(out / (number + '.png'), dpi=600, bbox_inches='tight')
    with Image.open(out / (number + '.png')) as im:
        # Encode the genuine rendered pixels at a 17.5 cm final width.
        dpi = im.width / (17.5 / 2.54)
        im.save(out / (number + '.tif'), compression='tiff_lzw', dpi=(dpi, dpi))
plt.savefig = savefig
runpy.run_path(script, run_name='__main__')
'''

    def render(job):
        script, number, _ = job
        with (OUT / (number + "_render.log")).open("w") as log:
            subprocess.run([sys.executable, "-c", wrapper, str(stage / script), str(OUT), number],
                           cwd=stage, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)

    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(render, jobs))
    for _, number, basename in jobs:
        shutil.copy2(OUT / (number + ".png"), stage / "paper/figures" / (basename + ".png"))
    subprocess.run(["bash", "build_docx.sh"], cwd=stage / "paper", env=env, check=True)
    manuscript = OUT / "molecular_informatics_submission_CURRENT_DRAFT.docx"
    shutil.copy2(stage / "paper/molecular_informatics_submission.docx", manuscript)
    finalize_document(manuscript)
    for number in range(1, 5):
        table_path = OUT / f"Table{number}.docx"
        subprocess.run(["pandoc", str(stage / f"paper/submission_tables/table{number}.tex"),
                        "-o", str(table_path)], check=True)
        number_table(table_path, number)
    toc = (ROOT / "paper/graphical_abstract_text.txt").read_text().strip().split("\n\n")
    toc_source = stage / "toc.md"
    toc_source.write_text("\n\n".join(toc) + "\n\n![](toc_graphic.png){width=5cm}\n")
    shutil.copy2(ROOT / "paper/figures/toc_graphic.png", stage / "toc_graphic.png")
    subprocess.run(["pandoc", str(toc_source), "--resource-path", str(stage),
                    "-o", str(OUT / "Graphical_abstract.docx")], check=True)
    for name in ["supporting_information.pdf", "cover_letter.txt",
                 "response_to_reviewers_round2.txt", "FIGURE_LEGENDS.txt", "abstract_plain.txt",
                 "graphical_abstract_text.txt"]:
        shutil.copy2(ROOT / "paper" / name, OUT / name)
    letter = OUT / "cover_letter.txt"
    fit = ("This study fits Molecular Informatics because it evaluates how computational "
           "summaries of kinase profiling data support compound comparisons in drug "
           "discovery, and provides empirical and formal guidance for interpreting "
           "those summaries.\n\n")
    letter.write_text(letter.read_text().replace("I am grateful", fit + "I am grateful", 1))
    validate_assets()
    print("Prepared draft assets:", OUT)


if __name__ == "__main__":
    main()
