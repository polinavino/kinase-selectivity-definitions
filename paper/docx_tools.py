"""OpenXML formatting and validation for the generated submission documents."""
from pathlib import Path
import re
from zipfile import ZipFile, ZIP_DEFLATED
import xml.etree.ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
NS = {"w": W, "m": M}


def read_docx(path):
    with ZipFile(path) as archive:
        return {name: archive.read(name) for name in archive.namelist()}


def write_docx(path, files):
    temporary = Path(str(path) + ".tmp")
    with ZipFile(temporary, "w", ZIP_DEFLATED) as archive:
        for name, data in files.items():
            archive.writestr(name, data)
    temporary.replace(path)


def prepare_main(path):
    files = read_docx(path)
    document = files["word/document.xml"].decode()
    captions = [p for p in re.findall(r'<w:p(?:\s[^>]*)?>.*?</w:p>', document, flags=re.S)
                if 'w:val="ImageCaption"' in p]
    assert len(captions) == 3
    for caption in captions:
        document = document.replace(caption, "", 1)
    heading = ('<w:p><w:pPr><w:pStyle w:val="Heading1"/></w:pPr>'
               '<w:r><w:t>Figure legends</w:t></w:r></w:p>')
    document = document.replace(heading, "")
    position = document.rfind("<w:sectPr")
    if position < 0:
        position = document.rfind("</w:body>")
    document = document[:position] + heading + "".join(captions) + document[position:]
    document = re.sub(r'<w:pgMar\b[^>]*/>',
                      '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" '
                      'w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>', document)
    # Pandoc's default reference has an implicit page size. Make geometry
    # explicit so Word, Pages and the journal converter use the same canvas.
    if '<w:pgSz' not in document:
        document = document.replace('<w:sectPr>', '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>', 1)
    if '<w:pgMar' not in document:
        document = document.replace('<w:sectPr>', '<w:sectPr><w:pgMar w:top="1440" w:right="1440" '
                                    'w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>', 1)
    content_types = files['[Content_Types].xml'].decode()
    if 'Extension="png"' not in content_types:
        content_types = content_types.replace('</Types>', '<Default Extension="png" ContentType="image/png"/></Types>')
    files['[Content_Types].xml'] = content_types.encode()
    maximum = 5943600  # 6.5 inches in EMUs
    for drawing in re.findall(r'<w:drawing>.*?</w:drawing>', document, flags=re.S):
        extent = re.search(r'<wp:extent cx="(\d+)" cy="(\d+)"', drawing)
        if extent and int(extent[1]) > maximum:
            cx, cy = int(extent[1]), int(extent[2])
            height = round(cy * maximum / cx)
            resized = re.sub(r'(<(?:wp:extent|a:ext) cx=")' + str(cx) +
                             r'(" cy=")' + str(cy) + '"',
                             rf'\g<1>{maximum}\g<2>{height}"', drawing)
            document = document.replace(drawing, resized, 1)
    files["word/document.xml"] = document.encode()
    write_docx(path, files)


def number_table(path, number):
    files = read_docx(path)
    document = files["word/document.xml"].decode()
    for paragraph in re.findall(r'<w:p(?:\s[^>]*)?>.*?</w:p>', document, flags=re.S):
        if 'w:val="TableCaption"' in paragraph:
            if re.search(r'Table \d+\.', paragraph):
                numbered = re.sub(r'Table \d+\.', f'Table {number}.', paragraph, count=1)
            else:
                numbered = re.sub(r'(<w:t\b[^>]*>)', rf'\g<1>Table {number}. ', paragraph, count=1)
            document = document.replace(paragraph, numbered, 1)
            break
    files["word/document.xml"] = document.encode()
    write_docx(path, files)


def validate_main(path):
    files = read_docx(path)
    document = ET.fromstring(files["word/document.xml"])
    paragraphs = ["".join(n.text or "" for n in p.iter()
                          if n.tag in {f"{{{W}}}t", f"{{{M}}}t"})
                  for p in document.findall(".//w:p", NS)]
    text = "\n".join(paragraphs)
    checks = dict(current_title="Analysis Choices and Ranking Stability" in text,
                  corresponding_author="Polina Vinogradova*" in text,
                  tables=len(document.findall(".//w:tbl", NS)),
                  figures=len(document.findall(".//w:drawing", NS)),
                  editable_equations=len(document.findall(".//m:oMath", NS)),
                  bibliography_entries=len(re.findall(r'^\[\d+\]', text, flags=re.M)),
                  figure_legends=[p for p in paragraphs if re.match(r'^Figure [123]\.', p)],
                  unresolved_labels=re.findall(r'\b(?:tab|fig|sec|eq):[\w:-]+|\?\?|CITATION NOT FOUND|Error! Reference', text))
    targets = {b.get(f"{{{W}}}name") for b in document.findall(".//w:bookmarkStart", NS)}
    checks["missing_internal_targets"] = [link.get(f"{{{W}}}anchor")
        for link in document.findall(".//w:hyperlink", NS)
        if link.get(f"{{{W}}}anchor") and link.get(f"{{{W}}}anchor") not in targets]
    assert checks["current_title"] and checks["corresponding_author"]
    assert checks["tables"] == 4 and checks["figures"] == 3
    assert checks["editable_equations"] > 100 and checks["bibliography_entries"] > 15
    assert len(checks["figure_legends"]) == 3 and not checks["unresolved_labels"]
    assert not checks["missing_internal_targets"]
    return checks
