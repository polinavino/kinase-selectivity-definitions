#!/usr/bin/env python3
"""Inject a centered page-number footer into a pandoc-generated .docx.

Pandoc does not emit page numbers. This adds word/footer1.xml with a PAGE
field, registers it in [Content_Types].xml and the document relationships, and
references it from the body sectPr. Idempotent: re-running on an already-patched
file is a no-op. Usage: python add_page_numbers.py <file.docx>
"""
import re
import shutil
import sys
import zipfile

RID = "rId1000"
FOOTER_PART = "word/footer1.xml"
FOOTER_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
    '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="end"/></w:r>'
    "</w:p></w:ftr>"
)
FOOTER_REF = f'<w:footerReference w:type="default" r:id="{RID}"/>'
CONTENT_TYPE = (
    '<Override PartName="/word/footer1.xml" '
    'ContentType="application/vnd.openxmlformats-officedocument.'
    'wordprocessingml.footer+xml"/>'
)
RELATIONSHIP = (
    f'<Relationship Type="http://schemas.openxmlformats.org/officeDocument/'
    f'2006/relationships/footer" Id="{RID}" Target="footer1.xml"/>'
)


def patch(path):
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        data = {n: z.read(n) for n in names}

    if FOOTER_PART in data:
        print("already has a footer; nothing to do")
        return

    doc = data["word/document.xml"].decode("utf-8")
    if "<w:sectPr>" not in doc:
        sys.exit("no <w:sectPr> found in document.xml; cannot attach footer")
    doc = doc.replace("<w:sectPr>", "<w:sectPr>" + FOOTER_REF, 1)
    data["word/document.xml"] = doc.encode("utf-8")

    ct = data["[Content_Types].xml"].decode("utf-8")
    ct = ct.replace("</Types>", CONTENT_TYPE + "</Types>")
    data["[Content_Types].xml"] = ct.encode("utf-8")

    rels = data["word/_rels/document.xml.rels"].decode("utf-8")
    rels = rels.replace("</Relationships>", RELATIONSHIP + "</Relationships>")
    data["word/_rels/document.xml.rels"] = rels.encode("utf-8")

    data[FOOTER_PART] = FOOTER_XML.encode("utf-8")

    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for name, payload in data.items():
            z.writestr(name, payload)
    shutil.move(tmp, path)
    print(f"added page-number footer to {path}")


if __name__ == "__main__":
    patch(sys.argv[1] if len(sys.argv) > 1 else "molecular_informatics_submission.docx")
