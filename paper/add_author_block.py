#!/usr/bin/env python3
"""Insert affiliation / ORCID / email paragraphs into a pandoc-generated .docx.

Pandoc renders achemso's \\author (the name) but ignores \\affiliation,
\\altaffiliation, and \\email, so the Word file loses them. The PDF keeps them
via achemso's title block; this restores them to the docx only, right after the
author name. Idempotent. Usage: python add_author_block.py <file.docx>
"""
import shutil
import sys
import zipfile

AUTHOR = "Polina Vinogradova"
AFFILIATION = "Independent Researcher"
ORCID_EMAIL = "ORCID: 0000-0003-3271-3841 · polina.vino@gmail.com"

ANCHOR = (
    f'<w:r><w:t xml:space="preserve">{AUTHOR}</w:t></w:r>\n'
    "    </w:p>"
)


def _para(text):
    return (
        '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
        f'<w:r><w:t xml:space="preserve">{text}</w:t></w:r></w:p>'
    )


def patch(path):
    with zipfile.ZipFile(path) as z:
        data = {n: z.read(n) for n in z.namelist()}

    doc = data["word/document.xml"].decode("utf-8")
    if AFFILIATION in doc:
        print("author block already present; nothing to do")
        return
    if ANCHOR not in doc:
        sys.exit("could not locate the author paragraph; aborting")

    block = ANCHOR + "\n" + _para(AFFILIATION) + "\n" + _para(ORCID_EMAIL)
    doc = doc.replace(ANCHOR, block, 1)
    data["word/document.xml"] = doc.encode("utf-8")

    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for name, payload in data.items():
            z.writestr(name, payload)
    shutil.move(tmp, path)
    print(f"added affiliation/ORCID/email block to {path}")


if __name__ == "__main__":
    patch(sys.argv[1] if len(sys.argv) > 1 else "molecular_informatics_submission.docx")
