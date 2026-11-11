"""Build current Word/PDF manuscripts, tables and submission text from source."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
import textwrap
from datetime import date

from docx_tools import number_table, prepare_main, validate_main

HERE = Path(__file__).resolve().parent


def run(command, **kwargs):
    return subprocess.run(command, cwd=HERE, check=True, text=True, **kwargs)


def plain(latex):
    return run(["pandoc", "-f", "latex", "-t", "plain", "--wrap=none"],
               input=latex, capture_output=True).stdout.strip()


def for_word(latex):
    # Pandoc can expose adjacent booktabs cmidrule arguments as cell text.
    # Rules are layout instructions, so omit them from the Word input only.
    latex = re.sub(r'([A-Za-z]+)\$_(\{[^$]*?\}|[0-9]+)\$',
                   lambda m: '$\\mathrm{' + m[1] + '}_' + m[2] + '$', latex)
    return re.sub(r'\\cmidrule(?:\([^)]*\))?\{[^}]*\}', '', latex)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--word-only", action="store_true")
    args = parser.parse_args()
    output = HERE / f"molecular_informatics_submission_revised_{date.today().isoformat()}.docx"
    version = 2
    while output.exists():
        output = HERE / f"molecular_informatics_submission_revised_{date.today().isoformat()}_{version}.docx"
        version += 1
    exports = HERE / 'submission_exports'
    exports.mkdir(exist_ok=True)
    source = (HERE / "main.tex").read_text()
    parts = [HERE / (name + ".tex") for name in re.findall(r'\\input\{([^}]+)\}', source)]
    full = "\n".join(path.read_text() for path in parts)
    abstract = " ".join(plain(re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', full, re.S)[1]).split())
    assert 600 <= len(abstract) <= 1000, len(abstract)
    (HERE / "abstract_plain.txt").write_text(textwrap.fill(abstract, width=78) + "\n")
    tables = re.findall(r'\\begin\{table\}.*?\\end\{table\}', full, re.S)
    assert len(tables) == 4
    for number, table in enumerate(tables, 1):
        wrapper = ("\\documentclass{article}\n\\usepackage[letterpaper,margin=1in]{geometry}\n"
                   "\\usepackage{booktabs,amsmath,amssymb}\n"
                   "\\pagestyle{empty}\n\\begin{document}\n" +
                   f"\\setcounter{{table}}{{{number-1}}}\n" + table + "\n\\end{document}\n")
        path = HERE / f"submission_tables/table{number}.tex"
        path.write_text(wrapper)
        table_output = exports / f"Table{number}.docx"
        run(["pandoc", "-f", "latex", "-o", str(table_output)], input=for_word(wrapper))
        number_table(table_output, number)
    legends = []
    for number, figure in enumerate(re.findall(r'\\begin\{figure\}.*?\\end\{figure\}', full, re.S), 1):
        caption = figure.split("\\caption{", 1)[1].rsplit("\\label", 1)[0].rstrip()
        assert caption.endswith("}")
        legends.append(f"Figure {number}. {plain(caption[:-1])}")
    assert len(legends) == 3
    (HERE / "FIGURE_LEGENDS.txt").write_text("Figure legends\n\n" + "\n\n".join(legends) + "\n")
    with tempfile.TemporaryDirectory(prefix="kinase-docx-") as temporary:
        bibliography = Path(temporary) / "references.bib"
        bibliography.write_text((HERE / "references.bib").read_text().replace("{et~al.}", "others"))
        word_source = source
        for path in parts:
            relative = str(path.relative_to(HERE).with_suffix(''))
            word_source = word_source.replace('\\input{'+relative+'}', path.read_text())
        word_tex = Path(temporary)/'main_word.tex'
        word_tex.write_text(for_word(word_source))
        run(["pandoc", str(word_tex), "--resource-path=.:sections", "--citeproc",
             "--bibliography", str(bibliography), "--csl", "angewandte-chemie.csl",
             "--metadata", "reference-section-title=References",
             "--lua-filter", "runin-headings.lua", "-o", str(output)])
        run(["python3", "add_author_block.py", str(output)])
        prepare_main(output)
        toc = (HERE / "graphical_abstract_text.txt").read_text().strip().split("\n\n")
        toc_body = " ".join(toc[-1].split())
        assert len(toc_body.split()) <= 80
        toc_source = Path(temporary) / "toc.tex"
        toc_source.write_text("\\section*{" + toc[0] + "}\n" + toc[1] + "\n\n" + toc_body +
                              "\n\n\\includegraphics[width=5cm]{figures/toc_graphic.png}\n")
        run(["pandoc", str(toc_source), "--resource-path=.", "-o", str(exports/'graphical_abstract.docx')])
    checks = validate_main(output)
    checks.update(document=output.name, exports_directory='submission_exports')
    checks.update(abstract_characters=len(abstract), toc_words=len(toc_body.split()))
    checks["source_hashes"] = {str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [HERE / "main.tex", *parts, HERE / "references.bib", HERE / "supporting_information.tex",
                     HERE / "response_to_reviewers_round2.txt", HERE / "cover_letter.txt"]}
    if not args.word_only:
        for document in ("main.tex", "supporting_information.tex"):
            run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", document],
                stdout=subprocess.DEVNULL)
        for number in range(1, 5):
            run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                 "-output-directory=submission_tables", f"submission_tables/table{number}.tex"],
                stdout=subprocess.DEVNULL)
    (HERE / "submission_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(f"Built current submission: {len(abstract)} abstract characters; 4 tables; "
          f"3 figures; {checks['bibliography_entries']} references; {len(toc_body.split())} TOC words.")
    print(f"New manuscript: {output}")


if __name__ == "__main__":
    main()
