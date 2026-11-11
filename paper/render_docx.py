"""Render current Word uploads with LibreOffice and verify the image count.

The resulting PDF is a visual QA artifact, not a replacement for the DOCX.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('documents', nargs='*', type=Path)
    parser.add_argument('--pdf-dir', type=Path,
                        help='Validate PDFs already exported by the authorized external renderer')
    args = parser.parse_args()
    executable = shutil.which('soffice') or '/Applications/LibreOffice.app/Contents/MacOS/soffice'
    if args.documents:
        documents = args.documents
    else:
        checks = json.loads((HERE/'submission_checks.json').read_text())
        exports = HERE/checks['exports_directory']
        documents = [HERE/checks['document'], *sorted(exports.glob('Table*.docx')),
                     exports/'graphical_abstract.docx']
    results = []
    with tempfile.TemporaryDirectory(prefix='kinase-render-') as temporary:
        temporary = Path(temporary)
        for document in documents:
            document = document.resolve()
            options = json.dumps({'UseLosslessCompression': {'type':'boolean','value':'true'},
                                  'ReduceImageResolution': {'type':'boolean','value':'false'}})
            if args.pdf_dir:
                rendered = args.pdf_dir/(document.stem+'.pdf')
            else:
                subprocess.run([executable, '-env:UserInstallation='+(temporary/'profile').as_uri(),
                                '--headless', '--convert-to', 'pdf:writer_pdf_Export:'+options,
                                '--outdir', str(temporary), str(document)], check=True)
                rendered = temporary/(document.stem+'.pdf')
            assert rendered.exists(), f'No rendered PDF for {document}'
            target = document.with_name(document.stem+'_rendered.pdf')
            shutil.copy2(rendered, target)
            listing = subprocess.run(['pdfimages','-list',str(target)],text=True,capture_output=True,check=True).stdout
            images = len(re.findall(r'^\s*\d+\s+\d+\s+image\s',listing,re.M))
            expected = 3 if document.stem.startswith('molecular_informatics_submission') else (1 if document.stem=='graphical_abstract' else 0)
            assert images == expected, f'{document.name}: rendered {images} images, expected {expected}'
            info = subprocess.run(['pdfinfo',str(target)],text=True,capture_output=True,check=True).stdout
            results.append(dict(document=str(document.relative_to(HERE)),pdf=str(target.relative_to(HERE)),
                                images=images,pages=int(re.search(r'Pages:\s*(\d+)',info)[1]),
                                document_sha256=hashlib.sha256(document.read_bytes()).hexdigest(),
                                pdf_sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    (HERE/'render_checks.json').write_text(json.dumps(results,indent=2)+'\n')
    print('All Word uploads rendered; image counts verified.')


if __name__=='__main__':
    main()
