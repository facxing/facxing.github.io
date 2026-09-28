#!/usr/bin/env python3
"""Render every figure PDF under assets/images/notes to a sibling PNG.

Why this is needed
------------------
The notes embed figures with plain HTML, e.g.

    <img src="/assets/images/notes/FeynmanDiag/WH_associated.pdf">

Browsers cannot display a PDF inside <img>; the element renders as a broken
image.  Pointing the same tag at a PNG fixes it, so every figure PDF gets a
same-named .png next to it and the references are switched to ".png".

Notes
-----
* Rendering uses ghostscript (`gs`) because `pdftoppm` / `pdftocairo` are not
  installed in this environment.
* Page 1 only -- all note figures are single page.
* The PDFs are deliberately kept on disk: they are the vector originals and
  stay useful for download or re-rendering at a different resolution.

Public API
----------
    collect_pdfs(root) -> list[pathlib.Path]
    render(pdf, dpi, force) -> bool
    main() -> int

Usage
-----
    python3 scripts/pdf2png.py                 # convert what is missing/stale
    python3 scripts/pdf2png.py --dpi 300       # higher resolution
    python3 scripts/pdf2png.py --force         # re-render everything
    python3 scripts/pdf2png.py --dry-run       # only report
"""

import argparse
import pathlib
import subprocess
import sys

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════
REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_ROOT = REPO_ROOT / "assets" / "images" / "notes"
DEFAULT_DPI = 200
GS = "gs"


def collect_pdfs(root):
    """Return every *.pdf under ``root``, in a deterministic order.

    Parameters
    ----------
    root : str or pathlib.Path
        Directory to search recursively.

    Returns
    -------
    list[pathlib.Path]
        Sorted list of PDF paths.
    """
    return sorted(pathlib.Path(root).rglob("*.pdf"))


def render(pdf, dpi=DEFAULT_DPI, force=False, dry_run=False):
    """Render page 1 of ``pdf`` to a sibling PNG.

    Parameters
    ----------
    pdf : pathlib.Path
        Source PDF.
    dpi : int
        Rasterisation resolution.
    force : bool
        Re-render even when the PNG is already newer than the PDF.
    dry_run : bool
        Do not write anything, only report what would happen.

    Returns
    -------
    bool
        True if the PNG is (or already was) up to date.
    """
    png = pdf.with_suffix(".png")
    if png.exists() and not force and png.stat().st_mtime >= pdf.stat().st_mtime:
        print("[ skip  ] %s  (up to date)" % png.relative_to(REPO_ROOT))
        return True
    if dry_run:
        print("[ would ] %s" % png.relative_to(REPO_ROOT))
        return True

    cmd = [GS, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
           "-dFirstPage=1", "-dLastPage=1",
           "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
           "-sDEVICE=png16m", "-r%d" % dpi,
           "-sOutputFile=%s" % png, str(pdf)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not png.exists():
        print("[ ERROR ] %s\n          %s"
              % (png.relative_to(REPO_ROOT), (proc.stderr or "").strip()[:200]))
        return False
    print("[  png  ] %-60s %6.1f KB" % (png.relative_to(REPO_ROOT),
                                        png.stat().st_size / 1024.0))
    return True


def main():
    """Entry point: convert every figure PDF under the notes asset tree."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default=str(DEFAULT_ROOT),
                    help="directory to scan (default: assets/images/notes)")
    ap.add_argument("--dpi", type=int, default=DEFAULT_DPI,
                    help="rasterisation resolution (default: %d)" % DEFAULT_DPI)
    ap.add_argument("--force", action="store_true",
                    help="re-render PNGs that already exist")
    ap.add_argument("--dry-run", action="store_true",
                    help="only list what would be converted")
    args = ap.parse_args()

    pdfs = collect_pdfs(args.root)
    if not pdfs:
        print("[ input ] no PDF found under %s" % args.root)
        return 0
    print("[ input ] %d PDF(s) under %s" % (len(pdfs), args.root))

    ok = 0
    for pdf in pdfs:
        ok += bool(render(pdf, dpi=args.dpi, force=args.force, dry_run=args.dry_run))
    print("[ done  ] %d/%d OK" % (ok, len(pdfs)))
    return 0 if ok == len(pdfs) else 1


if __name__ == "__main__":
    sys.exit(main())
