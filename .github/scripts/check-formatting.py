#!/usr/bin/env python3
"""Fail on the formatting BiocCheck flags in hand-written sources.

Mirrors BiocCheck's checkFormatting(): lines over 80 characters, tabs, and
leading-space indents that are not a multiple of 4. Only R/ and vignettes/
are checked. man/ and NAMESPACE are roxygen output, so their hits are fixed in
the roxygen comments or upstream, not here.

Usage: python3 .github/scripts/check-formatting.py  (from the package root)
"""
import pathlib
import re
import sys

MAX_CHARS = 80
SOURCES = {"R": ("*.R",), "vignettes": ("*.Rmd", "*.qmd")}


def stripFrontMatter(lines):
    """Blank a leading YAML block, which BiocCheck skips for tabs/indent."""
    if not lines or lines[0].strip() != "---":
        return lines
    for end in range(1, len(lines)):
        if lines[end].strip() == "---":
            return [""] * (end + 1) + lines[end + 1:]
    return lines


def checkFile(path):
    lines = path.read_text(encoding="utf-8").split("\n")
    problems = []
    for num, line in enumerate(lines, 1):
        if len(line) > MAX_CHARS:
            problems.append((num, f"{len(line)} characters (max {MAX_CHARS})"))
    body = stripFrontMatter(lines) if path.suffix != ".R" else lines
    for num, line in enumerate(body, 1):
        if "\t" in line:
            problems.append((num, "tab character"))
        indent = len(re.match(r"^ *", line).group(0))
        if indent % 4:
            problems.append((num, f"indent of {indent} spaces"))
    return problems


def main():
    failed = 0
    for folder, patterns in SOURCES.items():
        for pattern in patterns:
            for path in sorted(pathlib.Path(folder).glob(pattern)):
                for num, msg in sorted(checkFile(path)):
                    print(f"{path}:{num}: {msg}")
                    failed += 1
    if failed:
        print(f"\n{failed} formatting problem(s); see "
              "https://contributions.bioconductor.org/r-code.html")
        return 1
    print("Formatting OK: R/ and vignettes/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
