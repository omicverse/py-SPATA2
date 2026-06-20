from __future__ import annotations

import argparse
import csv
import inspect
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import spata2py


def parse_exports(namespace_path: Path) -> list[str]:
    text = namespace_path.read_text(encoding="utf-8")
    return re.findall(r"^export\(([^)]+)\)", text, flags=re.MULTILINE)


def category(symbol: str) -> str:
    prefixes = [
        "get",
        "set",
        "add",
        "create",
        "plot",
        "ggp",
        "compute",
        "contains",
        "identify",
        "remove",
        "read",
        "run",
        "as",
        "check",
        "is",
    ]
    for prefix in prefixes:
        if symbol.startswith(prefix):
            return prefix
    if symbol[:1].isupper():
        return "class-or-constructor"
    return "other"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--namespace", required=True, type=Path)
    parser.add_argument("--csv-out", required=True, type=Path)
    parser.add_argument("--md-out", required=True, type=Path)
    args = parser.parse_args()

    exports = parse_exports(args.namespace)
    public = {
        name
        for name, value in vars(spata2py).items()
        if not name.startswith("_") and (callable(value) or inspect.isclass(value))
    }

    rows = []
    for symbol in exports:
        implemented = symbol in public
        rows.append(
            {
                "r_export": symbol,
                "category": category(symbol),
                "python_symbol": symbol if implemented else "",
                "status": "implemented" if implemented else "pending",
            }
        )

    args.csv_out.parent.mkdir(parents=True, exist_ok=True)
    with args.csv_out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["r_export", "category", "python_symbol", "status"])
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    implemented = sum(row["status"] == "implemented" for row in rows)
    categories: dict[str, tuple[int, int]] = {}
    for row in rows:
        done, count = categories.get(row["category"], (0, 0))
        categories[row["category"]] = (done + (row["status"] == "implemented"), count + 1)

    lines = [
        "# SPATA2 Namespace Parity",
        "",
        "Source NAMESPACE: `theMILOlab/SPATA2` v3.1.4 (`NAMESPACE` on master at audit time)",
        "",
        f"- R exports audited: {total}",
        f"- Python symbols implemented with matching R export names: {implemented}",
        f"- Coverage: {implemented / total:.1%}",
        "",
        "This table is intentionally strict: a symbol is counted as implemented only",
        "when the Python package exposes the same exported SPATA2 name.",
        "",
        "## Coverage By Category",
        "",
        "| Category | Implemented | Total | Coverage |",
        "|---|---:|---:|---:|",
    ]
    for name, (done, count) in sorted(categories.items()):
        lines.append(f"| {name} | {done} | {count} | {done / count:.1%} |")

    lines.extend(
        [
            "",
            "## Implemented R-Compatible Symbols",
            "",
        ]
    )
    for row in rows:
        if row["status"] == "implemented":
            lines.append(f"- `{row['r_export']}`")

    lines.extend(
        [
            "",
            "## Pending Symbols",
            "",
            "See `references/spata2_namespace_parity.csv` for the full symbol-level table.",
        ]
    )
    args.md_out.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
