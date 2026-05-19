#!/usr/bin/env python3
import argparse
import csv
import re
import sys
from io import StringIO
from pathlib import Path
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen


GENE_ALIASES = ["gene", "locusid", "locus_id", "locus", "geneid"]
GROUP_ALIASES = ["group", "source", "organization", "org", "lab"]
SCORE_ALIASES = [
    "categoricalscore",
    "categorical_score",
    "classification",
    "score",
    "category",
]


def normalize(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def extract_sheet_parts(sheet_url: str) -> tuple[str, str]:
    match = re.search(r"/spreadsheets/d/([a-zA-Z0-9-_]+)", sheet_url)
    if not match:
        raise ValueError("Could not find spreadsheet ID in URL.")
    spreadsheet_id = match.group(1)

    parsed = urlparse(sheet_url)
    gid = "0"

    query_gid = parse_qs(parsed.query).get("gid")
    if query_gid and query_gid[0]:
        gid = query_gid[0]

    if parsed.fragment:
        frag_match = re.search(r"gid=(\d+)", parsed.fragment)
        if frag_match:
            gid = frag_match.group(1)

    return spreadsheet_id, gid


def build_export_url(sheet_url: str) -> str:
    spreadsheet_id, gid = extract_sheet_parts(sheet_url)
    return f"https://docs.google.com/spreadsheets/d/{spreadsheet_id}/export?format=csv&gid={gid}"


def download_csv_text(export_url: str) -> str:
    req = Request(
        export_url,
        headers={
            "User-Agent": "Mozilla/5.0",
            "Accept": "text/csv,application/octet-stream;q=0.9,*/*;q=0.8",
        },
    )
    with urlopen(req) as response:
        return response.read().decode("utf-8-sig")


def pick_column(headers: list[str], aliases: list[str], required_for: str) -> str:
    normalized_to_original = {normalize(h): h for h in headers}

    for alias in aliases:
        if alias in normalized_to_original:
            return normalized_to_original[alias]

    for alias in aliases:
        for norm_header, original in normalized_to_original.items():
            if alias in norm_header:
                return original

    raise ValueError(
        f"Could not detect the '{required_for}' column from headers: {headers}"
    )


def detect_header_row(rows: list[list[str]]) -> tuple[int, list[str]]:
    for idx, row in enumerate(rows):
        headers = [cell.strip() for cell in row]
        if not any(headers):
            continue

        try:
            pick_column(headers, aliases=GENE_ALIASES, required_for="Gene")
            pick_column(headers, aliases=GROUP_ALIASES, required_for="Group")
            pick_column(
                headers,
                aliases=SCORE_ALIASES,
                required_for="categorical_score",
            )
            return idx, headers
        except ValueError:
            continue

    preview = [
        [cell.strip() for cell in row]
        for row in rows[:5]
    ]
    raise ValueError(
        "Could not find a valid header row with Gene/Group/Score columns. "
        f"First rows preview: {preview}"
    )


def clean_score(raw_score: str) -> str:
    score = raw_score.strip()
    score_map = {
        "refuted": "Contradictory",
        "disputed": "Contradictory",
        "contradictory": "Contradictory",
        "limited": "Limited",
        "moderate": "Moderate",
        "strong": "Strong",
        "supportive": "Supportive",
        "definitive": "Definitive",
        "defintive": "Definitive",
    }
    return score_map.get(score.lower(), score)


def parse_wide_matrix(rows: list[list[str]]) -> list[dict[str, str]]:
    descriptor_idx = None

    for idx, row in enumerate(rows):
        normalized = [normalize(cell) for cell in row]
        if "categoricalscore" in normalized and "gene" in normalized:
            descriptor_idx = idx
            break

    if descriptor_idx is None:
        raise ValueError(
            "Could not detect wide-format descriptor row containing both "
            "'Categorical Score' and 'Gene'."
        )

    if descriptor_idx == 0:
        raise ValueError("Wide-format descriptor row found, but group row is missing.")

    group_row = [cell.strip() for cell in rows[descriptor_idx - 1]]
    descriptor_row = [cell.strip() for cell in rows[descriptor_idx]]
    descriptor_norm = [normalize(cell) for cell in descriptor_row]

    gene_indices = [i for i, value in enumerate(descriptor_norm) if value == "gene"]
    if not gene_indices:
        raise ValueError("Could not locate 'Gene' column in wide-format descriptor row.")
    gene_idx = gene_indices[0]

    score_columns: list[tuple[int, str]] = []
    for col_idx, value in enumerate(descriptor_norm):
        if value != "categoricalscore":
            continue

        group = ""
        search_idx = col_idx
        while search_idx >= 0:
            candidate = group_row[search_idx] if search_idx < len(group_row) else ""
            candidate = candidate.strip()
            if candidate:
                group = candidate
                break
            search_idx -= 1

        if group:
            score_columns.append((col_idx, group))

    if not score_columns:
        raise ValueError("Could not detect any group score columns in wide-format sheet.")

    standardized: list[dict[str, str]] = []
    seen = set()

    for raw_row in rows[descriptor_idx + 1:]:
        if not any(cell.strip() for cell in raw_row):
            continue

        gene = ""
        if gene_idx < len(raw_row):
            gene = raw_row[gene_idx].strip()
        if not gene and raw_row:
            gene = raw_row[0].strip().lstrip("*")
        if not gene:
            continue

        for col_idx, group in score_columns:
            score = raw_row[col_idx].strip() if col_idx < len(raw_row) else ""
            if not score:
                continue

            score = clean_score(score)

            key = (gene, group, score)
            if key in seen:
                continue
            seen.add(key)

            standardized.append(
                {
                    "Gene": gene,
                    "Group": group,
                    "categorical_score": score,
                }
            )

    standardized.sort(key=lambda r: (r["Gene"], r["Group"], r["categorical_score"]))
    return standardized


def standardize_rows(csv_text: str) -> list[dict[str, str]]:
    rows = list(csv.reader(StringIO(csv_text)))
    if not rows:
        raise ValueError("Downloaded sheet appears empty or has no header row.")

    try:
        header_idx, headers = detect_header_row(rows)
    except ValueError:
        return parse_wide_matrix(rows)

    gene_col = pick_column(
        headers,
        aliases=GENE_ALIASES,
        required_for="Gene",
    )
    group_col = pick_column(
        headers,
        aliases=GROUP_ALIASES,
        required_for="Group",
    )
    score_col = pick_column(
        headers,
        aliases=SCORE_ALIASES,
        required_for="categorical_score",
    )

    standardized: list[dict[str, str]] = []
    seen = set()

    for raw_row in rows[header_idx + 1:]:
        if not any(cell.strip() for cell in raw_row):
            continue

        row = {
            headers[i]: (raw_row[i] if i < len(raw_row) else "")
            for i in range(len(headers))
        }

        gene = (row.get(gene_col) or "").strip()
        group = (row.get(group_col) or "").strip()
        score = (row.get(score_col) or "").strip()

        if not gene or not group or not score:
            continue

        score = clean_score(score)

        key = (gene, group, score)
        if key in seen:
            continue
        seen.add(key)

        standardized.append(
            {
                "Gene": gene,
                "Group": group,
                "categorical_score": score,
            }
        )

    standardized.sort(key=lambda r: (r["Gene"], r["Group"], r["categorical_score"]))
    return standardized


def write_output(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Gene", "Group", "categorical_score"])
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download a Google Sheet tab as CSV and format it like criTRia_Dataset.csv "
            "(Gene, Group, categorical_score)."
        )
    )
    parser.add_argument(
        "--sheet-url",
        default="https://docs.google.com/spreadsheets/d/1DonSiPVjeQLsB8HoFzn50jFCnWvMPxVZ124SzS0s7r0/edit?gid=0#gid=0",
        help="Google Sheet URL for the tab you want to export (default: provided sheet).",
    )
    parser.add_argument(
        "--output",
        default="criTRia_Dataset.from_sheet.csv",
        help="Output CSV path (default: criTRia_Dataset.from_sheet.csv).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        export_url = build_export_url(args.sheet_url)
        csv_text = download_csv_text(export_url)
        rows = standardize_rows(csv_text)
        write_output(rows, Path(args.output))
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(
        f"Wrote {len(rows)} rows to {args.output} in format: "
        "Gene,Group,categorical_score"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
