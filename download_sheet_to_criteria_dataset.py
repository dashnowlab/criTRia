#!/usr/bin/env python3
"""Build criTRia_Dataset.csv from GenCC and the criTRia curations.

For every STRchive locus (plus any criTRia curation that is not a STRchive
locus), pull the gene-disease classifications for that gene from the GenCC
submissions export (which includes ClinGen). A
classification is only kept when its disease matches the locus disease, i.e.
one of its MONDO/OMIM/Orphanet IDs is listed for the locus in STRchive.
Some phenotypes were manually reviewed and determined similar enough to be
matched, see EXTRA_DISEASE_IDS.
"""
import argparse
import csv
import json
import sys
from io import StringIO
from pathlib import Path
from urllib.request import Request, urlopen

csv.field_size_limit(sys.maxsize)


CRITRIA_TSV_URL = "https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/criTRia-curations.tsv"
STRCHIVE_JSON_URL = "https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/STRchive-loci.json"
GENCC_TSV_URL = "https://thegencc.org/download/action/submissions-export-tsv"

# GenCC submitter -> Group name used in criTRia_Dataset.csv.
# Both PanelApp instances are reported as one "PanelApp" group.
GENCC_GROUPS = {
    "ClinGen": "ClinGen",
    "Ambry Genetics": "Ambry",
    "G2P": "G2P",
    "Genomics England PanelApp": "PanelApp",
    "PanelApp Australia": "PanelApp",
    "Illumina": "Illumina",
    "Labcorp Genetics (formerly Invitae)": "Labcorp",
    "Laboratory for Molecular Medicine": "Lab MM",
    "Myriad Women's Health": "Myriad",
    "Orphanet": "Orphanet",
}

# criTRia curation IDs (Disease_ID_Gene) that are named differently in STRchive.
CRITRIA_TO_STRCHIVE = {
    "ALS1_NIPA": "ALS1_NIPA1",
    "HFG-I_HOXA13": "HFG_HOXA13-I",
    "HFG-II_HOXA13": "HFG_HOXA13-II",
    "HFG-III_HOXA13": "HFG_HOXA13-III",
    "XLID, PHPX_SOX3": "XLID_SOX3",
}

# STRchive lumps all FMR1 diseases into FXS_FMR1, but criTRia curates FXS
# separately from FXTAS/POF1, so give each its own disease IDs.
LOCUS_OVERRIDES = {
    "FXS_FMR1": {
        "gene": "FMR1",
        "disease": "Fragile X syndrome",
        "disease_ids": {"MONDO:0010383", "OMIM:300624", "Orphanet:908"},
    },
    "FXTAS,POF1_FMR1": {
        "gene": "FMR1",
        "disease": "Fragile X-associated tremor/ataxia syndrome; premature ovarian failure 1",
        "disease_ids": {
            "MONDO:0010382", "MONDO:0010706",
            "OMIM:300623", "OMIM:311360",
            "Orphanet:93256", "Orphanet:642691",
        },
    },
}

# Other gene symbols GenCC uses for a locus gene.
GENE_SYMBOL_ALIASES: dict[str, set[str]] = {
    # GenCC files SCA8 under ATXN8, the gene on the opposite strand
    "ATXN8OS": {"ATXN8"},
}

# Extra disease IDs to accept as the same disease as a locus, e.g. when GenCC
# submitters use a lumped or newer MONDO term. Review gene_disease_matches.tsv
# (status == disease_mismatch) before adding anything here.
EXTRA_DISEASE_IDS: dict[str, set[str]] = {
    # Parent term of STRchive's MONDO:0100340 (Friedreich ataxia 1) (also in STRchive PR #523)
    "FRDA_FXN": {"MONDO:0100339"},
    # ClinGen lumped term covering XLID with GH deficiency and PHPX (also in STRchive PR #523)
    "XLID_SOX3": {"MONDO:0800474"},
    # ClinGen/G2P lumped term covering MED1 and pseudoachondroplasia
    "EDM1-PSACH_COMP": {"MONDO:0100593"},
    # Orphanet splits Machado-Joseph disease into clinical types 1-3
    "SCA3_ATXN3": {"MONDO:0017174", "MONDO:0017175", "MONDO:0017176"},
    # Gerstmann-Straussler-Scheinker syndrome, part of the locus disease (also in STRchive PR #523)
    "CJD_PRNP": {"MONDO:0007656", "OMIM:137440"},
    # Duchenne and Becker MD (obsolete MONDO) and ClinGen's lumped term
    "DMD_DMD": {"MONDO:0016899", "MONDO:0016106"},
    # Orphanet has no GIPC1-specific OPDM term
    "OPDM2_GIPC1": {"MONDO:0025193"},
    # ClinGen non-syndromic X-linked intellectual disability
    "FRAXE_AFF2": {"MONDO:0019181"},
    # G2P autism, vs STRchive's autism spectrum disorder (MONDO:0005258)
    "FRA7A_ZNF713": {"MONDO:0005260", "OMIM:209850"},
    # SCA27A (FGF14 coding variants) overlaps SCA27B clinically
    "SCA27B_FGF14": {"MONDO:0012247", "MONDO:0008654"},
    # General AR/AD progressive external ophthalmoplegia (parents of STRchive's POLG-specific terms)
    "CPEO_POLG": {"MONDO:0016810", "MONDO:0008003"},
    # Conotruncal heart malformations, a group that includes tetralogy of Fallot
    "TOF_TBX1": {"MONDO:0016581", "OMIM:217095"},
    # ClinGen lumped term; its NAXE evidence is all NAXE-related encephalopathy
    "NME_NAXE": {"MONDO:0044970"},
    # G2P intellectual disability from AFF3 loss of function (not KINSSHIP syndrome)
    "FRA2A_AFF3": {"MONDO:0001071"},
    # ClinGen lumped ARX into two terms; assign one to each polyalanine tract
    "EIEE1_ARX": {"MONDO:0100062"},  # genetic developmental and epileptic encephalopathy
    "PRTS_ARX": {"MONDO:0100148"},  # X-linked complex neurodevelopmental disorder
}

# Disease IDs that must never match a locus, even if STRchive lists them.
# A record carrying any of these IDs is reported with status "excluded".
EXCLUDED_DISEASE_IDS: dict[str, set[str]] = {
    # X-linked lissencephaly with abnormal genitalia (XLAG) is a different
    # disease from early infantile epileptic encephalopathy
    # (OMIM:300215 is removed from EIEE1 in STRchive PR #523)
    "EIEE1_ARX": {"OMIM:300215", "MONDO:0010268", "Orphanet:452"},
}

# Breaks ties between a group's matching records curated on the same date.
# Refuted and Disputed outrank everything, so a contradicting record is never
# hidden by a positive one from the same date.
SCORE_RANK = {
    "Refuted": 9,
    "Disputed": 8,
    "Definitive": 7,
    "Strong": 6,
    "Supportive": 5,
    "Moderate": 4,
    "Limited": 3,
    "No Known": 0,
}

REPORT_FIELDS = [
    "status", "locus", "gene", "locus_disease", "locus_disease_ids",
    "source", "submitter", "classification", "moi",
    "disease_id", "disease_name", "original_disease_id", "original_disease_name",
    "date",
]


def download_text(url: str) -> str:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req) as response:
        return response.read().decode("utf-8-sig")


def clean_score(raw_score: str) -> str:
    score = raw_score.strip()
    score_map = {
        "refuted": "Contradictory",
        "refuted evidence": "Contradictory",
        "disputed": "Contradictory",
        "disputed evidence": "Contradictory",
        "contradictory": "Contradictory",
        "limited": "Limited",
        "moderate": "Moderate",
        "strong": "Strong",
        "supportive": "Supportive",
        "definitive": "Definitive",
        "defintive": "Definitive",
        "no known": "No Known",
        "no known disease relationship": "No Known",
    }
    return score_map.get(score.lower(), score)


def rank_score(raw_score: str) -> int:
    score = raw_score.replace(" Evidence", "").replace(" Disease Relationship", "")
    return SCORE_RANK.get(score, -1)


def parse_strchive(json_text: str) -> dict[str, dict]:
    loci = {}
    for locus in json.loads(json_text):
        disease_ids = (
            {f"MONDO:{i}" for i in locus.get("mondo") or []}
            | {f"OMIM:{i}" for i in locus.get("omim") or []}
            | {f"Orphanet:{i}" for i in locus.get("orphanet") or []}
        )
        loci[locus["id"]] = {
            "gene": locus["gene"],
            "disease": locus.get("disease") or "",
            "disease_ids": disease_ids,
        }
    return loci


def parse_criTRia_tsv(tsv_text: str) -> list[dict[str, str]]:
    reader = csv.DictReader(StringIO(tsv_text), delimiter="\t")
    rows = []
    seen = set()
    for row in reader:
        gene = (row.get("Gene") or "").strip()
        disease_id = (row.get("Disease_ID") or "").strip()
        classification = (row.get("classification") or "").strip()
        if not gene or not disease_id or not classification:
            continue
        locus = f"{disease_id}_{gene}"
        locus = CRITRIA_TO_STRCHIVE.get(locus, locus)
        score = clean_score(classification)
        key = (locus, score)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"Gene": locus, "Group": "criTRia", "categorical_score": score})
    return rows


def parse_gencc_tsv(tsv_text: str) -> list[dict[str, str]]:
    records = []
    for row in csv.DictReader(StringIO(tsv_text, newline=""), delimiter="\t"):
        submitter = row["submitter_title"]
        records.append({
            "source": "GenCC",
            "submitter": submitter,
            "group": GENCC_GROUPS.get(submitter, ""),
            "gene": row["gene_symbol"],
            "classification": row["classification_title"],
            "moi": row["moi_title"],
            "disease_id": row["disease_curie"],
            "disease_name": row["disease_title"],
            "original_disease_id": row["disease_original_curie"],
            "original_disease_name": row["disease_original_title"],
            "date": row["submitted_as_date"][:10],
        })
    return records


def build_loci(strchive: dict[str, dict], criTRia_rows: list[dict[str, str]]) -> dict[str, dict]:
    loci = {locus_id: dict(info) for locus_id, info in strchive.items()}
    for locus_id, info in LOCUS_OVERRIDES.items():
        loci[locus_id] = dict(info)
    for locus_id, extra_ids in EXTRA_DISEASE_IDS.items():
        loci[locus_id]["disease_ids"] = loci[locus_id]["disease_ids"] | extra_ids
    for locus_id, excluded_ids in EXCLUDED_DISEASE_IDS.items():
        loci[locus_id]["disease_ids"] = loci[locus_id]["disease_ids"] - excluded_ids
        loci[locus_id]["excluded_ids"] = excluded_ids
    for row in criTRia_rows:
        if row["Gene"] not in loci:
            print(
                f"Warning: criTRia locus {row['Gene']} not found in STRchive; "
                "add it to CRITRIA_TO_STRCHIVE or LOCUS_OVERRIDES.",
                file=sys.stderr,
            )
    return loci


def match_records(loci: dict[str, dict], records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Return every record for a locus gene, labelled match, excluded or disease_mismatch."""
    by_gene: dict[str, list[dict[str, str]]] = {}
    for record in records:
        by_gene.setdefault(record["gene"], []).append(record)

    report = []
    for locus_id, locus in sorted(loci.items()):
        genes = {locus["gene"]} | GENE_SYMBOL_ALIASES.get(locus["gene"], set())
        for record in (r for gene in sorted(genes) for r in by_gene.get(gene, [])):
            ids = {record["disease_id"], record["original_disease_id"]} - {""}
            if ids & locus.get("excluded_ids", set()):
                status = "excluded"
            elif ids & locus["disease_ids"]:
                status = "match"
            else:
                status = "disease_mismatch"
            report.append({
                "status": status,
                "locus": locus_id,
                "locus_disease": locus["disease"],
                "locus_disease_ids": ";".join(sorted(locus["disease_ids"])),
                **record,
            })
    return report


def summarise_matches(report: list[dict[str, str]]) -> list[dict[str, str]]:
    """Collapse matched records to one score per locus and group.

    The most recent curation wins; records curated on the same date are ranked
    by SCORE_RANK.
    """
    best: dict[tuple[str, str], dict[str, str]] = {}
    for row in report:
        if row["status"] != "match" or not row["group"]:
            continue
        key = (row["locus"], row["group"])
        order = (row["date"], rank_score(row["classification"]))
        current = best.get(key)
        if current is None or order > (current["date"], rank_score(current["classification"])):
            best[key] = row
    return [
        {
            "Gene": locus,
            "Group": group,
            "categorical_score": clean_score(row["classification"]),
        }
        for (locus, group), row in best.items()
    ]


def write_csv(rows: list[dict[str, str]], fieldnames: list[str], output_path: Path, delimiter: str = ",") -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore", delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download GenCC and criTRia classifications and format them like "
            "criTRia_Dataset.csv (Gene, Group, categorical_score)."
        )
    )
    parser.add_argument(
        "--output",
        default="criTRia_Dataset.csv",
        help="Output CSV path (default: criTRia_Dataset.csv).",
    )
    parser.add_argument(
        "--report",
        default="gene_disease_matches.tsv",
        help=(
            "TSV listing every GenCC record for each locus gene and whether "
            "its disease matches the locus (default: gene_disease_matches.tsv)."
        ),
    )
    parser.add_argument("--criTRia-url", default=CRITRIA_TSV_URL, help="criTRia-curations TSV URL.")
    parser.add_argument(
        "--criTRia-copy",
        default="criTRia-curations.tsv",
        help="Where to save a copy of the downloaded criTRia curations (default: criTRia-curations.tsv).",
    )
    parser.add_argument("--strchive-url", default=STRCHIVE_JSON_URL, help="STRchive loci JSON URL.")
    parser.add_argument("--gencc-url", default=GENCC_TSV_URL, help="GenCC submissions TSV URL.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        criTRia_text = download_text(args.criTRia_url)
        criTRia_rows = parse_criTRia_tsv(criTRia_text)
        loci = build_loci(parse_strchive(download_text(args.strchive_url)), criTRia_rows)
        gencc = parse_gencc_tsv(download_text(args.gencc_url))
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    with open(args.criTRia_copy, "w", encoding="utf-8", newline="") as f:
        f.write(criTRia_text)
    print(f"Saved criTRia curations from {args.criTRia_url} to {args.criTRia_copy}")

    report = match_records(loci, gencc)
    write_csv(report, REPORT_FIELDS + ["group"], Path(args.report), delimiter="\t")

    rows = summarise_matches(report) + criTRia_rows
    rows.sort(key=lambda r: (r["Gene"], r["Group"], r["categorical_score"]))
    write_csv(rows, ["Gene", "Group", "categorical_score"], Path(args.output))
    print(f"Wrote {len(rows)} rows to {args.output} in format: Gene,Group,categorical_score")

    mismatches = [r for r in report if r["status"] == "disease_mismatch" and r["group"]]
    print(
        f"Wrote {len(report)} records to {args.report}; "
        f"{len(mismatches)} from included groups have a disease that does not match the locus."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
