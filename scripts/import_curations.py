# Import criTRia curations from a curator's Google Sheet (one locus per tab).
# Output format should match:
# https://github.com/dashnowlab/STRchive/blob/main/data/criTRia-curations.schema.json
# Example:
# https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/criTRia-curations.json
#
# Each tab holds one locus. Row 1 is a header of "Key: value" cells, row 2 is the
# evidence table header, and the remaining rows are evidence entries.
# Generated fields (scores, categories, references, publication dates) are left
# for check-curations.py to fill in, e.g.:
#   python scripts/import_curations.py --sheet 'https://docs.google.com/spreadsheets/d/<id>/edit' -o new-curations.json
#   python scripts/check-curations.py --json new-curations.json --out new-curations.json
import sys
import argparse
import json
import math
import re
import pandas as pd

def sheet_to_xlsx_url(sheet):
    """Turn a Google Sheets link or bare sheet ID into an xlsx export URL; leave file paths and other URLs unchanged.
    >>> sheet_to_xlsx_url('https://docs.google.com/spreadsheets/d/abc_123-X/edit?usp=sharing')
    'https://docs.google.com/spreadsheets/d/abc_123-X/export?format=xlsx'
    >>> sheet_to_xlsx_url('1DF2NPHb-CHq2XrvoyxNU_5yAuHo_szQr7SubEEfddMo')
    'https://docs.google.com/spreadsheets/d/1DF2NPHb-CHq2XrvoyxNU_5yAuHo_szQr7SubEEfddMo/export?format=xlsx'
    >>> sheet_to_xlsx_url('curations.xlsx')
    'curations.xlsx'
    """
    match = re.search(r'docs\.google\.com/spreadsheets/d/([A-Za-z0-9_-]+)', sheet)
    if match:
        sheet_id = match.group(1)
    elif re.fullmatch(r'[A-Za-z0-9_-]{30,}', sheet):
        sheet_id = sheet
    else:
        return sheet
    return f'https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=xlsx'

def parse_args():
    parser = argparse.ArgumentParser(description='Import criTRia curations from the Google Sheet (one locus per tab) into criTRia-curations JSON format.')
    parser.add_argument('--xlsx', '--sheet', dest='xlsx', required=True,
                        help='Curation spreadsheet: a Google Sheets link or sheet ID (must be viewable by anyone with the link), '
                             'or a local/remote xlsx file')
    parser.add_argument('-o', '--out', default='-',
                        help='Output JSON file path (default: stdout)')
    parser.add_argument('--merge',
                        help='Existing criTRia-curations JSON file. Imported loci replace records with the same Locus_ID; all others are kept.')
    parser.add_argument('--tabs', nargs='+',
                        help='Only import these tabs (default: all tabs except templates)')
    parser.add_argument('--source', default='criTRia',
                        help='Value for the Source field (default: criTRia)')
    parser.add_argument('--sop-version', default='criTRia v0',
                        help='Value for the SOP_version field (default: criTRia v0)')
    return parser.parse_args()

# Map header keys in the sheet to schema fields
HEADER_FIELDS = {
    'Disease_ID': 'Disease_ID',
    'Gene': 'Gene',
    'Locus_ID': 'Locus_ID',
    'Inheritance': 'Inheritance',
    'Curator': 'Curator',
    'Contributors': 'Contributors',
    'Date': 'Date',
    'Manual_evidence_level': 'Manual_evidence_level',
    'Summary': 'Description',
    'Description': 'Description',
}

# Map (Evidence Category, Evidence Type) in the sheet to the schema's Evidence type.
# The category is needed because e.g. "Cell Culture" means different things under Models and Rescue.
GENETIC_EVIDENCE = {
    ('singular evidence', 'probands'): 'Probands',
    ('collective evidence', 'allele'): 'Allele',
    ('collective evidence', 'computational'): 'Computational',
    ('collective evidence', 'segregation'): 'Segregation',
    ('statistics', 'case-control'): 'Case-control data',
    ('statistics', 'case-control data'): 'Case-control data',
}
EXPERIMENTAL_EVIDENCE = {
    ('function', 'biochemical function'): 'Biochemical function',
    ('function', 'protein interaction'): 'Protein interaction',
    ('function', 'regulatory impact'): 'Regulatory impact',
    ('functional alteration', 'patient cells'): 'Patient cells',
    ('functional alteration', 'non-patient cells'): 'Non-patient cells',
    ('models', 'non-human model organism'): 'Non-human model organism',
    ('models', 'cell culture'): 'Cell culture',
    ('rescue', 'human treatment'): 'Human treatment',
    ('rescue', 'non-human model organism'): 'Rescue in non-human model organism',
    ('rescue', 'cell culture'): 'Rescue in cell culture',
    ('rescue', 'patient cells'): 'Rescue in patient cells',
}

def clean(value):
    """Return a stripped string, or None for empty/NaN cells.
    >>> clean('  x ')
    'x'
    >>> clean(float('nan')) is None
    True
    >>> clean('') is None
    True
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    value = str(value).strip()
    return value or None

def is_placeholder(value):
    """True if a header value is still template text, e.g. '[AD, AR, XR, XD]', 'DISEASE1' or '10/DD/2026'.
    >>> is_placeholder('[AD, AR, XR, XD]')
    True
    >>> is_placeholder('DISEASE1_GENE1')
    True
    >>> is_placeholder('10/DD/2026')
    True
    >>> is_placeholder('OPDM_FAM193B')
    False
    """
    return bool(re.search(r'^\[.*\]$|DISEASE1|GENE1|MM/DD|DD/|/YYYY', value))

def parse_score(value, tab, evidence_type):
    """Convert a score cell to a number, keeping integers as ints."""
    value = clean(value)
    if value is None:
        return None
    try:
        score = float(value)
    except ValueError:
        sys.stderr.write(f"WARNING: tab '{tab}': non-numeric score '{value}' for {evidence_type}, setting to null\n")
        return None
    return int(score) if score.is_integer() else score

def normalize_citation(value):
    """Tidy separators in a citation cell while keeping multiple citations as one '; '-separated string.
    >>> normalize_citation('pmid:1;  pmid:2 ;')
    'pmid:1; pmid:2'
    >>> normalize_citation('PMID: 123, doi:10.1/x')
    'pmid:123; doi:10.1/x'
    """
    value = clean(value)
    if value is None:
        return None
    parts = [p.strip() for p in re.split(r'[;,\n]', value) if p.strip()]
    tidied = []
    for part in parts:
        match = re.match(r'^([A-Za-z]+)\s*:\s*(.+)$', part)
        if match:
            part = f'{match.group(1).lower()}:{match.group(2).strip()}'
        tidied.append(part)
    return '; '.join(tidied) or None

def parse_header(row, tab):
    """Parse the 'Key: value' cells in the first row of a tab."""
    record = {}
    for cell in row:
        cell = clean(cell)
        if cell is None:
            continue
        key, sep, value = cell.lstrip('#').partition(':')
        key = key.strip()
        if not sep:
            sys.stderr.write(f"WARNING: tab '{tab}': could not parse header cell '{cell}'\n")
            continue
        if key not in HEADER_FIELDS:
            sys.stderr.write(f"WARNING: tab '{tab}': unknown header field '{key}', skipping\n")
            continue
        value = clean(value)
        if value is not None and is_placeholder(value):
            value = None
        record[HEADER_FIELDS[key]] = value
    return record

def parse_tab(df, tab, source, sop_version):
    """Parse one tab into a curation record, or return None if it is an unfilled template."""
    record = parse_header(df.iloc[0], tab)
    if not record.get('Disease_ID') or not record.get('Gene'):
        sys.stderr.write(f"Skipping tab '{tab}': Disease_ID or Gene not filled in\n")
        return None

    expected_locus_id = f"{record['Disease_ID']}_{record['Gene']}"
    if not record.get('Locus_ID'):
        record['Locus_ID'] = expected_locus_id
    elif record['Locus_ID'] != expected_locus_id:
        sys.stderr.write(f"WARNING: tab '{tab}': Locus_ID '{record['Locus_ID']}' does not match '{expected_locus_id}'\n")

    # Find the evidence table header row
    header_idx = None
    for idx in range(len(df)):
        if clean(df.iat[idx, 0]) == 'Evidence Category':
            header_idx = idx
            break
    if header_idx is None:
        sys.stderr.write(f"WARNING: tab '{tab}': no 'Evidence Category' header row found\n")
        header_idx = len(df)
    columns = {clean(c): i for i, c in enumerate(df.iloc[header_idx]) if clean(c)} if header_idx < len(df) else {}

    genetic, experimental = [], []
    for idx in range(header_idx + 1, len(df)):
        row = df.iloc[idx]
        category = clean(row.iloc[columns.get('Evidence Category', 0)])
        sheet_type = clean(row.iloc[columns.get('Evidence Type', 1)])
        citation = normalize_citation(row.iloc[columns.get('Source', 2)])
        score_cell = row.iloc[columns.get('Score', 3)]
        detail = clean(row.iloc[columns.get('Details', 4)])

        # Skip empty template rows (template citation placeholder counts as empty)
        if citation and 'XXXXX' in citation:
            citation = None
        if citation is None and clean(score_cell) is None and detail is None:
            continue

        key = ((category or '').lower(), (sheet_type or '').lower())
        if key in GENETIC_EVIDENCE:
            evidence_type, target = GENETIC_EVIDENCE[key], genetic
        elif key in EXPERIMENTAL_EVIDENCE:
            evidence_type, target = EXPERIMENTAL_EVIDENCE[key], experimental
        else:
            sys.stderr.write(f"WARNING: tab '{tab}': unknown evidence '{category} / {sheet_type}' on row {idx + 1}, skipping\n")
            continue

        target.append({
            'Evidence type': evidence_type,
            'Score': parse_score(score_cell, tab, evidence_type),
            'Citation': citation,
            'Evidence detail': detail,
        })

    return {
        'Disease_ID': record['Disease_ID'],
        'Gene': record['Gene'],
        'Locus_ID': record['Locus_ID'],
        'Inheritance': record.get('Inheritance'),
        'Curator': record.get('Curator'),
        'Contributors': record.get('Contributors'),
        'Date': record.get('Date'),
        'Description': record.get('Description'),
        'Source': source,
        'SOP_version': sop_version,
        'Manual_evidence_level': record.get('Manual_evidence_level'),
        'genetic_evidence_details': genetic,
        'experimental_evidence_details': experimental,
    }

def main(args):
    xlsx = sheet_to_xlsx_url(args.xlsx)
    sys.stderr.write(f'Reading {xlsx}\n')
    sheets = pd.read_excel(xlsx, sheet_name=None, header=None, dtype=object)

    tabs = args.tabs or [t for t in sheets if 'template' not in t.lower()]
    missing = [t for t in tabs if t not in sheets]
    if missing:
        sys.exit(f"ERROR: tabs not found in spreadsheet: {', '.join(missing)}")

    curations = []
    for tab in tabs:
        record = parse_tab(sheets[tab], tab, args.source, args.sop_version)
        if record is not None:
            curations.append(record)
            sys.stderr.write(f"Imported '{tab}' as {record['Locus_ID']} "
                             f"({len(record['genetic_evidence_details'])} genetic, "
                             f"{len(record['experimental_evidence_details'])} experimental entries)\n")

    locus_ids = [c['Locus_ID'] for c in curations]
    duplicates = sorted({l for l in locus_ids if locus_ids.count(l) > 1})
    if duplicates:
        sys.exit(f"ERROR: duplicate Locus_ID across tabs: {', '.join(duplicates)}")

    if args.merge:
        with open(args.merge) as f:
            existing = json.load(f)
        imported = set(locus_ids)
        replaced = [c['Locus_ID'] for c in existing if c['Locus_ID'] in imported]
        if replaced:
            sys.stderr.write(f"Replacing existing records: {', '.join(replaced)}\n")
        curations = [c for c in existing if c['Locus_ID'] not in imported] + curations
        curations = sorted(curations, key=lambda x: (x['Gene'], x['Locus_ID']))

    out = json.dumps(curations, indent=2, ensure_ascii=False) + '\n'
    if args.out == '-':
        sys.stdout.write(out)
    else:
        with open(args.out, 'w') as f:
            f.write(out)
        sys.stderr.write(f'Wrote {len(curations)} curations to {args.out}\n')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main(parse_args())
