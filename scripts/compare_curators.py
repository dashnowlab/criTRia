# Compare criTRia curations from two curators (e.g. an individual vs the group) and write one
# row per locus scored by both, for inter-rater agreement (ICC on scores, weighted kappa on
# classifications) and citation overlap.
#
# Inputs are criTRia-curations JSON files as written by import_curations.py. Scores and
# classifications are calculated here with check-curations.py's summarize_curations, so the
# inputs do not need to have been through check-curations.py first.
#
# A locus is included only if BOTH curators have completed it, judged from the data:
# at least one evidence entry with a numeric Score. Entries left unscored simply don't count
# towards the totals. Loci with no scored evidence (blank, or classification-only via
# Manual_evidence_level) are excluded and listed on stderr, so they join automatically once scored.
# Entries with a source or details but no score are likely unfinished: they are warned about on
# stderr and listed in the <curator>_incomplete_entries column.
#
# Example:
#   python scripts/compare_curators.py \
#       --rater gabriel data/curators/curations-gabriel.json \
#       --rater group data/curators/STRchive-criTRia-curations.json \
#       -o data/curators/curator-comparison.tsv
import sys
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import pandas as pd

# check-curations.py has a hyphen in its name, so load it by path
_spec = importlib.util.spec_from_file_location('check_curations', Path(__file__).parent / 'check-curations.py')
check_curations = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_curations)

def parse_args():
    parser = argparse.ArgumentParser(description='Compare curations from two curators, one row per locus completed by both.')
    parser.add_argument('--rater', nargs=2, action='append', metavar=('NAME', 'JSON'), required=True,
                        help='Curator name (used as column prefix) and curations JSON file. Give exactly twice.')
    parser.add_argument('-o', '--out', default='-',
                        help='Output TSV file path (default: stdout)')
    args = parser.parse_args()
    if len(args.rater) != 2:
        parser.error('--rater must be given exactly twice')
    return args

def evidence_entries(locus):
    return locus.get('genetic_evidence_details', []) + locus.get('experimental_evidence_details', [])

def curation_status(locus):
    """Return None if the locus has at least one scored evidence entry, otherwise the reason it does not.
    Unscored entries alongside scored ones are fine; they don't count towards the totals.
    >>> curation_status({'genetic_evidence_details': [{'Score': 1}], 'experimental_evidence_details': [{'Score': None}]}) is None
    True
    >>> curation_status({'genetic_evidence_details': [], 'experimental_evidence_details': []})
    'no evidence entries'
    >>> curation_status({'genetic_evidence_details': [{'Evidence type': 'Probands', 'Score': None}], 'experimental_evidence_details': []})
    'no scored evidence'
    """
    entries = evidence_entries(locus)
    if not entries:
        return 'no evidence entries'
    if not any(isinstance(e.get('Score'), (int, float)) for e in entries):
        return 'no scored evidence'
    return None

def incomplete_entries(locus):
    """Evidence types that have a Citation or Evidence detail but no Score; likely unfinished.
    >>> incomplete_entries({'genetic_evidence_details': [{'Evidence type': 'Probands', 'Score': None, 'Citation': 'pmid:1'}, {'Evidence type': 'Allele', 'Score': 1, 'Citation': 'pmid:1'}], 'experimental_evidence_details': [{'Evidence type': 'Cell culture', 'Score': None, 'Citation': None, 'Evidence detail': None}]})
    ['Probands']
    """
    return [e.get('Evidence type') for e in evidence_entries(locus)
            if not isinstance(e.get('Score'), (int, float)) and (e.get('Citation') or e.get('Evidence detail'))]

def citations(locus):
    """Unique citation IDs from the evidence Citation fields, lowercased and sorted.
    >>> citations({'genetic_evidence_details': [{'Citation': 'PMID:2; pmid:1'}], 'experimental_evidence_details': [{'Citation': 'pmid:2'}, {'Citation': None}]})
    ['pmid:1', 'pmid:2']
    """
    found = set()
    for entry in evidence_entries(locus):
        for part in (entry.get('Citation') or '').split(';'):
            part = part.strip().lower()
            if part:
                found.add(part)
    return sorted(found)

def summarize(locus):
    """Score a locus, returning both the calculated classification and the final one (after Manual_evidence_level)."""
    manual = locus.get('Manual_evidence_level')
    scored = check_curations.summarize_curations({**copy.deepcopy(locus), 'Manual_evidence_level': None})
    supercategories = scored.get('supercategory_summary', {})
    cites = citations(locus)
    return {
        'curator': locus.get('Curator'),
        'date': locus.get('Date'),
        'inheritance': locus.get('Inheritance'),
        'total_score': scored.get('total_score'),
        'genetic_score': supercategories.get('Genetic Evidence'),
        'experimental_score': supercategories.get('Experimental Evidence'),
        'classification_calculated': scored.get('classification'),
        'manual_evidence_level': manual,
        'classification_final': manual or scored.get('classification'),
        'publication_count': scored.get('publication_count'),
        'publication_interval_years': scored.get('publication_interval_years'),
        'incomplete_entries': '; '.join(incomplete_entries(locus)),
        'n_citations': len(cites),
        'citations': '; '.join(cites),
    }

def main(args):
    (name_a, file_a), (name_b, file_b) = args.rater
    with open(file_a) as f:
        loci_a = {c['Locus_ID']: c for c in json.load(f)}
    with open(file_b) as f:
        loci_b = {c['Locus_ID']: c for c in json.load(f)}

    rows = []
    for locus_id in sorted(set(loci_a) | set(loci_b)):
        a, b = loci_a.get(locus_id), loci_b.get(locus_id)
        status_a = curation_status(a) if a else 'not curated'
        status_b = curation_status(b) if b else 'not curated'
        if status_a or status_b:
            # Only report loci that both curators have started, to keep the log short
            if a and b and not (status_a == 'no evidence entries' and status_b == 'no evidence entries'):
                sys.stderr.write(f'Excluding {locus_id}: {name_a} {status_a or "complete"}; {name_b} {status_b or "complete"}\n')
            continue

        sys.stderr.write(f'Scoring {locus_id}\n')
        for name, locus in ((name_a, a), (name_b, b)):
            incomplete = incomplete_entries(locus)
            if incomplete:
                sys.stderr.write(f"WARNING: {locus_id}: {name} has source/details but no score for {', '.join(incomplete)}; likely incomplete\n")
        row = {'Locus_ID': locus_id, 'Gene': b.get('Gene'), 'Disease_ID': b.get('Disease_ID')}
        summary_a, summary_b = summarize(a), summarize(b)
        row.update({f'{name_a}_{k}': v for k, v in summary_a.items()})
        row.update({f'{name_b}_{k}': v for k, v in summary_b.items()})
        shared = sorted(set(citations(a)) & set(citations(b)))
        row['n_shared_citations'] = len(shared)
        row['shared_citations'] = '; '.join(shared)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(args.out if args.out != '-' else sys.stdout, sep='\t', index=False)
    sys.stderr.write(f'Wrote {len(df)} loci completed by both {name_a} and {name_b}' + (f' to {args.out}\n' if args.out != '-' else '\n'))

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main(parse_args())
