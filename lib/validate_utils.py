import csv
import os
from typing import Dict


INVALID_ENTRIES_REPORT_FILENAME = 'invalid_entries.csv'


def read_report_from_dir(REPORT_DIR: str):
    issues: Dict[str, int]

    with open(os.path.join(REPORT_DIR, INVALID_ENTRIES_REPORT_FILENAME), 'r') as f:
        reader = csv.reader(f)
        issues = {c[0]: int(c[1]) for c in zip(*reader)} # type: ignore # pylance mistaken thinking `c` is a Tuple of 1 el
    return issues


def write_report_to_dir(issues: Dict[str, int], REPORT_DIR: str):
    with open(os.path.join(REPORT_DIR, INVALID_ENTRIES_REPORT_FILENAME), 'w') as f:
        w = csv.DictWriter(f, issues.keys())
        w.writeheader()
        w.writerow(issues)
