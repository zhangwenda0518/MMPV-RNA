#!/usr/bin/env python3
"""Consolidated unit tests for lineage_to_ranks() — 10 cases covering all known edge cases.

Usage:
    python test_lineage_ranks.py          # run all tests
    python -m pytest test_lineage_ranks.py -v  # verbose with pytest (if installed)
"""

import sys
from pathlib import Path

# Resolve import path relative to script location (tests/ -> pipeline dir)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import importlib
import virus_classifier
importlib.reload(virus_classifier)

fn = virus_classifier.lineage_to_ranks
RANK_NAMES = virus_classifier.RANK_NAMES

# ── Test cases ──
# Each case: (label, lineage_string, expected_rank_dict)
# Key: realm, kingdom, phylum, class, order, family, genus, species (lowercase from RANK_NAMES)
# Use "NA" for ranks expected to be absent

TEST_CASES = [
    # ── 1. ICTV format with subfamily ──
    ("ICTV: Zimmerviridae + Jacobvirinae subfamily",
     "Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes; Zimmerviridae; Jacobvirinae; Jouyvirus; Jouyvirus ev017",
     {"realm": "Duplodnaviria", "kingdom": "Heunggongvirae", "phylum": "Uroviricota",
      "class": "Caudoviricetes", "order": "NA", "family": "Zimmerviridae",
      "genus": "Jouyvirus", "species": "Jouyvirus ev017"}),

    # ── 2. Missing Order + Family (Class → Genus directly) ──
    ("ACVirus: Class→Genus, missing Order/Family",
     "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Jouyvirus;Jouyvirus ev017",
     {"realm": "Duplodnaviria", "kingdom": "Heunggongvirae", "phylum": "Uroviricota",
      "class": "Caudoviricetes", "order": "NA", "family": "NA",
      "genus": "Jouyvirus", "species": "Jouyvirus ev017"}),

    # ── 3. Realm→Family with 4 missing ranks + incertae sedis ──
    ("Metabuli: Realm→Family jump, 4 ranks missing, incertae sedis",
     "Viruses; Ribozyviria; Ribozyviria incertae sedis; Kolmioviridae; Deltavirus; Deltavirus italiense; Hepatitis delta virus",
     {"realm": "Ribozyviria", "kingdom": "NA", "phylum": "NA", "class": "NA",
      "order": "NA", "family": "Kolmioviridae", "genus": "Deltavirus",
      "species": "Deltavirus italiense"}),

    # ── 4. Complete 8-rank golden case ──
    ("CAT: complete 8-rank lineage",
     "Viruses; Varidnaviria; Bamfordvirae; Nucleocytoviricota; Megaviricetes; Imitervirales; Mimiviridae; Fadolivirus; Fadolivirus algeromassiliense",
     {"realm": "Varidnaviria", "kingdom": "Bamfordvirae", "phylum": "Nucleocytoviricota",
      "class": "Megaviricetes", "order": "Imitervirales", "family": "Mimiviridae",
      "genus": "Fadolivirus", "species": "Fadolivirus algeromassiliense"}),

    # ── 5. Short lineage (only to Class) ──
    ("VITAP: short, only to Class",
     "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes",
     {"realm": "Duplodnaviria", "kingdom": "Heunggongvirae", "phylum": "Uroviricota",
      "class": "Caudoviricetes", "order": "NA", "family": "NA",
      "genus": "NA", "species": "NA"}),

    # ── 6. Family name where Order should be (Peduoviridae after Class) ──
    ("Family where Order expected: Peduoviridae after Caudoviricetes",
     "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Peduoviridae;Lucasvirus;Lucasvirus L",
     {"realm": "Duplodnaviria", "kingdom": "Heunggongvirae", "phylum": "Uroviricota",
      "class": "Caudoviricetes", "order": "NA", "family": "Peduoviridae",
      "genus": "Lucasvirus", "species": "Lucasvirus L"}),

    # ── 7. NCBI format with subfamily, no species ──
    ("NCBI Zimmer: subfamily present, no species",
     "Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes; Zimmerviridae; Jacobvirinae; Jouyvirus",
     {"realm": "Duplodnaviria", "kingdom": "Heunggongvirae", "phylum": "Uroviricota",
      "class": "Caudoviricetes", "order": "NA", "family": "Zimmerviridae",
      "genus": "Jouyvirus", "species": "NA"}),

    # ── 8. Host lineage: should return all NA ──
    ("Metabuli: host lineage (all virus ranks should be NA)",
     "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Boreoeutheria; Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; Muridae; Murinae; Mus; Mus;Mus musculus",
     {"realm": "NA", "kingdom": "NA", "phylum": "NA", "class": "NA",
      "order": "NA", "family": "NA", "genus": "NA", "species": "Mus musculus"}),

    # ── 9. CAT dash placeholder (\"-\" in genus position) ──
    ("CAT: dash placeholder for missing genus",
     "Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; Picornavirales; Dicistroviridae; -; Triatovirus himetobi",
     {"realm": "Riboviria", "kingdom": "Orthornavirae", "phylum": "Pisuviricota",
      "class": "Pisoniviricetes", "order": "Picornavirales", "family": "Dicistroviridae",
      "genus": "NA", "species": "Triatovirus himetobi"}),

    # ── 10. Subgenus between Genus and Species (MSL38+) ──
    ("Subgenus: Potyvirus (genus) → Potyvirus yituberosi (species)",
     "Viruses;Riboviria;Orthornavirae;Pisuviricota;Stelpaviricetes;Patatavirales;Potyviridae;Potyvirus;Potyvirus yituberosi",
     {"realm": "Riboviria", "kingdom": "Orthornavirae", "phylum": "Pisuviricota",
      "class": "Stelpaviricetes", "order": "Patatavirales", "family": "Potyviridae",
      "genus": "Potyvirus", "species": "Potyvirus yituberosi"}),
]


def run_tests():
    """Run all test cases and report results."""
    print("=" * 72)
    print(f"Testing lineage_to_ranks() — {len(TEST_CASES)} cases")
    print("=" * 72)

    passed = 0
    failed = 0

    for label, lineage, expected in TEST_CASES:
        result = fn(lineage)
        ok = True
        mismatches = []

        for rank in RANK_NAMES:
            got = result.get(rank, "NA")
            exp = expected.get(rank, "NA")
            if got != exp:
                ok = False
                mismatches.append(f"  {rank}: expected='{exp}' got='{got}'")

        status = "PASS" if ok else "FAIL"
        if ok:
            passed += 1
        else:
            failed += 1

        print(f"\n[{status}] {label}")
        print(f"  Input: {lineage[:100]}{'...' if len(lineage) > 100 else ''}")

        result_str = " | ".join(f"{r}={result.get(r, 'NA')}" for r in RANK_NAMES)
        print(f"  Result: {result_str}")

        if mismatches:
            for m in mismatches:
                print(m)

    print(f"\n{'=' * 72}")
    print(f"Results: {passed} passed, {failed} failed out of {len(TEST_CASES)}")
    print(f"{'=' * 72}")

    return failed == 0


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
