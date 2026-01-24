---
id: ttd
title: "TTD - Therapeutic Target Database"
type: data-source
category: compounds
subcategory: compound-target-interactions
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [therapeutic-targets, drug-targets, druggability, clinical-trials, drug-development]
---

# TTD - Therapeutic Target Database

**Category:** [Compounds & Molecules](../../_index.md) > [Compound-Target Interactions](../_index.md)

## Overview

The Therapeutic Target Database (TTD) is a comprehensive resource providing information about therapeutic protein and nucleic acid targets, targeted diseases, pathway information, and corresponding drugs. TTD uniquely organizes druggability data across three perspectives: molecular interactions, human system profiles, and cell-based expression variations.

TTD classifies targets by development status (successful, clinical trial, preclinical, literature-reported) and includes extensive information on drug modalities including small molecules, antibodies, peptides, nucleic acid drugs, cell therapies, and gene therapies. This comprehensive coverage makes it valuable for understanding the therapeutic landscape.

The database integrates target information with disease associations (ICD-11 coded), pathway data, and cross-references to major databases including UniProt, DrugBank, ChEMBL, and PubChem.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Targets | 3,131 |
| Successful Targets | 426 |
| Clinical Trial Targets | 1,014 |
| Total Drugs | 39,862 |
| Approved Drugs | 2,895 |
| Clinical Trial Drugs | 11,796 |
| Last Update | 2024 |

## Primary Use Cases

1. **Target Assessment** - Evaluate druggability of therapeutic targets
2. **Drug-Target Mapping** - Link drugs to their molecular targets
3. **Clinical Landscape** - Track drug development status
4. **Disease-Target Links** - Identify targets for diseases
5. **Competitive Intelligence** - Monitor target development status

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| TTD Target ID | TTDT + 5 digits | TTDT00001 |
| TTD Drug ID | 6 characters | D0A9YA |
| UniProt ID | Accession | P15056 |
| ICD-11 Code | Disease code | 2B90.0 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://idrblab.net/ttd/ | Search and browse |
| Bulk Download | https://idrblab.net/ttd/full-data-download | TXT files |
| Target Search | Available | By name, gene, disease |

## Target Status Categories

| Status | Count | Description |
|--------|-------|-------------|
| Successful | 426 | FDA-approved drug exists |
| Clinical Trial | 1,014 | In clinical development |
| Preclinical | 212 | Preclinical/patented |
| Literature | 1,479 | Research targets |

## Drug Types

| Type | Description |
|------|-------------|
| Small Molecule | Traditional chemical drugs |
| Antibody | mAbs, ADCs, bispecifics |
| Peptide | Therapeutic peptides |
| Nucleic Acid | ASOs, siRNAs, mRNAs |
| Cell Therapy | CAR-T, stem cells |
| Gene Therapy | Gene delivery vectors |

## License

| Aspect | Value |
|--------|-------|
| License | Free for academic use |
| Commercial Use | Limited |
| Attribution Required | Yes |

## Related Resources

- [BindingDB](../bindingdb/_index.md) - Binding affinities
- [GtoPdb](../gtopdb/_index.md) - Pharmacology
- [DrugBank](../../2.2.pharmaceuticals/drugbank/_index.md) - Drug information

## See Also

- [Schema Documentation](./schema.md)

## References

1. Zhou Y, et al. (2024) "TTD: Therapeutic Target Database describing target druggability information." Nucleic Acids Res. 52(D1):D1465-D1477.
