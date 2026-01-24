# ChEMBL - LLM Context Reference

> Bioactivity database: 24M compound-target activity measurements from literature with drug-likeness annotations.

## Quick Reference

| Field | Value |
|-------|-------|
| **URL** | https://www.ebi.ac.uk/chembl |
| **Maintainer** | EMBL-EBI |
| **License** | CC BY-SA 3.0 |
| **Commercial OK** | Yes (share-alike) |
| **Update Freq** | Quarterly |
| **Version** | ChEMBL 36 (July 2025) |

## Content Summary

ChEMBL aggregates bioactivity data from medicinal chemistry literature, linking small molecules to protein targets with quantitative IC50/Ki/EC50 measurements. Each activity record traces from compound through assay to target with full provenance to source publication.

**Record Counts:**
- Activities: 24.3M
- Compounds: 2.9M distinct
- Targets: 17,803
- Publications: 99,142

## Key Identifiers

| ID Type | Format | Example |
|---------|--------|---------|
| ChEMBL ID | CHEMBL + digits | CHEMBL25 |
| InChI Key | 27-char hash | BSYNRYMUTXBXSQ-UHFFFAOYSA-N |

**Cross-references:** UniProt (targets), PubChem CID, DrugBank ID, PubMed ID

## Core Schema

### Primary Entity: MOLECULE_DICTIONARY

| Field | Type | Description |
|-------|------|-------------|
| chembl_id | VARCHAR(20) | Primary identifier |
| pref_name | VARCHAR(255) | Compound name |
| max_phase | NUMERIC | Clinical phase (0-4) |
| natural_product | SMALLINT | NP flag (0/1) |
| therapeutic_flag | SMALLINT | Is therapeutic drug |

### Key Relationships

```
MOLECULE_DICTIONARY --[molregno]--> ACTIVITIES
ACTIVITIES --[assay_id]--> ASSAYS
ASSAYS --[tid]--> TARGET_DICTIONARY
TARGET_DICTIONARY --[component_id]--> COMPONENT_SEQUENCES (UniProt)
```

## Access Methods

### API
```
https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json
```
Rate limit: None stated (be reasonable)

### Bulk Download
```
https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36_sqlite.tar.gz
```
Format: SQLite, MySQL, PostgreSQL, SDF, HDF5
Size: 5.2 GB (SQLite)

## Query Examples

### Get compound with activities
```sql
SELECT md.chembl_id, md.pref_name, a.standard_type,
       a.standard_value, a.standard_units, td.pref_name AS target
FROM molecule_dictionary md
JOIN activities a ON md.molregno = a.molregno
JOIN assays ay ON a.assay_id = ay.assay_id
JOIN target_dictionary td ON ay.tid = td.tid
WHERE md.chembl_id = 'CHEMBL25';
```

### Find potent compounds for target
```sql
SELECT md.chembl_id, a.pchembl_value, a.standard_type
FROM activities a
JOIN molecule_dictionary md ON a.molregno = md.molregno
JOIN assays ay ON a.assay_id = ay.assay_id
WHERE ay.tid = (SELECT tid FROM target_dictionary WHERE chembl_id = 'CHEMBL220')
  AND a.pchembl_value > 7
ORDER BY a.pchembl_value DESC;
```

## Sample Record

```json
{
  "molecule_chembl_id": "CHEMBL25",
  "pref_name": "ASPIRIN",
  "max_phase": 4,
  "molecule_properties": {
    "full_mwt": 180.16,
    "alogp": 1.31,
    "hba": 3,
    "hbd": 1,
    "psa": 63.60
  },
  "molecule_structures": {
    "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "standard_inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
  }
}
```

## Integration Notes

- **Primary use:** Drug-target bioactivity lookup, compound profiling
- **Best for:** Finding activity data, target selectivity analysis
- **Limitations:** Literature-derived (not standardized assays), confidence varies
- **Pairs with:** UniProt (targets), PubChem (structures), DrugBank (drug info)

---
*Size: ~5GB | Updated: January 2026*
