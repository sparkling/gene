# COCONUT Database Schema

**Document ID:** COCONUT-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** COCONUT January 2026 Release

---

## TL;DR

COCONUT (COlleCtion of Open NatUral producTs) contains 716,697 unique natural product structures from 70 collections with 70,896 unique organisms and 124,690 literature citations. The database uses a PostgreSQL backend with a Laravel REST API. Core entities are Molecules, Collections, Organisms, and Citations with full stereochemistry preservation for 42.6% of compounds.

---

## Database Statistics (January 2026)

| Entity | Count |
|--------|-------|
| **Unique Molecules** | 716,697 |
| **Total Collections** | 70 |
| **Unique Organisms** | 70,896 |
| **Organisms with IRI** | 63,924 (90.2%) |
| **Citations Mapped** | 124,690 |
| **Molecules with Organisms** | 305,138 (42.6%) |
| **Molecules with Citations** | 271,049 (37.8%) |
| **Distinct Geo Locations** | 2,678 |
| **Molecules with Geo Locations** | 63,473 (8.9%) |
| **Revoked Molecules** | 38,771 (5.1%) |

### NP-Likeness Distribution

- **Score Range:** -3.53 to 4.12
- **Mean NP-Likeness:** ~1.12
- **Higher scores indicate more "natural product-like" characteristics**

---

## Available Downloads (January 2026)

| File | Size | Description |
|------|------|-------------|
| coconut_sdf_2d_lite-01-2026.zip | 287.6 MB | 2D structures (lite) |
| coconut_sdf_2d-01-2026.zip | 691.7 MB | 2D structures (full metadata) |
| coconut_sdf_3d-01-2026.zip | 305.3 MB | 3D coordinates (not minimized) |
| coconut_csv_lite-01-2026.zip | 191 MB | CSV format (lite) |
| coconut_csv-01-2026.zip | 207.9 MB | CSV format (full) |
| coconut-dump-01-2026.sql | 31.91 GB | Complete PostgreSQL dump |

**Specialized Files (August 2024):**
- Fragment Analysis CSVs (functional groups, scaffolds)
- Drug Discovery TSV (SMILES, NP-likeness, SA scores, QED)
- Mass Spectrometry CSV (CompoundDiscoverer format)

**Download URL:** https://coconut.naturalproducts.net/download

**License:** CC0 (Creative Commons Zero) - completely unrestricted

---

## Database Architecture

### Technology Stack

- **Database:** PostgreSQL
- **Backend:** Laravel (PHP)
- **API:** RESTful with Laravel Sanctum authentication
- **Search:** Full-text + structure search
- **Hosting:** Figshare for static data, cloud for live instance

---

## Core Entity Relationship

```
MOLECULES (central)
    |
    +---> PROPERTIES (1:1)
    |
    +---> COLLECTIONS (N:M) via pivot
    |
    +---> ORGANISMS (N:M) via pivot
    |         |
    |         +---> GEO_LOCATIONS
    |         |
    |         +---> SAMPLE_LOCATIONS
    |
    +---> CITATIONS (N:M) via pivot
    |
    +---> LICENSES (via collections)
```

---

## Schema Tables

### MOLECULES

**Description:** Central table for natural product structures

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| identifier | VARCHAR | COCONUT ID (CNP0000001 format) |
| canonical_smiles | TEXT | Canonical SMILES notation |
| standard_inchi | TEXT | Standard InChI string |
| standard_inchi_key | VARCHAR(27) | InChI key |
| name | VARCHAR | Compound name (if known) |
| cas | VARCHAR | CAS registry number |
| molecular_formula | VARCHAR | Molecular formula |
| molecular_weight | DECIMAL | Molecular weight |
| exact_mass | DECIMAL | Monoisotopic mass |
| is_active | BOOLEAN | Active/revoked status |
| has_stereo | BOOLEAN | Has stereochemistry |
| stereo_defined | BOOLEAN | Stereochemistry fully defined |
| created_at | TIMESTAMP | Record creation date |
| updated_at | TIMESTAMP | Last modification date |

**Indexes:** PRIMARY KEY (id), UNIQUE (identifier), INDEX (standard_inchi_key)

---

### PROPERTIES

**Description:** Computed molecular properties (1:1 with molecules)

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| molecule_id | INTEGER | Foreign key to MOLECULES |
| np_likeness_score | DECIMAL | Natural product likeness (-4 to +4) |
| alogp | DECIMAL | Calculated LogP |
| topological_polar_surface_area | DECIMAL | TPSA |
| hydrogen_bond_acceptors | INTEGER | HBA count |
| hydrogen_bond_donors | INTEGER | HBD count |
| rotatable_bond_count | INTEGER | Rotatable bonds |
| heavy_atom_count | INTEGER | Heavy atom count |
| ring_count | INTEGER | Number of rings |
| aromatic_ring_count | INTEGER | Aromatic rings |
| fraction_csp3 | DECIMAL | Fraction of sp3 carbons |
| lipinski_rule_of_five_violations | INTEGER | Ro5 violations |
| qed_score | DECIMAL | QED drug-likeness |
| synthetic_accessibility_score | DECIMAL | SA score (1-10, lower = easier) |
| murcko_framework | VARCHAR | Murcko scaffold SMILES |

---

### COLLECTIONS

**Description:** Source databases and datasets

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| title | VARCHAR | Collection name |
| description | TEXT | Description |
| identifier | VARCHAR | Unique identifier |
| url | VARCHAR | Source URL |
| doi | VARCHAR | Collection DOI |
| citation_text | TEXT | How to cite |
| license_id | INTEGER | Foreign key to LICENSES |
| molecule_count | INTEGER | Number of molecules |
| created_at | TIMESTAMP | Creation date |
| updated_at | TIMESTAMP | Last update |

**Major Collections:**

| Collection | Compound Count |
|------------|----------------|
| FooDB | 61,899 |
| InterBioScreen Ltd | 58,550 |
| CMAUP | 43,358 |
| KNApSaCK | ~40,000 |
| ChEBI NPs | 14,357 |
| GNPS | 6,170 |

---

### COLLECTION_MOLECULE (Pivot)

**Description:** Many-to-many relationship between collections and molecules

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| collection_id | INTEGER | Foreign key to COLLECTIONS |
| molecule_id | INTEGER | Foreign key to MOLECULES |
| external_id | VARCHAR | ID in source collection |
| created_at | TIMESTAMP | Record creation |

---

### ORGANISMS

**Description:** Source organisms (biological origin)

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| name | VARCHAR | Scientific name |
| iri | VARCHAR | IRI (Internationalized Resource Identifier) |
| rank | VARCHAR | Taxonomic rank |
| molecule_count | INTEGER | Associated molecule count |
| ncbi_id | INTEGER | NCBI Taxonomy ID |
| ott_id | INTEGER | Open Tree of Life ID |
| gbif_id | INTEGER | GBIF Taxonomic ID |

**Note:** 90.2% of organisms have an IRI (63,924 of 70,896)

---

### ORGANISM_MOLECULE (Pivot)

**Description:** Many-to-many relationship between organisms and molecules

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| organism_id | INTEGER | Foreign key to ORGANISMS |
| molecule_id | INTEGER | Foreign key to MOLECULES |
| reference_doi | VARCHAR | Literature reference |
| created_at | TIMESTAMP | Record creation |

---

### SAMPLE_LOCATIONS

**Description:** Where in the organism the compound is found

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| organism_molecule_id | INTEGER | FK to ORGANISM_MOLECULE |
| name | VARCHAR | Sample location (leaf, root, etc.) |

---

### GEO_LOCATIONS

**Description:** Geographic locations where organisms are found

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| name | VARCHAR | Location name |
| latitude | DECIMAL | Latitude coordinate |
| longitude | DECIMAL | Longitude coordinate |
| country | VARCHAR | Country name |
| region | VARCHAR | Region/state |

**Statistics:** 2,678 distinct locations, 63,473 molecules with geo data

---

### CITATIONS

**Description:** Literature references

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| doi | VARCHAR | Digital Object Identifier |
| title | VARCHAR | Publication title |
| authors | TEXT | Author list |
| citation_text | TEXT | Full citation |
| year | INTEGER | Publication year |
| journal | VARCHAR | Journal name |
| pubmed_id | INTEGER | PubMed ID |

---

### CITATION_MOLECULE (Pivot)

**Description:** Many-to-many relationship between citations and molecules

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| citation_id | INTEGER | Foreign key to CITATIONS |
| molecule_id | INTEGER | Foreign key to MOLECULES |

---

### LICENSES

**Description:** Data licenses for collections

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| name | VARCHAR | License name |
| url | VARCHAR | License URL |
| spdx_id | VARCHAR | SPDX identifier |
| allows_commercial | BOOLEAN | Commercial use allowed |

---

### REPORTS

**Description:** User-submitted curation reports

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| user_id | INTEGER | Submitting user |
| molecule_id | INTEGER | Related molecule |
| title | VARCHAR | Report title |
| canonical_smiles | TEXT | Reported structure |
| doi | VARCHAR | Supporting reference |
| status | VARCHAR | pending/approved/rejected |
| created_at | TIMESTAMP | Submission date |
| updated_at | TIMESTAMP | Last update |

---

### USERS

**Description:** User accounts for authenticated access

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key |
| first_name | VARCHAR | First name |
| last_name | VARCHAR | Last name |
| username | VARCHAR | Unique username |
| email | VARCHAR | Email address |
| affiliation | VARCHAR | Institution |
| email_verified_at | TIMESTAMP | Email verification |
| created_at | TIMESTAMP | Account creation |

---

## REST API Specification

### Base URL
```
https://coconut.naturalproducts.net/api
```

### Authentication

- **Login:** `POST /api/auth/login`
  - Body: `{ "email": "...", "password": "..." }`
  - Response: `{ "access_token": "...", "token_type": "Bearer" }`

- **Register:** `POST /api/auth/register`
  - Body: `{ "first_name", "last_name", "username", "affiliation", "email", "password", "password_confirmation" }`

- **Logout:** `GET /api/auth/logout` (requires Bearer token)

### Molecule Endpoints

#### Get Molecules
```
GET /api/molecules
```

Query parameters:
- `page`: Page number
- `limit`: Results per page (default 15)
- `sort`: Sort field
- `filter`: Filter criteria

#### Search Molecules
```
POST /api/molecules/search
```

Request body:
```json
{
  "query": {
    "filters": [
      {
        "field": "properties.np_likeness_score",
        "operator": ">=",
        "value": 1.5
      },
      {
        "field": "molecular_weight",
        "operator": "between",
        "value": [200, 500]
      }
    ],
    "scopes": ["active"],
    "includes": ["properties", "organisms", "collections"]
  },
  "page": 1,
  "limit": 50
}
```

#### Mutate Molecules (Restricted)
```
POST /api/molecules/mutate
```
Requires authentication and appropriate permissions.

### Collection Endpoints

```
GET /api/collections
POST /api/collections/search
POST /api/collections/mutate (restricted)
```

### Organism Endpoints

```
GET /api/organisms
POST /api/organisms/search
POST /api/organisms/mutate (restricted)
```

Search fields: `name`, `iri`, `rank`, `molecule_count`

### Citation Endpoints

```
GET /api/citations
POST /api/citations/search
POST /api/citations/mutate (restricted)
```

Search fields: `doi`, `title`, `authors`, `citation_text`

### Properties Endpoints

```
GET /api/properties
POST /api/properties/search
POST /api/properties/mutate (restricted)
```

### Reports Endpoints

```
POST /api/reports/search
POST /api/reports
PATCH /api/reports/{id} (restricted)
DELETE /api/reports (restricted)
```

### Special Endpoints

#### BioSchemas (Schema.org)
```
GET /api/schemas/bioschemas/{id}
```
Returns Schema.org MolecularEntity JSON-LD.

---

## Sample API Response

### Molecule Object

```json
{
  "id": 12345,
  "identifier": "CNP0012345",
  "name": "Quercetin",
  "canonical_smiles": "Oc1cc(O)c2c(=O)c(O)c(-c3ccc(O)c(O)c3)oc2c1",
  "standard_inchi": "InChI=1S/C15H10O7/c16-7-4-10(19)12-11(5-7)22-15(14(21)13(12)20)6-1-2-8(17)9(18)3-6/h1-5,16-19,21H",
  "standard_inchi_key": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
  "molecular_formula": "C15H10O7",
  "molecular_weight": 302.236,
  "exact_mass": 302.0427,
  "is_active": true,
  "has_stereo": false,
  "stereo_defined": true,

  "properties": {
    "np_likeness_score": 2.34,
    "alogp": 1.54,
    "topological_polar_surface_area": 131.36,
    "hydrogen_bond_acceptors": 7,
    "hydrogen_bond_donors": 5,
    "rotatable_bond_count": 1,
    "heavy_atom_count": 22,
    "ring_count": 3,
    "aromatic_ring_count": 3,
    "lipinski_rule_of_five_violations": 0,
    "qed_score": 0.43,
    "synthetic_accessibility_score": 3.2
  },

  "collections": [
    {
      "id": 5,
      "title": "FooDB",
      "external_id": "FDB012345"
    },
    {
      "id": 12,
      "title": "KNApSaCK",
      "external_id": "C00001234"
    }
  ],

  "organisms": [
    {
      "id": 456,
      "name": "Allium cepa",
      "rank": "species",
      "iri": "http://www.wikidata.org/entity/Q23195",
      "ncbi_id": 4679,
      "sample_locations": ["bulb", "leaf"],
      "geo_locations": [
        {"name": "Central Asia", "country": null}
      ]
    }
  ],

  "citations": [
    {
      "id": 789,
      "doi": "10.1016/j.foodchem.2020.127123",
      "title": "Quercetin in food sources...",
      "year": 2020
    }
  ]
}
```

---

## CSV File Structure

### Lite Format

| Column | Description |
|--------|-------------|
| coconut_id | COCONUT identifier |
| canonical_smiles | SMILES notation |
| inchikey | InChI key |
| molecular_formula | Formula |
| molecular_weight | MW |

### Full Format (Additional Columns)

| Column | Description |
|--------|-------------|
| name | Compound name |
| cas | CAS number |
| np_likeness_score | NP-likeness |
| alogp | LogP |
| tpsa | Topological PSA |
| hba | H-bond acceptors |
| hbd | H-bond donors |
| rotatable_bonds | Rotatable bonds |
| qed_score | QED |
| sa_score | Synthetic accessibility |
| collections | Pipe-separated list |
| organisms | Pipe-separated list |

---

## SDF File Structure

### 2D SDF Properties

```
> <coconut_id>
CNP0012345

> <name>
Quercetin

> <inchi>
InChI=1S/C15H10O7/c16-7-4-10(19)...

> <inchikey>
REFJWTPEDVJJIY-UHFFFAOYSA-N

> <molecular_formula>
C15H10O7

> <molecular_weight>
302.236

> <np_likeness_score>
2.34

> <collections>
FooDB|KNApSaCK|CMAUP

> <organisms>
Allium cepa|Brassica oleracea
```

### 3D SDF

Same properties as 2D but with 3D coordinates generated by RDKit (not energy-minimized).

---

## Integration Notes for Knowledge Base

### Key Identifiers for Cross-Referencing

| Identifier | Format | Example | Cross-ref to |
|------------|--------|---------|--------------|
| COCONUT ID | CNP + 7 digits | CNP0012345 | Internal |
| InChI Key | 27 characters | REFJWTPEDVJJIY-UHFFFAOYSA-N | PubChem, ChEMBL, LOTUS |
| CAS Number | Variable | 117-39-5 | SciFinder, commercial |
| Collection External ID | Varies | FDB012345 | Source database |
| NCBI Taxonomy ID | Integer | 4679 | NCBI, UniProt |
| Wikidata IRI | URL | http://www.wikidata.org/entity/Q23195 | Wikidata, LOTUS |
| DOI | doi:... | 10.1016/j.foodchem.2020.127123 | CrossRef, PubMed |

### Data Quality Indicators

- **is_active:** True = current, False = revoked/deprecated
- **stereo_defined:** True = complete stereochemistry
- **Organism IRI coverage:** 90.2% have validated IRIs
- **Citation coverage:** 37.8% of molecules have literature references

### Licensing

**License:** CC0 (Creative Commons Zero)
**Attribution Required:** No (but appreciated)
**Commercial Use:** Unrestricted

---

## Query Examples

### PostgreSQL Queries (SQL Dump)

#### Find Highly NP-Like Compounds

```sql
SELECT m.identifier, m.name, m.canonical_smiles, p.np_likeness_score
FROM molecules m
JOIN properties p ON m.id = p.molecule_id
WHERE p.np_likeness_score > 2.0
AND m.is_active = true
ORDER BY p.np_likeness_score DESC
LIMIT 100;
```

#### Find Compounds from Specific Organism

```sql
SELECT m.identifier, m.name, o.name AS organism
FROM molecules m
JOIN organism_molecule om ON m.id = om.molecule_id
JOIN organisms o ON om.organism_id = o.id
WHERE o.name LIKE 'Curcuma%'
AND m.is_active = true;
```

#### Find Compounds with Drug-Like Properties

```sql
SELECT m.identifier, m.name, p.*
FROM molecules m
JOIN properties p ON m.id = p.molecule_id
WHERE p.lipinski_rule_of_five_violations = 0
AND p.qed_score > 0.5
AND m.is_active = true
LIMIT 100;
```

### API Query Examples

#### Search by NP-Likeness (cURL)

```bash
curl -X POST https://coconut.naturalproducts.net/api/molecules/search \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "filters": [
        {"field": "properties.np_likeness_score", "operator": ">=", "value": 2.0}
      ],
      "includes": ["properties", "organisms"]
    },
    "limit": 50
  }'
```

#### Search by Organism (cURL)

```bash
curl -X POST https://coconut.naturalproducts.net/api/organisms/search \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "filters": [
        {"field": "name", "operator": "like", "value": "%Curcuma%"}
      ],
      "includes": ["molecules"]
    },
    "limit": 10
  }'
```

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation |

---

## References

1. Venkata C, et al. (2024) "COCONUT 2.0: a comprehensive overhaul and curation of the collection of open natural products database." Nucleic Acids Res. gkae1063.

2. COCONUT Documentation: https://steinbeck-lab.github.io/coconut/

3. API Documentation: https://coconut.naturalproducts.net/api-documentation

4. GitHub Repository: https://github.com/Steinbeck-Lab/coconut
