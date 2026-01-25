# Data Source Ontology Context for SPARQL Generation

Use this context to generate SPARQL queries for the Gene data source knowledge graph.

## Prefixes (Always Include)

```sparql
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX ds: <https://gene.ai/ontology/datasource#>
PREFIX tax: <https://gene.ai/taxonomy/>
PREFIX inst: <https://gene.ai/data/source/>
```

## Core Classes

| Class | Description | Key Properties |
|-------|-------------|----------------|
| `ds:DataSource` | Individual data source/database | id, title, description, tier, status, website |
| `ds:Category` | Top-level domain (9 total) | prefLabel, notation |
| `ds:Subcategory` | Second-level classification (44 total) | prefLabel, notation, broader |
| `ds:Version` | Release version of a source | versionNumber, releaseDate, size, recordCount |
| `ds:License` | Licensing terms | licenseName, licenseType, commercialUseAllowed |
| `ds:AccessMethod` | How to access data | accessMethodType, apiEndpoint, authenticationRequired |
| `ds:Format` | Data format | formatType |
| `ds:CrossReference` | Link between sources | referencesSource, coveragePercentage |

## Key Object Properties

| Property | Domain | Range | Description |
|----------|--------|-------|-------------|
| `ds:belongsToCategory` | DataSource/Subcategory | Category | Parent category |
| `ds:belongsToSubcategory` | DataSource | Subcategory | Parent subcategory |
| `ds:tier` | DataSource | skos:Concept | Tier classification (1, 2, 3) |
| `ds:status` | HierarchyLevel | skos:Concept | Documentation status |
| `ds:hasLicense` | DataSource | License | License terms |
| `ds:hasAccessMethod` | DataSource | AccessMethod | Access methods |
| `ds:hasFormat` | DataSource/AccessMethod | Format | Data formats |
| `ds:hasVersion` | DataSource | Version | Available versions |
| `ds:currentVersion` | DataSource | Version | Current/latest version |
| `ds:relatedTo` | DataSource | DataSource | Related sources (symmetric) |
| `ds:complementedBy` | DataSource | DataSource | Complementary sources |
| `ds:hasCrossReference` | DataSource | CrossReference | Cross-references |
| `ds:licenseType` | License | skos:Concept | License type classification |
| `ds:accessMethodType` | AccessMethod | skos:Concept | Access method type |

## Key Datatype Properties

| Property | Domain | Range | Description |
|----------|--------|-------|-------------|
| `ds:id` | HierarchyLevel | xsd:string | Unique identifier |
| `ds:title` | - | xsd:string | Human-readable name |
| `ds:description` | - | xsd:string | Detailed description |
| `ds:website` | DataSource | xsd:anyURI | Official website |
| `ds:tag` | - | xsd:string | Keyword tags |
| `ds:commercialUseAllowed` | License | xsd:boolean | Commercial use permitted |
| `ds:attributionRequired` | License | xsd:boolean | Attribution required |
| `ds:authenticationRequired` | AccessMethod | xsd:boolean | Auth required |
| `ds:apiEndpoint` | AccessMethod | xsd:anyURI | API base URL |
| `ds:recordCount` | Version | xsd:long | Number of records |
| `ds:sizeBytes` | Version | xsd:long | Size in bytes |

## Taxonomy Concepts (for Filtering)

### Tiers (tax:TierScheme)
- `tax:Tier1` - Primary authoritative sources
- `tax:Tier2` - Secondary supplementary sources
- `tax:Tier3` - Specialized niche sources

### Status (tax:StatusScheme)
- `tax:Draft` - Documentation incomplete
- `tax:Active` - Recommended for use
- `tax:Deprecated` - No longer recommended

### Access Methods (tax:AccessMethodScheme)
- `tax:RestApi` - REST API
- `tax:SparqlEndpoint` - SPARQL endpoint
- `tax:GraphqlApi` - GraphQL API
- `tax:FtpDownload` - FTP bulk download
- `tax:S3Download` - AWS S3 download
- `tax:DirectDownload` - HTTP direct download

### License Types (tax:LicenseTypeScheme)
- `tax:PublicDomain` - No restrictions (CC0)
- `tax:OpenAccess` - Open with attribution
- `tax:CreativeCommons` - CC licensed (BY, NC, SA variants)
- `tax:AcademicOnly` - Academic/research only
- `tax:CommercialRestricted` - Requires commercial license
- `tax:Proprietary` - Subscription required

### Categories (tax:DataSourceTaxonomy)
- `tax:GeneticsGenomics` (01) - Genetics and Genomics
- `tax:CompoundsMolecules` (02) - Compounds and Molecules
- `tax:DiseasesPhenotypes` (03) - Diseases and Phenotypes
- `tax:PathwaysNetworks` (04) - Pathways and Networks
- `tax:TraditionalMedicine` (05) - Traditional Medicine
- `tax:NutritionFood` (06) - Nutrition and Food
- `tax:ProteinsMolecularBiology` (07) - Proteins and Molecular Biology
- `tax:LiteratureKnowledge` (08) - Literature and Knowledge
- `tax:Microbiome` (09) - Microbiome

### Selected Subcategories
- `tax:TraditionalChineseMedicine` (5.1) - TCM databases
- `tax:VariantRepositories` (1.1) - Genetic variant databases
- `tax:Pharmaceuticals` (2.2) - Drug databases
- `tax:CancerGenomics` (1.6) - Cancer-specific databases
- `tax:ProteinSequencesAnnotations` (7.1) - Protein databases

## Query Patterns

### List sources in a category
```sparql
SELECT ?source ?title ?tier WHERE {
  ?source a ds:DataSource ;
          ds:title ?title ;
          ds:belongsToSubcategory/ds:belongsToCategory ?cat ;
          ds:tier ?tierConcept .
  ?tierConcept skos:notation ?tier .
  ?cat skos:prefLabel "Traditional Medicine"@en .
}
```

### Find by access method
```sparql
SELECT ?source ?title WHERE {
  ?source a ds:DataSource ;
          ds:title ?title ;
          ds:hasAccessMethod ?access .
  ?access ds:accessMethodType tax:RestApi .
}
```

### Get source details
```sparql
SELECT ?title ?desc ?website ?license ?tier WHERE {
  ?source ds:id "drugbank" ;
          ds:title ?title ;
          ds:description ?desc ;
          ds:website ?website ;
          ds:tier/skos:notation ?tier .
  OPTIONAL { ?source ds:hasLicense/ds:licenseName ?license }
}
```

### Filter by license type
```sparql
SELECT ?source ?title WHERE {
  ?source a ds:DataSource ;
          ds:title ?title ;
          ds:hasLicense ?lic .
  ?lic ds:commercialUseAllowed true .
}
```

### Count by category
```sparql
SELECT ?catLabel (COUNT(?source) as ?count) WHERE {
  ?source a ds:DataSource ;
          ds:belongsToSubcategory/ds:belongsToCategory ?cat .
  ?cat skos:prefLabel ?catLabel .
}
GROUP BY ?catLabel
ORDER BY DESC(?count)
```

### Find related sources
```sparql
SELECT ?related ?title WHERE {
  ?source ds:id "uniprot" ;
          (ds:relatedTo|ds:complementedBy) ?related .
  ?related ds:title ?title .
}
```

## SHACL Validation Rules (Constraints)

- DataSource MUST have: id, title, belongsToSubcategory, tier, status
- DataSource SHOULD have: description, website, hasLicense, hasAccessMethod
- Category MUST have: prefLabel, notation (2 digits)
- Subcategory MUST have: prefLabel, notation (X.Y format), broader
- License MUST have: licenseName, licenseType

## Query Generation Guidelines

1. **Always use proper prefixes** - Include all relevant prefixes at the start
2. **Use OPTIONAL for nullable fields** - description, license details, etc.
3. **Filter with taxonomy concepts** - Use tax: URIs for tiers, status, access methods
4. **Navigate hierarchy** - Use property paths for category/subcategory traversal
5. **Return human-readable values** - Select labels, not just URIs
6. **Order results sensibly** - By tier, alphabetically, or by relevance
7. **Limit results** - Use LIMIT for large result sets
8. **Handle text search** - Use FILTER with regex or CONTAINS for keyword search
