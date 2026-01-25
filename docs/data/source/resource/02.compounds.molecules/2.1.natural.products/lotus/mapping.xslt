<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  LOTUS -> Natural Products Unified Schema Mapping
  ============================================================
  Source: ./schema.md (LOTUS - Linked Open Total Unified Structure-organism)
  Target: ../schema.json (2.1 Natural Products Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ wikidataId                  │ wikidata_id         │ Direct       │
  │ name                        │ name                │ Direct       │
  │ canonicalSmiles             │ canonical_smiles    │ Direct       │
  │ isomericSmiles              │ isomeric_smiles     │ Direct       │
  │ inchi                       │ inchi               │ Direct       │
  │ inchiKey (P235)             │ inchi_key           │ Direct       │
  │ molecularFormula            │ molecular_formula   │ Direct       │
  │ molecularWeight             │ molecular_weight    │ Direct       │
  │ taxon (P703)                │ organism_sources    │ Object array │
  │ ncbiTaxonId (P685)          │ ncbi_taxon_id       │ Direct       │
  │ nplScore                    │ npl_score           │ Direct       │
  │ deepSmiles                  │ deep_smiles         │ Direct       │
  │ reference (P248)            │ references          │ Object array │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" or "null" → null
  - Wikidata properties use Q-identifiers

  NOTES:
  - LOTUS is integrated with Wikidata using SPARQL
  - Wikidata property IDs: P235 (InChIKey), P703 (found in taxon), P248 (stated in)
  - License: CC0 (Creative Commons Zero)
-->
<xsl:stylesheet version="3.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:fn="http://www.w3.org/2005/xpath-functions"
    xmlns:map="http://www.w3.org/2005/xpath-functions/map"
    xmlns:array="http://www.w3.org/2005/xpath-functions/array"
    xmlns:local="http://local.functions"
    exclude-result-prefixes="xs fn map array local">

  <xsl:output method="json" indent="yes"/>

  <!-- ============================================================
       Entry Point: Parse JSON input
       ============================================================ -->
  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <!-- ============================================================
       Main Record Transformation
       ============================================================ -->
  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== IDENTIFIERS ========== -->

      <!-- Wikidata ID (Q-identifier) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='wikidataId'])">
          <fn:null key="wikidata_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="wikidata_id">
            <xsl:value-of select="fn:string[@key='wikidataId']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <!-- Compound name -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <!-- Canonical SMILES (P233) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='canonicalSmiles'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='canonicalSmiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Isomeric SMILES (P2017) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='isomericSmiles'])">
          <fn:null key="isomeric_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="isomeric_smiles">
            <xsl:value-of select="fn:string[@key='isomericSmiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- InChI (P234) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchi'])">
          <fn:null key="inchi"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi">
            <xsl:value-of select="fn:string[@key='inchi']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- InChI Key (P235) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchiKey'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='inchiKey']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- DeepSMILES (LOTUS-specific) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='deepSmiles'])">
          <fn:null key="deep_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="deep_smiles">
            <xsl:value-of select="fn:string[@key='deepSmiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MOLECULAR PROPERTIES ========== -->

      <!-- Molecular formula -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='molecularFormula'])">
          <fn:null key="molecular_formula"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_formula">
            <xsl:value-of select="fn:string[@key='molecularFormula']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Molecular weight -->
      <xsl:choose>
        <xsl:when test="not(fn:number[@key='molecularWeight']) or fn:number[@key='molecularWeight'] &lt;= 0">
          <fn:null key="molecular_weight"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='molecularWeight']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Natural product likeness score -->
      <xsl:choose>
        <xsl:when test="not(fn:number[@key='nplScore'])">
          <fn:null key="npl_score"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="npl_score">
            <xsl:value-of select="fn:number[@key='nplScore']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ORGANISM SOURCES (P703: found in taxon) ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='taxa']/fn:map">
          <fn:array key="organism_sources">
            <xsl:for-each select="fn:array[@key='taxa']/fn:map">
              <fn:map>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='taxonName']"/>
                </fn:string>
                <xsl:if test="fn:number[@key='ncbiTaxonId']">
                  <fn:number key="ncbi_taxon_id">
                    <xsl:value-of select="fn:number[@key='ncbiTaxonId']"/>
                  </fn:number>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="organism_sources"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Primary NCBI taxon ID -->
      <xsl:choose>
        <xsl:when test="fn:number[@key='ncbiTaxonId']">
          <fn:number key="ncbi_taxon_id">
            <xsl:value-of select="fn:number[@key='ncbiTaxonId']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="ncbi_taxon_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== REFERENCES (P248: stated in) ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='references']/fn:map">
          <fn:array key="references">
            <xsl:for-each select="fn:array[@key='references']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='doi']">
                  <fn:string key="doi">
                    <xsl:value-of select="fn:string[@key='doi']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:number[@key='pmid']">
                  <fn:number key="pmid">
                    <xsl:value-of select="fn:number[@key='pmid']"/>
                  </fn:number>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="references"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">LOTUS</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='wikidataId']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Utility Functions
       ============================================================ -->

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
