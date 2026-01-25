<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  COCONUT -> Natural Products Unified Schema Mapping
  ============================================================
  Source: ./schema.md (COCONUT Natural Products Database)
  Target: ../schema.json (2.1 Natural Products Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ coconut_id                  │ compound_id         │ Direct       │
  │ coconut_id                  │ coconut_id          │ Direct       │
  │ name                        │ name                │ Direct       │
  │ canonical_smiles            │ canonical_smiles    │ Direct       │
  │ isomeric_smiles             │ isomeric_smiles     │ Direct       │
  │ inchi                       │ inchi               │ Direct       │
  │ inchi_key                   │ inchi_key           │ Direct       │
  │ molecular_formula           │ molecular_formula   │ Direct       │
  │ molecular_weight            │ molecular_weight    │ Direct       │
  │ organisms                   │ organism_sources    │ Array        │
  │ sugar_free_smiles           │ sugar_free_smiles   │ Direct       │
  │ qed_drug_likeliness         │ qed_drug_likeness   │ Direct       │
  │ lipinski_rule_of_5          │ lipinski_violations │ Direct       │
  │ murko_framework             │ murko_framework     │ Direct       │
  │ sources                     │ references          │ Array        │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - name: "" → null
  - smiles fields: "" or "-" → null
  - molecular_weight: -1 or 0 → null

  NOTES:
  - COCONUT uses CNP prefix for compound IDs (e.g., CNP0123456)
  - qed_drug_likeliness maps to qed_drug_likeness (spelling normalization)
  - Organism data includes NCBI taxonomy IDs
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

      <!-- Primary compound identifier -->
      <fn:string key="compound_id">
        <xsl:value-of select="fn:string[@key='coconut_id']"/>
      </fn:string>

      <!-- COCONUT-specific identifier -->
      <fn:string key="coconut_id">
        <xsl:value-of select="fn:string[@key='coconut_id']"/>
      </fn:string>

      <!-- ========== BASIC PROPERTIES ========== -->

      <!-- Compound name with null handling -->
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

      <!-- Canonical SMILES -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='canonical_smiles'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='canonical_smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Isomeric SMILES -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='isomeric_smiles'])">
          <fn:null key="isomeric_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="isomeric_smiles">
            <xsl:value-of select="fn:string[@key='isomeric_smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- InChI -->
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

      <!-- InChI Key -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='inchi_key'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='inchi_key']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Sugar-free SMILES (COCONUT-specific) -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='sugar_free_smiles'])">
          <fn:null key="sugar_free_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="sugar_free_smiles">
            <xsl:value-of select="fn:string[@key='sugar_free_smiles']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== MOLECULAR PROPERTIES ========== -->

      <!-- Molecular formula -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='molecular_formula'])">
          <fn:null key="molecular_formula"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_formula">
            <xsl:value-of select="fn:string[@key='molecular_formula']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Molecular weight -->
      <xsl:choose>
        <xsl:when test="not(fn:number[@key='molecular_weight']) or fn:number[@key='molecular_weight'] &lt;= 0">
          <fn:null key="molecular_weight"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="molecular_weight">
            <xsl:value-of select="fn:number[@key='molecular_weight']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DRUG-LIKENESS PROPERTIES ========== -->

      <!-- QED drug-likeness (spelling normalized) -->
      <xsl:choose>
        <xsl:when test="not(fn:number[@key='qed_drug_likeliness'])">
          <fn:null key="qed_drug_likeness"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="qed_drug_likeness">
            <xsl:value-of select="fn:number[@key='qed_drug_likeliness']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Lipinski violations -->
      <xsl:choose>
        <xsl:when test="not(fn:number[@key='lipinski_rule_of_5'])">
          <fn:null key="lipinski_violations"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:number key="lipinski_violations">
            <xsl:value-of select="fn:number[@key='lipinski_rule_of_5']"/>
          </fn:number>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Murko framework -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='murko_framework'])">
          <fn:null key="murko_framework"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="murko_framework">
            <xsl:value-of select="fn:string[@key='murko_framework']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ORGANISM SOURCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='organisms']/fn:map">
          <fn:array key="organism_sources">
            <xsl:for-each select="fn:array[@key='organisms']/fn:map">
              <fn:map>
                <fn:string key="name">
                  <xsl:value-of select="fn:string[@key='name']"/>
                </fn:string>
                <xsl:if test="fn:number[@key='taxonomy_id']">
                  <fn:number key="ncbi_taxon_id">
                    <xsl:value-of select="fn:number[@key='taxonomy_id']"/>
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

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">COCONUT</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='coconut_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <!-- ============================================================
       Utility Functions
       ============================================================ -->

  <!-- Check for null/empty values -->
  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
