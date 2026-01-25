<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ChEBI -> Chemical Ontology Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Chemical Entities of Biological Interest)
  Target: ../schema.json (2.6 Chemical Ontology and Classification Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ CHEBI_ID                    │ chebi_id            │ Direct       │
  │ CHEBI_ID                    │ compound_id         │ Direct       │
  │ NAME                        │ name                │ Direct       │
  │ SMILES                      │ canonical_smiles    │ Direct       │
  │ InChI                       │ inchi               │ Direct       │
  │ InChIKey                    │ inchi_key           │ Direct       │
  │ FORMULA                     │ molecular_formula   │ Direct       │
  │ DEFINITION                  │ definition          │ Direct       │
  │ SYNONYM                     │ synonyms            │ Array        │
  │ MASS                        │ mass                │ Direct       │
  │ CHARGE                      │ charge              │ Direct       │
  │ STAR                        │ star                │ Direct       │
  │ IS_A                        │ is_a                │ Array        │
  │ HAS_ROLE                    │ has_role            │ Array        │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - ChEBI uses CHEBI: prefix for IDs (e.g., CHEBI:17234)
  - Star rating: 1=automatic, 2=basic, 3=expert-curated
  - Includes ontology relationships (is_a, has_role)
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

  <xsl:template name="xsl:initial-template">
    <xsl:param name="input" as="xs:string"/>
    <xsl:variable name="source" select="json-to-xml($input)"/>
    <xsl:apply-templates select="$source/fn:map"/>
  </xsl:template>

  <xsl:template match="fn:map">
    <fn:map>
      <!-- ========== IDENTIFIERS ========== -->

      <fn:string key="chebi_id">
        <xsl:value-of select="fn:string[@key='CHEBI_ID']"/>
      </fn:string>

      <fn:string key="compound_id">
        <xsl:value-of select="fn:string[@key='CHEBI_ID']"/>
      </fn:string>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='NAME'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='NAME']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='DEFINITION'])">
          <fn:null key="definition"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="definition">
            <xsl:value-of select="fn:string[@key='DEFINITION']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CHEMICAL STRUCTURE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='SMILES'])">
          <fn:null key="canonical_smiles"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="canonical_smiles">
            <xsl:value-of select="fn:string[@key='SMILES']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='InChI'])">
          <fn:null key="inchi"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi">
            <xsl:value-of select="fn:string[@key='InChI']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='InChIKey'])">
          <fn:null key="inchi_key"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="inchi_key">
            <xsl:value-of select="fn:string[@key='InChIKey']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='FORMULA'])">
          <fn:null key="molecular_formula"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="molecular_formula">
            <xsl:value-of select="fn:string[@key='FORMULA']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== PHYSICAL PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='MASS']">
          <fn:number key="mass">
            <xsl:value-of select="fn:number[@key='MASS']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="mass"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:number[@key='CHARGE']">
          <fn:number key="charge">
            <xsl:value-of select="fn:number[@key='CHARGE']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="charge"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CURATION ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='STAR']">
          <fn:number key="star">
            <xsl:value-of select="fn:number[@key='STAR']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="star"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SYNONYMS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='SYNONYM']/fn:string">
          <fn:array key="synonyms">
            <xsl:for-each select="fn:array[@key='SYNONYM']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="synonyms"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ONTOLOGY RELATIONSHIPS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='IS_A']/fn:string">
          <fn:array key="is_a">
            <xsl:for-each select="fn:array[@key='IS_A']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="fn:array[@key='HAS_ROLE']/fn:string">
          <fn:array key="has_role">
            <xsl:for-each select="fn:array[@key='HAS_ROLE']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="has_role"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='IS_CONJUGATE_ACID_OF'])">
          <fn:null key="is_conjugate_acid_of"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="is_conjugate_acid_of">
            <xsl:value-of select="fn:string[@key='IS_CONJUGATE_ACID_OF']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='IS_ENANTIOMER_OF'])">
          <fn:null key="is_enantiomer_of"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="is_enantiomer_of">
            <xsl:value-of select="fn:string[@key='IS_ENANTIOMER_OF']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">ChEBI</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='CHEBI_ID']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
