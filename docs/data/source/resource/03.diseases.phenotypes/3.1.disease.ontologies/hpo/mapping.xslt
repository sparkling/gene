<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  HPO -> Disease Ontologies Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Human Phenotype Ontology)
  Target: ../schema.json (3.1 Disease Ontologies Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ id                          │ id                  │ Direct       │
  │ name                        │ name                │ Direct       │
  │ def                         │ definition          │ Extract      │
  │ synonym                     │ synonyms            │ Array        │
  │ is_a                        │ is_a                │ Array        │
  │ xref                        │ xref                │ Array        │
  │ comment                     │ hpo_comment         │ Direct       │
  │ created_by                  │ created_by          │ Direct       │
  │ creation_date               │ creation_date       │ Direct       │
  │ subset                      │ subset              │ Array        │
  │ alt_id                      │ alt_id              │ Array        │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - HPO uses HP:####### format for IDs
  - Definition may include quotes and citations
  - Used for rare disease diagnosis
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

      <fn:string key="id">
        <xsl:value-of select="fn:string[@key='id']"/>
      </fn:string>

      <!-- ========== BASIC PROPERTIES ========== -->

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

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='def'])">
          <fn:null key="definition"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="definition">
            <xsl:value-of select="fn:string[@key='def']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='comment'])">
          <fn:null key="hpo_comment"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="hpo_comment">
            <xsl:value-of select="fn:string[@key='comment']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CURATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='created_by'])">
          <fn:null key="created_by"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="created_by">
            <xsl:value-of select="fn:string[@key='created_by']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='creation_date'])">
          <fn:null key="creation_date"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="creation_date">
            <xsl:value-of select="fn:string[@key='creation_date']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SYNONYMS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='synonym']/fn:string">
          <fn:array key="synonyms">
            <xsl:for-each select="fn:array[@key='synonym']/fn:string">
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
        <xsl:when test="fn:array[@key='is_a']/fn:string">
          <fn:array key="is_a">
            <xsl:for-each select="fn:array[@key='is_a']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="is_a"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== SUBSETS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='subset']/fn:string">
          <fn:array key="subset">
            <xsl:for-each select="fn:array[@key='subset']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="subset"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ALTERNATIVE IDS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='alt_id']/fn:string">
          <fn:array key="alt_id">
            <xsl:for-each select="fn:array[@key='alt_id']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="alt_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== CROSS-REFERENCES ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='xref']/fn:string">
          <fn:array key="xref">
            <xsl:for-each select="fn:array[@key='xref']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="xref"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">HPO</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
