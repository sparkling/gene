<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  Allen Brain Atlas -> Mental Health/Neurological Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Allen Brain Atlas)
  Target: ../schema.json (3.7 Mental Health/Neurological Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ gene_symbol                 │ gene_symbol             │ Direct       │
  │ structure_id                │ structure_id            │ Direct       │
  │ structure_acronym           │ structure_acronym       │ Direct       │
  │ structure_name              │ brain_regions           │ Array        │
  │ donor_id                    │ donor_id                │ Direct       │
  │ expression_level            │ expression_level        │ Direct       │
  │ data_type                   │ data_type               │ Enum         │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - Allen Brain Atlas provides genome-wide expression across ~900 brain structures
  - Human brain samples from 6 donors
  - Data types: microarray, RNA-seq, ISH
  - Structure IDs reference Allen Brain Ontology
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
      <!-- ========== GENE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='gene_symbol'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='gene_symbol']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== STRUCTURE IDENTIFIERS ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='structure_id']">
          <fn:number key="structure_id">
            <xsl:value-of select="fn:number[@key='structure_id']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="structure_id"/>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='structure_acronym'])">
          <fn:null key="structure_acronym"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="structure_acronym">
            <xsl:value-of select="fn:string[@key='structure_acronym']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BRAIN REGIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='structure_name']/fn:string">
          <fn:array key="brain_regions">
            <xsl:for-each select="fn:array[@key='structure_name']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:when test="fn:string[@key='structure_name'] and not(local:is-null(fn:string[@key='structure_name']))">
          <fn:array key="brain_regions">
            <fn:string><xsl:value-of select="fn:string[@key='structure_name']"/></fn:string>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="brain_regions"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DONOR ID ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='donor_id'])">
          <fn:null key="donor_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="donor_id">
            <xsl:value-of select="fn:string[@key='donor_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== EXPRESSION LEVEL ========== -->

      <xsl:choose>
        <xsl:when test="fn:number[@key='expression_level']">
          <fn:number key="expression_level">
            <xsl:value-of select="fn:number[@key='expression_level']"/>
          </fn:number>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="expression_level"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DATA TYPE ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='data_type'])">
          <fn:null key="data_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="data_type">
            <xsl:value-of select="fn:string[@key='data_type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">Allen Brain Atlas</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="concat(fn:string[@key='gene_symbol'], '_', fn:number[@key='structure_id'])"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
