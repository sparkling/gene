<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  OMIM -> Phenotype Databases Unified Schema Mapping
  ============================================================
  Source: ./schema.md (Online Mendelian Inheritance in Man)
  Target: ../schema.json (3.2 Phenotype Databases Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────────┬──────────────┐
  │ Source Field                │ Target Field            │ Transform    │
  ├─────────────────────────────┼─────────────────────────┼──────────────┤
  │ MIM_Number                  │ mim_number              │ Direct       │
  │ MIM_Number                  │ disease_id              │ Prefix       │
  │ Title                       │ disease_name            │ Direct       │
  │ MIM_Entry_Type              │ entry_type              │ Enum         │
  │ Gene_Symbols                │ gene_symbol             │ Direct       │
  │ Allelic_Variants            │ allelic_variants        │ Array        │
  └─────────────────────────────┴─────────────────────────┴──────────────┘

  NOTES:
  - OMIM uses 6-digit MIM numbers with prefix
  - Entry types: * (gene), # (phenotype), + (combined), % (unknown)
  - Includes allelic variants causing disease
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

      <fn:string key="disease_id">
        <xsl:value-of select="concat('OMIM:', fn:string[@key='MIM_Number'])"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='MIM_Number'])">
          <fn:null key="mim_number"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="mim_number">
            <xsl:value-of select="fn:string[@key='MIM_Number']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Title'])">
          <fn:null key="disease_name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="disease_name">
            <xsl:value-of select="fn:string[@key='Title']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='MIM_Entry_Type'])">
          <fn:null key="entry_type"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="entry_type">
            <xsl:value-of select="fn:string[@key='MIM_Entry_Type']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== GENE INFORMATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='Gene_Symbols'])">
          <fn:null key="gene_symbol"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="gene_symbol">
            <xsl:value-of select="fn:string[@key='Gene_Symbols']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== ALLELIC VARIANTS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='Allelic_Variants']/fn:string">
          <fn:array key="allelic_variants">
            <xsl:for-each select="fn:array[@key='Allelic_Variants']/fn:string">
              <fn:string><xsl:value-of select="."/></fn:string>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="allelic_variants"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">OMIM</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='MIM_Number']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
