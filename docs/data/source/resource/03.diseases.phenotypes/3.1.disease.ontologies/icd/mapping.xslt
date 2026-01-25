<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  ICD -> Disease Ontologies Unified Schema Mapping
  ============================================================
  Source: ./schema.md (International Classification of Diseases)
  Target: ../schema.json (3.1 Disease Ontologies Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ code                        │ id                  │ Direct       │
  │ title                       │ name                │ Direct       │
  │ icd10_code                  │ icd10_code          │ Direct       │
  │ icd11_code                  │ icd11_code          │ Direct       │
  │ block_code                  │ block_code          │ Direct       │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - All string fields: "" → null

  NOTES:
  - ICD-10 codes: A00-Z99 format
  - ICD-11 codes: 4-6 character alphanumeric
  - WHO standard for mortality and morbidity statistics
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
        <xsl:value-of select="fn:string[@key='code']"/>
      </fn:string>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='icd10_code'])">
          <fn:null key="icd10_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="icd10_code">
            <xsl:value-of select="fn:string[@key='icd10_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='icd11_code'])">
          <fn:null key="icd11_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="icd11_code">
            <xsl:value-of select="fn:string[@key='icd11_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='block_code'])">
          <fn:null key="block_code"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="block_code">
            <xsl:value-of select="fn:string[@key='block_code']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='title'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='title']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">ICD</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='code']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
