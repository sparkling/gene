<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  RxNorm -> Pharmaceuticals Unified Schema Mapping
  ============================================================
  Source: ./schema.md (RxNorm - NLM Drug Vocabulary)
  Target: ../schema.json (2.2 Pharmaceuticals Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ rxcui                       │ drug_id             │ Direct       │
  │ rxcui                       │ rxcui               │ Direct       │
  │ name                        │ name                │ Direct       │
  │ tty                         │ tty                 │ Direct       │
  │ sab                         │ sab                 │ Direct       │
  │ rela                        │ rela                │ Direct       │
  │ dose_form                   │ dose_form           │ Direct       │
  │ strength                    │ strength            │ Direct       │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - name: "" → null
  - All fields: empty values → null

  NOTES:
  - RxNorm uses RXCUI as primary identifier
  - TTY = Term Type (IN, SCD, SBD, GPCK, BPCK, etc.)
  - SAB = Source Abbreviation
  - RELA = Relationship Attribute
  - RxNorm is the standard vocabulary for clinical drugs in US
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

      <fn:string key="drug_id">
        <xsl:value-of select="fn:string[@key='rxcui']"/>
      </fn:string>

      <fn:string key="rxcui">
        <xsl:value-of select="fn:string[@key='rxcui']"/>
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

      <!-- ========== RXNORM ATTRIBUTES ========== -->

      <!-- Term Type -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='tty'])">
          <fn:null key="tty"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="tty">
            <xsl:value-of select="fn:string[@key='tty']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Source Abbreviation -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='sab'])">
          <fn:null key="sab"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="sab">
            <xsl:value-of select="fn:string[@key='sab']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- Relationship Attribute -->
      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='rela'])">
          <fn:null key="rela"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="rela">
            <xsl:value-of select="fn:string[@key='rela']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== DOSAGE INFORMATION ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='dose_form'])">
          <fn:null key="dose_form"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="dose_form">
            <xsl:value-of select="fn:string[@key='dose_form']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='strength'])">
          <fn:null key="strength"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="strength">
            <xsl:value-of select="fn:string[@key='strength']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">RxNorm</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='rxcui']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
