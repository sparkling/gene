<?xml version="1.0" encoding="UTF-8"?>
<!--
  ============================================================
  DailyMed -> Pharmaceuticals Unified Schema Mapping
  ============================================================
  Source: ./schema.md (DailyMed FDA Drug Label Database)
  Target: ../schema.json (2.2 Pharmaceuticals Unified Schema)
  XSLT Version: 3.0

  TRANSFORMATION SUMMARY:
  ┌─────────────────────────────┬─────────────────────┬──────────────┐
  │ Source Field                │ Target Field        │ Transform    │
  ├─────────────────────────────┼─────────────────────┼──────────────┤
  │ set_id                      │ set_id              │ Direct       │
  │ product_name                │ name                │ Direct       │
  │ ndc                         │ ndc                 │ Direct       │
  │ spl_section                 │ spl_sections        │ Object array │
  │ unii                        │ unii                │ Direct       │
  │ application_number          │ application_number  │ Direct       │
  │ (derived)                   │ approval_status     │ Set approved │
  └─────────────────────────────┴─────────────────────┴──────────────┘

  NULL HANDLING:
  - product_name: "" → null
  - ndc: "" or "-" → null

  NOTES:
  - DailyMed contains FDA Structured Product Labeling (SPL)
  - set_id is a UUID that uniquely identifies the SPL document
  - NDC is National Drug Code (10-11 digits)
  - UNII is FDA Unique Ingredient Identifier
  - All products in DailyMed are FDA-listed (approved/marketed)
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

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='set_id'])">
          <fn:null key="set_id"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="set_id">
            <xsl:value-of select="fn:string[@key='set_id']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='ndc'])">
          <fn:null key="ndc"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="ndc">
            <xsl:value-of select="fn:string[@key='ndc']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='unii'])">
          <fn:null key="unii"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="unii">
            <xsl:value-of select="fn:string[@key='unii']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='application_number'])">
          <fn:null key="application_number"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="application_number">
            <xsl:value-of select="fn:string[@key='application_number']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== BASIC PROPERTIES ========== -->

      <xsl:choose>
        <xsl:when test="local:is-null(fn:string[@key='product_name'])">
          <fn:null key="name"/>
        </xsl:when>
        <xsl:otherwise>
          <fn:string key="name">
            <xsl:value-of select="fn:string[@key='product_name']"/>
          </fn:string>
        </xsl:otherwise>
      </xsl:choose>

      <!-- All DailyMed products are FDA-listed -->
      <fn:string key="approval_status">approved</fn:string>

      <!-- ========== SPL SECTIONS ========== -->

      <xsl:choose>
        <xsl:when test="fn:array[@key='spl_section']/fn:map">
          <fn:array key="spl_sections">
            <xsl:for-each select="fn:array[@key='spl_section']/fn:map">
              <fn:map>
                <xsl:if test="fn:string[@key='loinc_code']">
                  <fn:string key="loinc_code">
                    <xsl:value-of select="fn:string[@key='loinc_code']"/>
                  </fn:string>
                </xsl:if>
                <xsl:if test="fn:string[@key='title']">
                  <fn:string key="content">
                    <xsl:value-of select="fn:string[@key='title']"/>
                  </fn:string>
                </xsl:if>
              </fn:map>
            </xsl:for-each>
          </fn:array>
        </xsl:when>
        <xsl:otherwise>
          <fn:null key="spl_sections"/>
        </xsl:otherwise>
      </xsl:choose>

      <!-- ========== METADATA ========== -->

      <fn:map key="_source">
        <fn:string key="database">DailyMed</fn:string>
        <fn:string key="version">
          <xsl:value-of select="fn:string[@key='source_version']"/>
        </fn:string>
        <fn:string key="original_id">
          <xsl:value-of select="fn:string[@key='set_id']"/>
        </fn:string>
      </fn:map>

    </fn:map>
  </xsl:template>

  <xsl:function name="local:is-null" as="xs:boolean">
    <xsl:param name="value" as="xs:string?"/>
    <xsl:sequence select="not($value) or $value = '' or $value = '-' or $value = 'NA' or $value = 'N/A' or $value = 'null'"/>
  </xsl:function>

</xsl:stylesheet>
